import hashlib
import random
import pysam

import varlock.bdiff as bdiff
import varlock.common as common
import varlock.iters as iters
import varlock.po as po

from varlock.cigar import Cigar
from varlock.fasta_index import FastaIndex


class Mutator:
    def __init__(
            self,
            fai: FastaIndex,
            rnd=random.SystemRandom(),
            verbose: bool = False
    ):
        """
        :param rnd: (secure) random generator
        :param verbose:
        """
        self._fai = fai
        self._rnd = rnd
        self._verbose = verbose
        self._bam_checksum = None
        
        self.prev_alignment = None
    
    def __is_before_index(self, alignment, index):
        """
        Check if alignment is mapped before index.
        :return: True if alignment end is mapped before index
        """
        # reference_end points to one past the last aligned residue
        ref_end = alignment.reference_end - 1
        alignment_end = self._fai.pos2index(alignment.reference_name, ref_end)
        return alignment_end < index
    
    def __is_after_index(self, alignment, index):
        """
        Check if alignment is mapped after index.
        :return: True if alignment start is mapped after index
        """
        alignment_start = self._fai.pos2index(alignment.reference_name, alignment.reference_start)
        return alignment_start > index
    
    def __write_alignment(self, bam_file, alignment):
        """
        :param alignment: pysam.AlignedSegment
        :return:
        """
        are_placed = self.prev_alignment is not None and common.is_placed_alignment(
            alignment) and common.is_placed_alignment(self.prev_alignment)
        # FIXME is case of placed after unplaced alignment valid ?
        if are_placed and alignment.reference_start < self.prev_alignment.reference_start:
            # safety check
            raise IndexError('Unordered write: %s after %s' % (
                self.alignment2str(alignment),
                self.alignment2str(self.prev_alignment)
            ))
        
        bam_file.write(alignment)
        self.alignment_counter += 1
        self.prev_alignment = alignment
        
        if self._verbose and self.alignment_counter % 10000 == 0:
            print("%d alignments processed" % self.alignment_counter)
    
    def __init_counters(self):
        self.alignment_counter = 0  # written alignments
        self.diff_counter = 0  # processed diff records
        self.mut_counter = 0  # all mutations
        self.vac_counter = 0  # read vac records (variants)
        
        self.covering_counter = 0  # mutated alignments
        self.max_coverage = 0  # maximum alignments overlapping single SNV
    
    def mutate(
            self,
            mut_bam_file: pysam.AlignmentFile,
            vac_iter: iters.VacIterator,
            bam_iter: iters.FullBamIterator,
            secret: bytes
    ):
        """
        :param mut_bam_file:
        :param vac_iter:
        :param bam_iter:
        :param secret:
        :return: BdiffIO object
        BDIFF records are written only if bases at VAC position differ after mutation.
        """
        self.__init_counters()
        
        alignment_queue = []
        alignment = next(bam_iter)
        vac = next(vac_iter)
        
        if self._verbose:
            print('first vac: %s' % vac)
            print('first alignment: %s' % self.alignment2str(alignment))
        
        bdiff_io = bdiff.BdiffIO()
        while True:
            if vac is None or alignment is None:
                # finish
                if self._verbose:
                    if vac is None:
                        print("EOF VAC")
                    else:
                        print("EOF BAM")
                
                # last mutation
                if vac is not None:
                    # noinspection PyTypeChecker
                    self.__mutate_pos(bdiff_io, alignment_queue, vac)
                
                for snv_alignment in alignment_queue:
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(mut_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF VAC case only)
                while alignment is not None:
                    self.__encrypt_unmapped(alignment, secret)
                    self.__write_alignment(mut_bam_file, alignment)
                    alignment = next(bam_iter)
                
                if self._verbose:
                    print('last alignment: %s' % self.alignment2str(self.prev_alignment))
                    print("total of %d alignments processed" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, vac.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceeding alignments are written
                    self.__write_alignment(mut_bam_file, alignment)
                else:
                    # append to queue
                    alignment_queue.append(po.SnvAlignment(alignment, None))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, vac.index):
                self.__mutate_pos(bdiff_io, alignment_queue, vac)
                # done with this vac, read next
                prev_vac = vac
                vac = next(vac_iter)
                if vac is not None:
                    # could be end of VAC file
                    self.__write_done(mut_bam_file, alignment_queue, vac.index)
                    self.__set_seq_positions(alignment_queue, vac.ref_pos)
                elif self._verbose:
                    print('last vac: %s' % prev_vac)
            
            else:  # alignment is covering vac position
                # find sequence position of vac
                seq_pos = common.ref_pos2seq_pos(alignment, vac.ref_pos)
                alignment_queue.append(po.SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.covering_counter += 1
                # noinspection PyAttributeOutsideInit
                # TODO this is not the real coverage
                self.max_coverage = max(len(alignment_queue), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.vac_counter = vac_iter.counter
        # noinspection PyAttributeOutsideInit
        self.diff_counter = bdiff_io.snv_count + bdiff_io.indel_count
        return bdiff_io
    
    def unmutate(
            self,
            bam_iter: iters.BamIterator,
            bdiff_iter: iters.BdiffIterator,
            out_bam_file,
            secret: bytes
    ):
        self.__init_counters()
        
        alignment_queue = []
        alignment = next(bam_iter)
        diff = next(bdiff_iter)
        if self._verbose:
            print('first diff: %s' % diff)
            print('first alignment: %s' % self.alignment2str(alignment))
        
        while True:
            if diff is None or alignment is None:
                # finish
                if self._verbose:
                    if diff is None:
                        print("EOF DIFF")
                    else:
                        print("EOF BAM")
                
                # last mutation
                if diff is not None:
                    # noinspection PyTypeChecker
                    self.__unmutate_pos(alignment_queue, diff.mut_map)
                
                for snv_alignment in alignment_queue:
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(out_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF DIFF case only)
                while alignment is not None:
                    self.__encrypt_unmapped(alignment, secret)
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = next(bam_iter)
                
                if self._verbose:
                    print('last alignment: %s' % self.alignment2str(self.prev_alignment))
                    print("total of %d alignments processed" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, diff.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceeding alignments are written
                    self.__write_alignment(out_bam_file, alignment)
                else:
                    # append to queue
                    alignment_queue.append(po.SnvAlignment(alignment, None))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, diff.index):
                self.__unmutate_pos(alignment_queue, diff.mut_map)
                # done with this diff, read next
                prev_diff = diff
                diff = next(bdiff_iter)
                if diff is not None:
                    # could be end of DIFF file
                    self.__write_done(out_bam_file, alignment_queue, diff.index)
                    self.__set_seq_positions(alignment_queue, diff.ref_pos)
                elif self._verbose:
                    print('last diff: %s' % prev_diff)
            
            else:  # alignment is covering diff position
                # find sequence position of vac
                seq_pos = common.ref_pos2seq_pos(alignment, diff.ref_pos)
                alignment_queue.append(po.SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.covering_counter += 1
                # noinspection PyAttributeOutsideInit
                # TODO this is not the real coverage
                self.max_coverage = max(len(alignment_queue), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.diff_counter = bdiff_iter.counter
    
    def __unmutate_pos(self, alignment_queue, mut_map):
        for snv_alignment in alignment_queue:
            if snv_alignment.pos is not None:
                # alignment has vac to mutate
                self.__mutate_alignment(snv_alignment.alignment, snv_alignment.pos, mut_map)
                # done with current vac
                snv_alignment.pos = None
    
    def alignment2str(self, alignment):
        if common.is_placed_alignment(alignment):
            ref_name = alignment.reference_name
            ref_start = alignment.reference_start
            # unplaced alignment
            return '#%d %s:%d' % (self._fai.pos2index(ref_name, ref_start), ref_name, ref_start)
        else:
            return alignment.query_sequence
    
    @staticmethod
    def __encrypt_unmapped(alignment, secret: bytes):
        """
        Stream cipher encryption / decryption.
        alignment + secret => encrypted_alignment
        encrypted_alignment + secret => alignment
        :param alignment:
        :param secret:
        :return: encrypter/decrypted alignment
        """
        if alignment.is_unmapped:
            # use 64B long hash (encrypts 256 bases)
            sha512 = hashlib.sha512()
            # TODO maybe include query quality?
            sha512.update(secret + alignment.query_name.encode())
            
            mut_seq = common.stream_cipher(alignment.query_sequence, sha512.digest())
            alignment.query_sequence = mut_seq
    
    def __mutate_pos(self, bdiff_io: bdiff.BdiffIO, alignment_queue, vac: po.GenomicPosition):
        """
        Mutate alignments at the variant allele count position.
        :param alignment_queue: alignments together with vac positions
        :param vac: variant allele count list
        :return: True if NS mutation has occured
        """
        is_mutated = False
        # get pileup of bases from alignments at vac mapping position
        base_pileup = common.get_base_pileup(alignment_queue)
        alt_ac = common.count_bases(base_pileup)
        
        if isinstance(vac, po.VacSnvRecord):
            mut_map = common.create_mut_map(alt_ac=alt_ac, ref_ac=vac.ac, rnd=self._rnd)
            
            for snv_alignment in alignment_queue:
                if snv_alignment.pos is not None:
                    # alignment has vac to mutate
                    is_mutated |= self.__mutate_alignment(snv_alignment.alignment, snv_alignment.pos, mut_map)
                    # done with current vac
                    snv_alignment.pos = None
            
            if is_mutated:
                # at least one alignment has been mutated
                self.__write_diff(
                    bdiff_io=bdiff_io,
                    index=vac.index,
                    mut_map=mut_map
                )
        
        elif isinstance(vac, po.VacIndelRecord):
            # TODO isinstance(vac, VacIndelRecord)
            raise ValueError("VAC indel not implemented yet!")
        else:
            raise ValueError("%s is not VAC record instance" % type(vac).__name__)
        
        return is_mutated
    
    def __mutate_alignment(self, alignment, seq_pos, mut_map):
        """
        Mutate alignment by mutation map at SNV position.
        :param alignment: pysam.AlignedSegment
        :param seq_pos: position of SNV in aligned sequence
        :param mut_map: mutation map
        :return: True if base has been mutated
        """
        is_mutated = False
        # alignment is mapped at snv position
        snv_base = common.get_base(alignment, seq_pos)
        mut_base = mut_map[snv_base]
        
        if snv_base != mut_base:
            # base has been mutated to another base
            self.mut_counter += 1
            is_mutated = True
            common.set_base(alignment, seq_pos, mut_base)
            Cigar.validate(alignment, seq_pos)
        
        return is_mutated
    
    # noinspection PyMethodMayBeStatic
    def __write_diff(self, bdiff_io: bdiff.BdiffIO, index: int, mut_map: dict):
        """
        :param bdiff_io:
        :param index:
        :param mut_map:
        mut_map.key: original base
        mut_map.value: mutated base
        """
        bdiff_io.write_snv(index, (mut_map['A'], mut_map['T'], mut_map['G'], mut_map['C']))
    
    def __write_done(self, out_bam_file, alignment_queue, index):
        """
        Write snv alignments while their end reference position is before index.
        :param out_bam_file:
        :param alignment_queue:
        :param index:
        :return:
        """
        tmp_alignment_queue = alignment_queue[:]
        for snv_alignment in tmp_alignment_queue:
            is_mapped = not snv_alignment.alignment.is_unmapped
            if is_mapped and not self.__is_before_index(snv_alignment.alignment, index):
                # break to keep alignments in order
                # otherwise shorter alignment could be written before longer alignment
                break
            
            self.__write_alignment(out_bam_file, snv_alignment.alignment)
            # remove written alignment
            alignment_queue.remove(snv_alignment)
    
    @staticmethod
    def __set_seq_positions(alignment_queue, ref_pos):
        """
        Set SNV position derived from reference position for each alignment.
        :param alignment_queue:
        :param ref_pos:
        :return:
        """
        for snv_alignment in alignment_queue:
            # TODO optimize, skip when ref_pos is greater than alignment reference end
            snv_pos = common.ref_pos2seq_pos(snv_alignment.alignment, ref_pos)
            if snv_pos is not None:
                # alignment is mapped at another new_snv position
                snv_alignment.pos = snv_pos
