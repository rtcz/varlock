import hashlib
import random

# import common
# import varlock.iterator as iters
#
# from .cigar import Cigar
# from .diff import Diff
# from .fasta_index import FastaIndex
# from .po import SnvAlignment, GenomicPosition, VacSnvRecord


from varlock.cigar import Cigar
import varlock.common as common
from varlock.diff import Diff
from varlock.fasta_index import FastaIndex
import varlock.iterator as iters
from varlock.po import SnvAlignment, GenomicPosition, VacSnvRecord


class Mutator:
    def __init__(
            self,
            bam_file,
            rnd=random.SystemRandom(),
            verbose: bool = False
    ):
        """
        :param rnd: (secure) random generator
        :param verbose:
        """
        self.rnd = rnd
        self.verbose = verbose
        self.fai = FastaIndex(bam_file)
        self.bam_file = bam_file
        self._bam_checksum = None
        
        self.prev_alignment = None
    
    def __is_before_index(self, alignment, index):
        """
        Check if alignment is mapped before index.
        :return: True if alignment end is mapped before index
        """
        # reference_end points to one past the last aligned residue
        ref_end = alignment.reference_end - 1
        alignment_end = self.fai.pos2index(alignment.reference_name, ref_end)
        return alignment_end < index
    
    def __is_after_index(self, alignment, index):
        """
        Check if alignment is mapped after index.
        :return: True if alignment start is mapped after index
        """
        alignment_start = self.fai.pos2index(alignment.reference_name, alignment.reference_start)
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
        
        if self.verbose and self.alignment_counter % 10000 == 0:
            print("%d alignments processed" % self.alignment_counter)
    
    def __init_counters(self):
        self.alignment_counter = 0  # written alignments
        self.diff_counter = 0  # processed diff records
        self.mut_counter = 0  # all mutations
        self.vac_counter = 0  # read vac records (variants)
        
        self.mut_alignment_counter = 0  # mutated alignments
        self.max_coverage = 0  # maximum alignments overlapping single SNV
    
    def resolve_diff_range(
            self,
            diff_file,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
    ):
        checksum, diff_start_index, diff_end_index, secret = Diff.read_header(diff_file)
        self.validate_checksum(checksum)
        
        return self.resolve_range(
            diff_start_index,
            diff_end_index,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
        )
    
    def resolve_range(
            self,
            diff_start_index,
            diff_end_index,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
    ):
        """
        Resolve between DIFF range and user specified range.
        :param diff_start_index:
        :param diff_end_index:
        :param start_ref_name:
        :param start_ref_pos:
        :param end_ref_name:
        :param end_ref_pos:
        :return: tuple (start_index, end_index)
        """
        if start_ref_name is not None or end_ref_name is not None:
            # apply user range
            start_index = self.fai.resolve_start_index(start_ref_name, start_ref_pos)
            end_index = self.fai.resolve_end_index(end_ref_name, end_ref_pos)
            
            diff_range = self.fai.index2pos(diff_start_index) + self.fai.index2pos(diff_end_index)  # type: tuple
            
            if start_index < diff_start_index:
                args = self.fai.index2pos(start_index) + diff_range  # type: tuple
                # noinspection PyStringFormat
                raise ValueError("Start position [%s, %d] must be within DIFF range [%s, %d]:[%s, %d]." % args)
            if end_index > diff_end_index:
                args = self.fai.index2pos(end_index) + diff_range  # type: tuple
                # noinspection PyStringFormat
                raise ValueError("End position [%s, %d] must be within DIFF range [%s, %d]:[%s, %d]." % args)
            return start_index, end_index
        
        else:
            # apply diff range
            return diff_start_index, diff_end_index
    
    def bam_checksum(self):
        if self._bam_checksum is None:
            if self.verbose:
                print("Calculating BAM's checksum")
            
            self._bam_checksum = common.filename_checksum(self.bam_file.filename)
        
        return self._bam_checksum
    
    def validate_checksum(self, check_sum):
        if self.bam_checksum() != check_sum:
            # print(bin2hex(self.bam_checksum()))
            # print(bin2hex(check_sum))
            raise ValueError("Invalid checksum.")
    
    # TODO refactor iterators as parameters
    def unmutate(
            self,
            diff_file,
            out_bam_file,
            start_ref_name: str = None,
            start_ref_pos: int = None,
            end_ref_name: str = None,
            end_ref_pos: int = None,
            include_unmapped: bool = False,
            unmapped_only: bool = False
    ):
        """
        Unmutate BAM file in range specified by DIFF file or by parameters.
        :param diff_file: binary file
        :param out_bam_file: pysam.AlignmentFile
        :param start_ref_name: inclusive
        :param start_ref_pos: 0-based, inclusive
        :param end_ref_name: inclusive
        :param end_ref_pos: 0-based, inclusive
        :param include_unmapped: Include all unplaced unmapped reads.
        :param unmapped_only: Only unmapped reads - both placed and unplaced.
         Overrides other parameters.
         
         When range is supplied partialy covered reads are also included,
         but only snv's within range are unmutated.
        """
        self.__init_counters()
        
        Diff.validate(diff_file)
        start_index, end_index = self.resolve_diff_range(
            diff_file,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
        )
        secret = Diff.read_header(diff_file)[3]
        
        alignment_queue = []
        
        if unmapped_only:
            bam_iter = iters.UnmappedBamIterator(self.bam_file)
        elif include_unmapped:
            bam_iter = iters.RangedBamIterator(self.bam_file, start_index, end_index)
        else:
            # mapped only
            bam_iter = iters.MappedBamIterator(self.bam_file, start_index, end_index)
        
        alignment = next(bam_iter)
        diff_iter = iters.DiffIterator(diff_file, self.fai, start_index, end_index)
        diff = next(diff_iter)
        
        if self.verbose:
            print('first diff: %s' % diff)
            print('first alignment: %s' % self.alignment2str(alignment))
        
        while True:
            if diff is None or alignment is None:
                # finish
                if self.verbose:
                    if diff is None:
                        print("EOF DIFF")
                    else:
                        print("EOF BAM")
                
                # last mutation
                if diff is not None:
                    # noinspection PyTypeChecker
                    self.__unmutate_overlap(alignment_queue, diff.mut_map)
                
                for snv_alignment in alignment_queue:
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(out_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF DIFF case only)
                while alignment is not None:
                    self.__encrypt_unmapped(alignment, secret)
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = next(bam_iter)
                
                if self.verbose:
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
                    alignment_queue.append(SnvAlignment(alignment, None))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, diff.index):
                self.__unmutate_overlap(alignment_queue, diff.mut_map)
                # done with this diff, read next
                prev_diff = diff
                diff = next(diff_iter)
                if diff is not None:
                    # could be end of DIFF file
                    self.__write_done(out_bam_file, alignment_queue, diff.index)
                    self.__set_seq_positions(alignment_queue, diff.ref_pos)
                elif self.verbose:
                    print('last diff: %s' % prev_diff)
            
            else:  # alignment is covering diff position
                # find sequence position of vac
                seq_pos = common.ref_pos2seq_pos(alignment, diff.ref_pos)
                alignment_queue.append(SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.mut_alignment_counter += 1
                # noinspection PyAttributeOutsideInit
                # TODO this is not the real coverage
                self.max_coverage = max(len(alignment_queue), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.diff_counter = diff_iter.counter
    
    def __unmutate_overlap(self, alignment_queue, mut_map):
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
            return '#%d %s:%d' % (self.fai.pos2index(ref_name, ref_start), ref_name, ref_start)
        else:
            return alignment.query_sequence
    
    # TODO refactor, iterators as parameters
    def mutate(self, vac_filename: str, out_bam_file, out_diff_file, secret: bytes):
        """
        Mutate BAM file SNVs.
        :param vac_filename:
        input variant allele count file
        :param out_bam_file: pysam.AlignmentFile
        output bam file
        :param out_diff_file: binary file
        output diff file
        :param secret: Secret key written into DIFF used for unmapped alignment encryption.
        """
        self.__init_counters()
        
        alignment_queue = []
        
        bam_iter = iters.FullBamIterator(self.bam_file)
        alignment = next(bam_iter)
        
        vac_iter = iters.VacIterator(vac_filename, self.fai)
        vac = next(vac_iter)
        
        if self.verbose:
            print('first vac: %s' % vac)
            print('first alignment: %s' % self.alignment2str(alignment))
        
        Diff.write_header(
            diff_file=out_diff_file,
            start_index=self.fai.first_index(),
            end_index=self.fai.last_index(),
            secret=secret
        )
        while True:
            if vac is None or alignment is None:
                # finish
                if self.verbose:
                    if vac is None:
                        print("EOF VAC")
                    else:
                        print("EOF BAM")
                
                # last mutation
                if vac is not None:
                    # noinspection PyTypeChecker
                    self.__mutate_overlap(out_diff_file, alignment_queue, vac)
                
                for snv_alignment in alignment_queue:
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(out_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF VAC case only)
                while alignment is not None:
                    self.__encrypt_unmapped(alignment, secret)
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = next(bam_iter)
                
                if self.verbose:
                    print('last alignment: %s' % self.alignment2str(self.prev_alignment))
                    print("total of %d alignments processed" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, vac.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceeding alignments are written
                    self.__write_alignment(out_bam_file, alignment)
                else:
                    # append to queue
                    alignment_queue.append(SnvAlignment(alignment, None))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, vac.index):
                self.__mutate_overlap(out_diff_file, alignment_queue, vac)
                # done with this vac, read next
                prev_vac = vac
                vac = next(vac_iter)
                if vac is not None:
                    # could be end of VAC file
                    self.__write_done(out_bam_file, alignment_queue, vac.index)
                    self.__set_seq_positions(alignment_queue, vac.ref_pos)
                elif self.verbose:
                    print('last vac: %s' % prev_vac)
            
            else:  # alignment is covering vac position
                # find sequence position of vac
                seq_pos = common.ref_pos2seq_pos(alignment, vac.ref_pos)
                alignment_queue.append(SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.mut_alignment_counter += 1
                # noinspection PyAttributeOutsideInit
                # TODO this is not the real coverage
                self.max_coverage = max(len(alignment_queue), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.vac_counter = vac_iter.counter
    
    @staticmethod
    def __encrypt_unmapped(alignment, secret):
        """
        Stream cipher encryption method both for encryption and decryption.
        alignment + secret => encrypted_alignment
        encrypted_alignment + secret => alignment
        :param alignment:
        :param secret:
        :return: encrypter/decrypted alignment
        """
        # TODO query_qualities ?
        if alignment.is_unmapped:
            # use 64B long hash (encrypts 256 bases)
            sha512 = hashlib.sha512()
            sha512.update(secret + alignment.query_name.encode())
            mut_seq = common.stream_cipher(alignment.query_sequence, sha512.digest())
            alignment.query_sequence = mut_seq
    
    def __mutate_overlap(self, out_diff_file, alignment_queue, vac: GenomicPosition):
        """
        Mutate alignments at position of the variant
        :param alignment_queue: alignments together with vac positions
        :param vac: variant allele count list
        :return: True if NS mutation has occured
        """
        is_mutated = False
        # get pileup of bases from alignments at vac mapping position
        base_pileup = common.get_base_pileup(alignment_queue)
        alt_ac = common.count_bases(base_pileup)
        
        if isinstance(vac, VacSnvRecord):
            mut_map = common.create_mut_map(alt_ac=alt_ac, ref_ac=vac.ac, rnd=self.rnd)
            
            for snv_alignment in alignment_queue:
                if snv_alignment.pos is not None:
                    # alignment has vac to mutate
                    is_mutated |= self.__mutate_alignment(snv_alignment.alignment, snv_alignment.pos, mut_map)
                    # done with current vac
                    snv_alignment.pos = None
            
            if is_mutated:
                # at least one alignment has been mutated
                self.__write_diff(
                    out_diff_file=out_diff_file,
                    index=vac.index,
                    mut_map=mut_map
                )
        
        # TODO isinstance(vac, VacIndelRecord)
        
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
    
    def __write_diff(self, out_diff_file, index, mut_map):
        self.diff_counter += 1
        mut_tuple = (mut_map['A'], mut_map['T'], mut_map['G'], mut_map['C'])
        Diff.write_record(out_diff_file, index, mut_tuple)
    
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
