import hashlib

import pysam

import varlock_src.bdiff as bdiff
import varlock_src.common as cmn
import varlock_src.iters as iters
import varlock_src.po as po
from varlock_src.fasta_index import FastaIndex
from varlock_src.random import VeryRandom
from varlock_src.variant import AlignedVariant


class Mutator:
    def __init__(
            self,
            fai: FastaIndex,
            verbose: bool = False
    ):
        """
        :param verbose:
        """
        self._fai = fai
        self._verbose = verbose
        self._prev_alignment = None
    
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
        are_placed = self._prev_alignment is not None and cmn.is_placed_alignment(alignment) and cmn.is_placed_alignment(self._prev_alignment)
        # FIXME is case of placed after unplaced alignment valid ?
        if are_placed:
            # get the reference or pass without checking
            try:
                ref_name = alignment.reference_name
                ref_name_prev = self._prev_alignment.reference_name
            except ValueError:
                ref_name = None
                ref_name_prev = None
            # unplaced alignment / wrong order
            if ref_name is not None and ref_name_prev is not None and self._fai.pos2index(ref_name, alignment.reference_start) < self._fai.pos2index(ref_name_prev, self._prev_alignment.reference_start):
                # safety check
                raise IndexError('Unordered write: %s after %s' % (
                    self.alignment2str(alignment),
                    self.alignment2str(self._prev_alignment)
                ))
        
        bam_file.write(alignment)
        self.alignment_counter += 1
        self._prev_alignment = alignment
        
        if self._verbose and self.alignment_counter % 1000 == 0:
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
            vac_iter: iters.VariantIterator,
            bam_iter: iters.FullBamIterator,
            secret: bytes,
            rnd: VeryRandom
    ) -> bdiff.BdiffIO:
        """
        :param mut_bam_file:
        :param vac_iter:
        :param bam_iter:
        :param secret:
        :param rnd: (secure) random generator
        :return: BdiffIO object
        BDIFF records are written only if bases at VAC position differ after mutation.
        """
        self.__init_counters()
        
        alignment_queue = []
        alignment = next(bam_iter)

        # TODO optimization: seek to first alignment covering vac to skip all preceding records
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
                    self.__mutate_pos(bdiff_io, alignment_queue, vac, rnd)
                
                for snv_alignment in alignment_queue:
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(mut_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF VAC case only)
                while alignment is not None:
                    self.__encrypt_unmapped(alignment, secret)
                    self.__write_alignment(mut_bam_file, alignment)
                    alignment = next(bam_iter)
                
                if self._verbose:
                    print("total of %d alignments processed" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, vac.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceding alignments are written
                    self.__write_alignment(mut_bam_file, alignment)
                else:
                    # append to queue
                    alignment_queue.append(AlignedVariant(alignment))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, vac.index):
                # apply the variant to alignment
                self.__mutate_pos(bdiff_io, alignment_queue, vac, rnd)
                # done with this vac, read next
                prev_vac = vac
                vac = next(vac_iter)
                if self._verbose and vac_iter.counter % 100000 == 0:
                    print('%d VAC records done' % vac_iter.counter)
                # not the end of VAC file
                if vac is not None:
                    # write alignments to mutated BAM
                    self.__write_before_index(mut_bam_file, alignment_queue, vac.index)
                    # update alignment queue
                    for i in range(len(alignment_queue)):
                        if not alignment_queue[i].alignment.is_unmapped:
                            alignment_variant = cmn.vac_aligned_variant(alignment_queue[i].alignment, vac, True)
                            if alignment_variant is not None:
                                alignment_queue[i] = alignment_variant
                
                elif self._verbose:
                    print('last vac: %s' % prev_vac)
            
            else:  # alignment is covering vac position
                alignment_variant = cmn.vac_aligned_variant(alignment, vac, True)
                if alignment_variant is not None:
                    alignment_queue.append(alignment_variant)
                alignment = next(bam_iter)
                self.covering_counter += 1
        
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
            secret: bytes = None
    ):
        """
        :param bam_iter:
        :param bdiff_iter:
        :param out_bam_file:
        :param secret: secret key used for encryptiuon of unmapped alignments
            Must be present when iterating over unmapped alignments.
        :return:
        """
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
                    if alignment is None:
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
                    if alignment is not None:
                        print('last alignment: %s' % self.alignment2str(self._prev_alignment))
                    print("total of %d alignments processed" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, diff.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceding alignments are written
                    self.__write_alignment(out_bam_file, alignment)
                else:
                    # append to queue
                    alignment_queue.append(AlignedVariant(alignment))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, diff.index):
                self.__unmutate_pos(alignment_queue, diff.mut_map)
                # done with this diff, read next
                prev_diff = diff
                diff = next(bdiff_iter)
                if self._verbose and bdiff_iter.counter % 100000 == 0:
                    print('%d DIFF records done' % bdiff_iter.counter)
                if diff is not None:
                    # not the end of DIFF file
                    self.__write_before_index(out_bam_file, alignment_queue, diff.index)
                    # update alignment queue
                    for i in range(len(alignment_queue)):
                        if not alignment_queue[i].alignment.is_unmapped:
                            alignment_variant = cmn.vac_aligned_variant(alignment_queue[i].alignment, diff, False)
                            if alignment_variant is not None:
                                alignment_queue[i] = alignment_variant
                elif self._verbose:
                    print('last diff: %s' % prev_diff)
            
            else:  # alignment is covering diff position
                alignment_variant = cmn.vac_aligned_variant(alignment, diff, False)
                if alignment_variant is not None:
                    alignment_queue.append(alignment_variant)
                alignment = next(bam_iter)
                self.covering_counter += 1
        
        # noinspection PyAttributeOutsideInit
        self.diff_counter = bdiff_iter.counter
    
    def __unmutate_pos(self, variant_queue: list, mut_map: dict):
        for variant in variant_queue:
            if variant.is_present():
                # alignment has vac to mutate
                self._mutate_variant(variant, mut_map)
                # done with current vac
                variant.clear()
    
    def alignment2str(self, alignment) -> str:
        if alignment is None:
            return "None"
        if cmn.is_placed_alignment(alignment):
            ref_name = alignment.reference_name
            ref_start = alignment.reference_start
            # unplaced alignment
            return '#%d %s:%d' % (self._fai.pos2index(ref_name, ref_start), ref_name, ref_start)
        else:
            return alignment.query_sequence
    
    @staticmethod
    def __encrypt_unmapped(alignment: pysam.AlignedSegment, secret: bytes):
        """
        Stream cipher encryption / decryption.
        alignment + secret => encrypted_alignment
        encrypted_alignment + secret => alignment
        :param alignment:
        :param secret:
        :return: encrypter/decrypted alignment
        """
        if alignment.is_unmapped:
            if secret is None:
                raise ValueError('Secret key must be present when unmapped alignments are iterated.')
            
            # use 64B long hash (encrypts 256 bases)
            sha512 = hashlib.sha512()
            sha512.update(secret + alignment.query_name.encode())
            mut_seq = cmn.stream_cipher(alignment.query_sequence, sha512.digest())

            # change and preserve quality TODO: maybe do something else with the quality?
            quality = alignment.query_qualities
            alignment.query_sequence = mut_seq
            alignment.query_qualities = quality
    
    def __mutate_pos(
            self,
            bdiff_io: bdiff.BdiffIO,
            variant_queue: list,
            vac: po.VariantPosition,
            rnd: VeryRandom
    ) -> bool:
        """
        Mutate alignments at the variant allele count position.
        :param variant_queue: list of SnvAlignment
        :param vac: variant allele count list
        :return: True if NS mutation has occured
        """
        if isinstance(vac, po.SnvOccurrence):
            is_mutated = self._mutate_snv_pos(bdiff_io, variant_queue, vac, rnd)
        elif isinstance(vac, po.IndelOccurrence):
            is_mutated = self._mutate_indel_pos(bdiff_io, variant_queue, vac, rnd)
        else:
            raise ValueError("%s is not VAC record instance" % type(vac).__name__)
        
        return is_mutated
    
    def _mutate_snv_pos(
            self,
            bdiff_io: bdiff.BdiffIO,
            variant_queue: list,
            vac: po.SnvOccurrence,
            rnd: VeryRandom
    ) -> bool:
        is_mutated = False
        variant_seqs = cmn.variant_seqs(variant_queue)
        
        alt_freqs = cmn.base_freqs(variant_seqs)
        mut_map = cmn.snv_mut_map(alt_freqs=alt_freqs, ref_freqs=vac.freqs, rnd=rnd)
        
        for variant in variant_queue:  # type: AlignedVariant
            if variant.is_present():
                # alignment has vac to mutate
                is_mutated |= self._mutate_variant(variant, mut_map)
                # done with current vac
                variant.clear()
        
        if is_mutated:
            # at least one alignment has been mutated
            bdiff_io.write_snv(vac.index, vac.ref_id, mut_map)
        
        return is_mutated
    
    def _mutate_indel_pos(
            self,
            bdiff_io: bdiff.BdiffIO,
            variant_queue: list,
            vac: po.IndelOccurrence,
            rnd: VeryRandom
    ) -> bool:
        is_mutated = False
        variant_seqs = cmn.variant_seqs(variant_queue)
        
        alt_freq_map = cmn.freq_map(variant_seqs)
        mut_map = cmn.indel_mut_map(
            alt_freq_map=alt_freq_map,
            ref_freq_map=dict(zip(vac.seqs, vac.freqs)),
            rnd=rnd
        )
        
        for variant in variant_queue:  # type: AlignedVariant
            if variant.is_present():
                # alignment has vac to mutate
                is_mutated |= self._mutate_variant(variant, mut_map)
                # done with current vac
                variant.clear()
        
        if is_mutated:
            bdiff_io.write_indel(vac.index, vac.ref_seq, mut_map)
        
        return is_mutated
    
    def _mutate_variant(self, variant: AlignedVariant, mut_map: dict) -> bool:
        """
        Mutate alignment by mutation map at SNV position.
        :param variant:
        :param mut_map: mutation map
        :return: True if variant mutation is non-synonymous
        """
        is_mutated = False
        
        # alignment is mapped at variant position
        mut_seq = mut_map[variant.seq]
        
        if variant.seq != mut_seq:
            # non-synonymous mutation
            self.mut_counter += 1
            is_mutated = True
            variant.seq = mut_seq
        
        return is_mutated
    
    def __write_before_index(self, out_bam_file, variant_queue, index):
        """
        Write snv alignments while their end reference position is before index.
        :param out_bam_file:
        :param variant_queue:
        :param index:
        :return:
        """
        # TODO optimize, search for index, then cut the array
        tmp_queue = variant_queue[:]
        for variant in tmp_queue:
            is_mapped = not variant.alignment.is_unmapped
            if is_mapped and not self.__is_before_index(variant.alignment, index):
                # break to keep alignments in order
                # otherwise shorter alignment could be written before longer alignment
                break
            
            self.__write_alignment(out_bam_file, variant.alignment)
            # remove written alignment
            variant_queue.remove(variant)
