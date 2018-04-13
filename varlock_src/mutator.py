import hashlib
import typing

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
    
    def __write_alignment(self, bam_file, alignment: pysam.AlignedSegment, is_mutated: bool):
        """
        :param alignment: pysam.AlignedSegment
        :return:
        """
        bam_file.write(alignment)
        self.alignment_counter += 1
        self.alignment_mut_counter += is_mutated
        
        if self._verbose and self.alignment_counter % 10000 == 0:
            print("%d alignments done" % self.alignment_counter)
    
    def __init_counters(self):
        self.alignment_counter = 0  # written alignments
        self.diff_counter = 0  # processed diff records
        self.mut_counter = 0  # all mutations
        self.alignment_mut_counter = 0  # mutated alignments
        self.vac_counter = 0  # read vac records (variants)
        
        self.covering_counter = 0  # alignments with variant positions
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
                
                for variant_alignment in alignment_queue:  # type: AlignedVariant
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(mut_bam_file, variant_alignment.alignment, variant_alignment.is_mutated)
                
                # write remaining alignments (in EOF VAC case only)
                while alignment is not None:
                    self.__encrypt_unmapped(alignment, secret)
                    self.__write_alignment(mut_bam_file, alignment, False)
                    alignment = next(bam_iter)
                
                if self._verbose:
                    print("total of %d alignments done" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, vac.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceding alignments are written
                    self.__write_alignment(mut_bam_file, alignment, False)
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
                if self._verbose and vac_iter.counter % 10000 == 0:
                    print('%d VAC records done' % vac_iter.counter)
                # not the end of VAC file
                if vac is not None:
                    # write alignments to mutated BAM
                    self.__write_before_index(mut_bam_file, alignment_queue, vac.index)
                    # update alignment queue
                    for i in range(len(alignment_queue)):
                        alignment_queue[i] = self.create_aligned_variant(
                            alignment_queue[i].alignment,
                            vac,
                            True,
                            alignment_queue[i].is_mutated
                        )
                
                elif self._verbose:
                    print('last vac: %s' % prev_vac)
            
            else:  # alignment is covering vac position
                alignment_queue.append(self.create_aligned_variant(alignment, vac, True))
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
                
                for variant_alignment in alignment_queue:
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(out_bam_file, variant_alignment.alignment, variant_alignment.is_mutated)
                
                # write remaining alignments (in EOF DIFF case only)
                while alignment is not None:
                    self.__encrypt_unmapped(alignment, secret)
                    self.__write_alignment(out_bam_file, alignment, False)
                    alignment = next(bam_iter)
                
                if self._verbose:
                    print("total of %d alignments done" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, diff.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceding alignments are written
                    self.__write_alignment(out_bam_file, alignment, False)
                else:
                    # append to queue
                    alignment_queue.append(AlignedVariant(alignment))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, diff.index):
                self.__unmutate_pos(alignment_queue, diff.mut_map)
                # done with this diff, read next
                prev_diff = diff
                diff = next(bdiff_iter)
                if self._verbose and bdiff_iter.counter % 10000 == 0:
                    print('%d DIFF records done' % bdiff_iter.counter)
                if diff is not None:
                    # not the end of DIFF file
                    self.__write_before_index(out_bam_file, alignment_queue, diff.index)
                    # update alignment queue
                    for i in range(len(alignment_queue)):
                        alignment_queue[i] = self.create_aligned_variant(
                            alignment_queue[i].alignment,
                            diff,
                            False,
                            alignment_queue[i].is_mutated
                        )
                
                elif self._verbose:
                    print('last diff: %s' % prev_diff)
            
            else:  # alignment is covering diff position
                alignment_queue.append(self.create_aligned_variant(alignment, diff, False))
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
    
    # def alignment2str(self, alignment) -> str:
    #     if alignment is None:
    #         return "None"
    #     if cmn.is_placed_alignment(alignment):
    #         ref_name = alignment.reference_name
    #         ref_start = alignment.reference_start
    #         return '#%d %s:%d' % (self._fai.pos2index(ref_name, ref_start), ref_name, ref_start)
    #     else:
    #         # unplaced alignment
    #         return alignment.query_sequence
    
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
        for variant in tmp_queue:  # type: AlignedVariant
            is_mapped = not variant.alignment.is_unmapped
            if is_mapped and not self.__is_before_index(variant.alignment, index):
                # break to keep alignments in order
                # otherwise shorter alignment could be written before longer alignment
                break
            
            self.__write_alignment(out_bam_file, variant.alignment, variant.is_mutated)
            # remove written alignment
            variant_queue.remove(variant)
    
    # TODO move back to commons?
    def create_aligned_variant(
            self,
            alignment: pysam.AlignedSegment,
            vac: typing.Union[po.VariantPosition, po.VariantDiff],
            vac_occurrence: bool,
            is_mutated: bool = False
    ):
        """
        Factory that creates AlignedVariant from pysam alignment and VAC record.
        Alignment and VAC are assumed to be of the same reference.
        :param alignment: aligned read
        :param vac: variant position of diff record if vac_occurrence is False
        :param vac_occurrence: if we compare to Snv/IndelOccurrence (encrypt) or to Snv/IndelDiff (decrpyt)
        :param is_mutated:
        """
        assert alignment.reference_name == vac.ref_name
        
        SnpCompare = po.SnvOccurrence
        IndelCompare = po.IndelOccurrence
        if not vac_occurrence:
            SnpCompare = po.SnvDiff
            IndelCompare = po.IndelDiff
        
        if alignment.is_unmapped or (vac.ref_pos >= alignment.reference_end or vac.ref_pos < alignment.reference_start):
            # variant is unmapped or either after or before the alignment
            variant = AlignedVariant(alignment, is_mutated=is_mutated)
        else:
            pos = cmn.ref_pos2seq_pos(alignment, vac.ref_pos)
            if pos is None:
                # message = "WARNING: reference position %d on alignment with range <%d,%d> not found (possible deletion in bam read)" \
                #           % (vac.ref_pos, alignment.reference_start, alignment.reference_end - 1)
                # print(message)
                variant = AlignedVariant(alignment, is_mutated=is_mutated)
            elif isinstance(vac, SnpCompare):
                variant = AlignedVariant(alignment, pos, is_mutated=is_mutated)
            elif isinstance(vac, IndelCompare):
                words = vac.seqs if vac_occurrence else list(vac.mut_map.values())
                end_pos = pos + cmn.max_match_cigar(alignment, pos, vac.ref_pos, words, vac.ref_seq)
                if end_pos > pos:
                    # indel was found
                    if vac_occurrence:
                        assert len(vac.seqs)
                    
                    if alignment.query_name == 'ERR015528.28046840' and alignment.reference_start == 356649:
                        print('YYY')
                        print('vac_pos: %s' % vac.ref_pos)
                        print('vac_seqs: %s' % vac.seqs)
                        print('ref_seq: %s' % vac.ref_seq)
                        # print('infer_query_length %d' % alignment.infer_query_length())
                        # print('query length %d' % len(alignment.query_sequence))
                        # print(alignment.cigarstring)
                        # print(alignment.query_sequence)
                        print('YYY')
                        # exit(0)
                    
                    variant = AlignedVariant(alignment, pos, end_pos, vac.ref_seq, is_mutated)
                else:
                    # match not found or at least one variant exceeds alignment end
                    variant = AlignedVariant(alignment, is_mutated=is_mutated)
            else:
                raise ValueError("%s is not %s/%s record instance" % (
                    type(vac).__name__, type(SnpCompare).__name__, type(IndelCompare).__name__,))
        
        return variant
