import hashlib

import pysam

import varlock_src.bdiff as bdiff
import varlock_src.common as cmn
import varlock_src.iters as iters
import varlock_src.po as po
from varlock_src.fasta_index import FastaIndex
from varlock_src.random import VeryRandom
from varlock_src.variant import AlignmentAllele


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
        
        # TODO optimization: seek to first alignment covering variant to skip all preceding records
        variant = next(vac_iter)  # type: po.VariantOccurrence
        
        if self._verbose:
            print('first variant: %s' % variant)
        
        bdiff_io = bdiff.BdiffIO()
        while True:
            if variant is None or alignment is None:
                # finish
                if self._verbose:
                    if variant is None:
                        print("EOF VAC")
                    else:
                        print("EOF BAM")
                
                # last mutation
                if variant is not None:
                    # noinspection PyTypeChecker
                    self.__mutate_pos(bdiff_io, alignment_queue, variant, rnd)
                
                for variant_alignment in alignment_queue:  # type: AlignmentAllele
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
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, variant.pos.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceding alignments are written
                    self.__write_alignment(mut_bam_file, alignment, False)
                else:
                    # append to queue
                    alignment_queue.append(AlignmentAllele(alignment))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, variant.pos.index):
                # apply the variant to alignment
                self.__mutate_pos(bdiff_io, alignment_queue, variant, rnd)
                # done with this variant, read next
                prev_vac = variant
                variant = next(vac_iter)
                if self._verbose and vac_iter.counter % 10000 == 0:
                    print('%d VAC records done' % vac_iter.counter)
                # not the end of VAC file
                if variant is not None:
                    # write alignments to mutated BAM
                    self.__write_before_index(mut_bam_file, alignment_queue, variant.pos.index)
                    # update alignment queue
                    for i in range(len(alignment_queue)):
                        alignment_queue[i] = AlignmentAllele(
                            alignment_queue[i].alignment,
                            variant,
                            alignment_queue[i].is_mutated
                        )
                
                elif self._verbose:
                    print('last variant: %s' % prev_vac)
            
            else:  # alignment is covering variant position
                alignment_queue.append(AlignmentAllele(alignment, variant))
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
        variant = next(bdiff_iter)  # type: po.VariantDiff
        
        if self._verbose:
            print('first variant: %s' % variant)
        
        while True:
            if variant is None or alignment is None:
                # finish
                if self._verbose:
                    if variant is None:
                        print("EOF DIFF")
                    if alignment is None:
                        print("EOF BAM")
                
                # last mutation
                if variant is not None:
                    # noinspection PyTypeChecker
                    self.__unmutate_pos(alignment_queue, variant.mut_map)
                
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
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, variant.pos.index):
                self.__encrypt_unmapped(alignment, secret)
                if len(alignment_queue) == 0:
                    # good to go, all preceding alignments are written
                    self.__write_alignment(out_bam_file, alignment, False)
                else:
                    # append to queue
                    alignment_queue.append(AlignmentAllele(alignment))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, variant.pos.index):
                self.__unmutate_pos(alignment_queue, variant.mut_map)
                # done with this variant, read next
                prev_diff = variant
                variant = next(bdiff_iter)
                if self._verbose and bdiff_iter.counter % 10000 == 0:
                    print('%d DIFF records done' % bdiff_iter.counter)
                if variant is not None:
                    # not the end of DIFF file
                    self.__write_before_index(out_bam_file, alignment_queue, variant.pos.index)
                    # update alignment queue
                    for i in range(len(alignment_queue)):
                        alignment_queue[i] = AlignmentAllele(
                            alignment_queue[i].alignment,
                            variant,
                            alignment_queue[i].is_mutated
                        )
                
                elif self._verbose:
                    print('last variant: %s' % prev_diff)
            
            else:  # alignment is covering variant position
                alignment_queue.append(AlignmentAllele(alignment, variant))
                alignment = next(bam_iter)
                self.covering_counter += 1
        
        # noinspection PyAttributeOutsideInit
        self.diff_counter = bdiff_iter.counter
    
    def __unmutate_pos(self, allele_queue: list, mut_map: dict):
        for allele in allele_queue:  # type: AlignmentAllele
            if allele.is_known:
                # alignment has vac to mutate
                self._mutate_allele(allele, mut_map)
                # done with current vac
                # TODO
                # variant.clear()
    
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
            allele_queue: list,
            variant: po.VariantOccurrence,
            rnd: VeryRandom
    ) -> bool:
        """
        Mutate alignments at the variant allele count position.
        :param allele_queue:
        :param variant: variant allele count list
        :return: True if NS mutation has occured
        """
        if variant.is_type(po.VariantType.SNV):
            is_mutated = self._mutate_snv_pos(bdiff_io, allele_queue, variant, rnd)
        elif variant.is_type(po.VariantType.INDEL):
            is_mutated = self._mutate_indel_pos(bdiff_io, allele_queue, variant, rnd)
        else:
            raise ValueError("unnknown variant type")
        
        return is_mutated
    
    def _mutate_snv_pos(
            self,
            bdiff_io: bdiff.BdiffIO,
            allele_queue: list,
            variant: po.VariantOccurrence,
            rnd: VeryRandom
    ) -> bool:
        is_mutated = False
        variant_seqs = cmn.variant_seqs(allele_queue)
        
        alt_freqs = cmn.base_freqs(variant_seqs)
        mut_map = cmn.snv_mut_map(alt_freqs=alt_freqs, ref_freqs=variant.freqs, rnd=rnd)
        
        for allele in allele_queue:  # type: AlignmentAllele
            if allele.is_known:
                # alignment has vac to mutate
                is_mutated |= self._mutate_allele(allele, mut_map)
                # done with current vac
                # TODO
                # variant.clear()
        
        if is_mutated:
            # at least one alignment has been mutated
            # TODO ref_allele_id
            bdiff_io.write_snv(variant.pos.index, cmn.BASES.index(variant.ref_allele), mut_map)
        
        return is_mutated
    
    def _mutate_indel_pos(
            self,
            bdiff_io: bdiff.BdiffIO,
            allele_queue: list,
            variant: po.VariantOccurrence,
            rnd: VeryRandom
    ) -> bool:
        is_mutated = False
        variant_seqs = cmn.variant_seqs(allele_queue)
        
        alt_freq_map = cmn.freq_map(variant_seqs)
        
        # if vac.ref_pos == 106700:
        #     print(alt_freq_map)
        #     print(dict(zip(vac.seqs, vac.freqs)))
        #     exit(0)
        
        mut_map = cmn.indel_mut_map(
            alt_freq_map=alt_freq_map,
            ref_freq_map=dict(zip(variant.alleles, variant.freqs)),
            rnd=rnd
        )
        
        for allele in allele_queue:  # type: AlignmentAllele
            if allele.is_known:
                # alignment has vac to mutate
                is_mutated |= self._mutate_allele(allele, mut_map)
                # done with current vac
                # TODO
                # variant.clear()
        
        if is_mutated:
            bdiff_io.write_indel(variant.pos.index, variant.ref_allele, mut_map)
        
        return is_mutated
    
    def _mutate_allele(self, variant: AlignmentAllele, mut_map: dict) -> bool:
        """
        Mutate alignment by mutation map at SNV position.
        :param variant:
        :param mut_map: mutation map
        :return: True if variant mutation is non-synonymous
        """
        is_mutated = False
        
        mut_seq = mut_map[variant.allele]
        if variant.allele != mut_seq:
            # non-synonymous mutation
            self.mut_counter += 1
            is_mutated = True
            variant.allele = mut_seq
        
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
        for variant in tmp_queue:  # type: AlignmentAllele
            is_mapped = not variant.alignment.is_unmapped
            if is_mapped and not self.__is_before_index(variant.alignment, index):
                # break to keep alignments in order
                # otherwise shorter alignment could be written before longer alignment
                break
            
            self.__write_alignment(out_bam_file, variant.alignment, variant.is_mutated)
            # remove written alignment
            variant_queue.remove(variant)
