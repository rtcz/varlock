import hashlib
import typing

import numpy as np
import pysam

import src.bdiff as bdiff
import src.common as cmn
import src.iters as iters
import src.po as po
from src.alignment import AlleleAlignment, pileup_alleles
from src.fasta_index import FastaIndex
from src.random import VeryRandom


class Mutator:
    MIN_ALLELE_RATIO = 0.2

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

        if self._verbose and self.alignment_counter % 100000 == 0:
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
            variant_iter: iters.VariantIterator,
            bam_iter: iters.FullBamIterator,
            secret: bytes,
            rnd: VeryRandom
    ) -> bdiff.BdiffIO:
        """
        :param mut_bam_file:
        :param variant_iter:
        :param bam_iter:
        :param secret:
        :param rnd: (secure) random generator
        :return: BdiffIO object
        BDIFF records are written only if bases at variant position differ after mutation.
        """
        self.__init_counters()

        alignment_queue = []
        alignment = next(bam_iter)
        # TODO optimization: seek to first alignment covering variant to skip all preceding records
        variant = variant_iter.next_after()

        if self._verbose:
            print('first variant: %s' % variant)

        # # TODO temporary code - remove
        # self._test_counter = 0
        # self._pos_file = open('/home/hekel/projects/python/varlock/analysis/exome/analysis/vac_cov.pos', 'w')
        # self._freqs_file = open('/home/hekel/projects/python/varlock/analysis/exome/analysis/freqs_v5.tsv', 'w')
        # self._bed_file_lines = open('in/chr20.bed', 'rt').readlines()
        # self._bed_id = 0

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
                    self._mutate_pos(bdiff_io, alignment_queue, variant, rnd)

                for allele_alignment in alignment_queue:  # type: AlleleAlignment
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(mut_bam_file, allele_alignment.alignment, allele_alignment.is_mutated)

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
                    alignment_queue.append(AlleleAlignment(alignment))

                alignment = next(bam_iter)

            elif self.__is_after_index(alignment, variant.pos.index):
                # apply the variant to alignment
                # if len(alignment_queue) == 0:
                #     self.test_counter += 1

                self._mutate_pos(bdiff_io, alignment_queue, variant, rnd)
                # done with this variant, read next
                prev_variant = variant
                variant = variant_iter.next_after()
                if self._verbose and variant_iter.counter % 100000 == 0:
                    print('%d VAC records done' % variant_iter.counter)
                # not the end of VAC file
                if variant is not None:
                    # write alignments to mutated BAM
                    self.__write_before_index(mut_bam_file, alignment_queue, variant.pos.index)
                    # update alignment queue with new variant
                    for i in range(len(alignment_queue)):
                        # if alignment.reference_name != variant.pos.ref_name:
                        #     print(f'alignment {alignment.reference_name}:{alignment.reference_start} '
                        #           f'has foreign variant {variant}')
                        #     print(f'prev variant {prev_variant}')

                        alignment_queue[i] = AlleleAlignment(
                            alignment_queue[i].alignment,
                            variant,
                            alignment_queue[i].is_mutated
                        )

                elif self._verbose:
                    print('last variant: %s' % prev_variant)

            else:  # alignment is covering variant position
                alignment_queue.append(AlleleAlignment(alignment, variant))
                alignment = next(bam_iter)
                self.covering_counter += 1

        # noinspection PyAttributeOutsideInit
        self.vac_counter = variant_iter.counter
        # noinspection PyAttributeOutsideInit
        self.diff_counter = bdiff_io.snv_count + bdiff_io.indel_count

        # print('TEST: %d' % self._test_counter)

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

                for allele_alignment in alignment_queue:  # type: AlleleAlignment
                    # unmapped alignments in queue are already encrypted
                    self.__write_alignment(out_bam_file, allele_alignment.alignment, allele_alignment.is_mutated)

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
                    alignment_queue.append(AlleleAlignment(alignment))

                alignment = next(bam_iter)

            elif self.__is_after_index(alignment, variant.pos.index):
                self.__unmutate_pos(alignment_queue, variant.mut_map)
                # done with this variant, read next
                prev_variant = variant
                variant = next(bdiff_iter)
                if self._verbose and bdiff_iter.counter % 10000 == 0:
                    print('%d DIFF records done' % bdiff_iter.counter)
                if variant is not None:
                    # not the end of DIFF file
                    self.__write_before_index(out_bam_file, alignment_queue, variant.pos.index)
                    # update alignment queue
                    for i in range(len(alignment_queue)):
                        alignment_queue[i] = AlleleAlignment(
                            alignment_queue[i].alignment,
                            variant,
                            alignment_queue[i].is_mutated
                        )

                elif self._verbose:
                    print('last variant: %s' % prev_variant)

            else:  # alignment is covering variant position
                alignment_queue.append(AlleleAlignment(alignment, variant))

                alignment = next(bam_iter)

                self.covering_counter += 1

        # noinspection PyAttributeOutsideInit
        self.diff_counter = bdiff_iter.counter

    def __unmutate_pos(self, allele_queue: list, mut_map: dict):
        for allele in allele_queue:  # type: AlleleAlignment
            self._mutate_allele(allele, mut_map)

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

            # change and preserve quality
            # TODO: maybe something else with the quality?
            quality = alignment.query_qualities
            alignment.query_sequence = mut_seq
            alignment.query_qualities = quality

    def _mutate_pos(
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
        is_mutated = False
        if any(allele.is_known for allele in allele_queue):

            # TODO temp code
            # BEGIN temp code
            # self._pos_file.write('%d\n' % (variant.pos.ref_pos + 1))
            # END temp code

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
        pileup = pileup_alleles(allele_queue)
        private_counts = np.array(cmn.base_counts(pileup))
        public_counts = np.array(variant.freqs)

        private_freqs = private_counts / private_counts.sum()
        public_freqs = public_counts / public_counts.sum()

        # private_ids = (-private_freqs).argsort()

        personal_alleles = []
        for i in range(len(cmn.BASES)):
            if private_freqs[i] >= self.MIN_ALLELE_RATIO:
                personal_alleles.append(cmn.BASES[i])

        if len(personal_alleles) == 1:
            # homozygote
            from_allele_a = from_allele_b = personal_alleles[0]
        elif len(personal_alleles) == 2:
            # heterozygote
            from_allele_a = personal_alleles[0]
            from_allele_b = personal_alleles[1]
        else:
            # unable to asses
            # self._test_counter += 1
            return False

        prob_matrix = np.outer(public_freqs, public_freqs)
        draw_id = rnd.multirand_index(prob_matrix.reshape(prob_matrix.size))
        to_allele_a = cmn.BASES[int(draw_id / 4)]  # row_id
        to_allele_b = cmn.BASES[draw_id % 4]  # col_id

        homo2homo = False
        homo2hetero = False
        hetero2homo = False
        hetero2hetero = False

        is_masked = False
        mut_map = None
        rev_mut_map = None
        if to_allele_a == to_allele_b:
            # to homo
            if from_allele_a == from_allele_b:
                # from homo
                homo2homo = True
                if from_allele_a != to_allele_a:
                    # homo to other homo
                    mut_map = dict(zip(cmn.BASES, cmn.BASES))
                    mut_map[from_allele_a] = to_allele_a
                    mut_map[to_allele_a] = from_allele_a
                else:
                    # homo to same homo
                    pass
            else:
                # TODO not reversible yet
                hetero2homo = True
                mut_map = dict(zip(cmn.BASES, cmn.BASES))
                mut_map[from_allele_a] = to_allele_a
                mut_map[from_allele_b] = to_allele_a

        else:
            # to hetero
            if from_allele_a == from_allele_b:
                # from homo
                homo2hetero = True

                mut_map_a = {from_allele_a: to_allele_a}
                mut_map_b = {from_allele_a: to_allele_b}

                # find remaining bases that must be mapped
                from_bases = []
                for i in range(len(cmn.BASES)):
                    if cmn.BASES[i] != from_allele_a and private_freqs[i] != 0:
                        from_bases.append(cmn.BASES[i])

                # find remaining mappings
                to_bases = list(cmn.BASES)
                to_bases.remove(to_allele_a)
                to_bases.remove(to_allele_b)

                if len(from_bases) <= len(to_bases):
                    # all bases can be mapped
                    for i in range(len(from_bases)):
                        mut_map_a[from_bases[i]] = to_bases[i]
                        mut_map_b[from_bases[i]] = to_bases[i]

                    for aligned_allele in allele_queue:  # type: AlleleAlignment
                        # select random allele from the hetero pair
                        if rnd.random() > 0.5:
                            is_masked |= self._mutate_allele(aligned_allele, mut_map_a)
                        else:
                            is_masked |= self._mutate_allele(aligned_allele, mut_map_b)

                    rev_mut_map = dict(zip(cmn.BASES, cmn.BASES))
                    rev_mut_map[to_allele_a] = from_allele_a
                    rev_mut_map[to_allele_b] = from_allele_a

                    for base in cmn.BASES:
                        if base in mut_map_a:
                            rev_mut_map[mut_map_a[base]] = base
                        elif base in mut_map_b:
                            rev_mut_map[mut_map_b[base]] = base

                else:
                    # TODO can not be mapped by the current method
                    # self._test_counter += 1
                    return False

            else:
                hetero2hetero = True
                is_synonymous = from_allele_a == to_allele_a and from_allele_b == to_allele_b or \
                                from_allele_a == to_allele_b and from_allele_b == to_allele_a

                if not is_synonymous:
                    mut_map = dict(zip(cmn.BASES, cmn.BASES))
                    # select the least changes alternative
                    if from_allele_a == to_allele_b or from_allele_b == to_allele_a:
                        mut_map[from_allele_a] = to_allele_b
                        mut_map[from_allele_b] = to_allele_a

                        mut_map[to_allele_a] = from_allele_b
                        mut_map[to_allele_b] = from_allele_a
                    else:
                        mut_map[from_allele_a] = to_allele_a
                        mut_map[from_allele_b] = to_allele_b

                        mut_map[to_allele_a] = from_allele_a
                        mut_map[to_allele_b] = from_allele_b

        if mut_map is not None and not homo2hetero:
            for aligned_allele in allele_queue:  # type: AlleleAlignment
                is_masked |= self._mutate_allele(aligned_allele, mut_map)

            if hetero2homo:
                # temp solution, use dummy dict
                rev_mut_map = dict(zip(cmn.BASES, cmn.BASES))
            else:
                rev_mut_map = {value: key for key, value in mut_map.items()}

        if is_masked:
            # at least one alignment has been mutated
            bdiff_io.write_snv(variant.pos.index, cmn.BASES.index(variant.ref_allele), rev_mut_map)

        return is_masked

    def _mutate_indel_pos(
            self,
            bdiff_io: bdiff.BdiffIO,
            allele_queue: list,
            variant: po.VariantOccurrence,
            rnd: VeryRandom
    ) -> bool:
        # # TODO
        # return False

        pileup = pileup_alleles(allele_queue)
        alt_freq_map = cmn.freq_map(pileup)

        mut_map = cmn.indel_mut_map(
            private_freq_map=alt_freq_map,
            public_freq_map=dict(zip(variant.alleles, variant.freqs)),
            rnd=rnd
        )

        # # BEGIN temp code
        # personal_allele = sorted(alt_freq_map.keys(), key=lambda key: alt_freq_map[key], reverse=True)[0]
        # # personal allele is mutated
        # ref_id = variant.alleles.index(variant.ref_allele)
        # max_id = variant.alleles.index(personal_allele)  # type: int
        # personal_freq = variant.freqs[max_id] / sum(variant.freqs)
        #
        # is_mutated = personal_allele != mut_map[personal_allele]
        # is_reference = ref_id == max_id
        #
        # record = '%d\t%f\t%d\t%d\n' % (
        #     variant.pos.index, personal_freq, is_reference, is_mutated
        # )
        # self._freqs_file.write(record)
        # # END temp code

        is_masked = False
        for allele in allele_queue:  # type: AlleleAlignment
            is_masked |= self._mutate_allele(allele, mut_map)

        if is_masked:
            bdiff_io.write_indel(variant.pos.index, variant.ref_allele, mut_map)

        return is_masked

    def _mutate_allele(self, alignment: AlleleAlignment, mut_map: dict) -> bool:
        """
        Mutate alignment by mutation map at SNV position.
        :param alignment:
        :param mut_map: mutation map
        :return: True if variant mutation is non-synonymous
        """
        is_mutated = False

        if alignment.allele is not None:
            # alignment with allele (defined by VAC)
            # mapping must exist
            assert alignment.allele in mut_map
            mut_allele = mut_map[alignment.allele]

            if alignment.allele != mut_allele:
                # non-synonymous mutation
                # print(mut_seq)

                self.mut_counter += 1
                #
                # print(self.mut_counter)
                # print(alignment.variant.pos.index)

                is_mutated = True
                alignment.allele = mut_allele

        return is_mutated

    def __write_before_index(self, out_bam_file, alignment_queue: typing.List[AlleleAlignment], index):
        """
        Write snv alignments while their end reference position is before index.
        :param out_bam_file:
        :param variant_queue:
        :param index:
        :return:
        """
        # TODO optimize, search for index, then cut the array
        tmp_queue = alignment_queue[:]
        for variant in tmp_queue:  # type: AlleleAlignment
            is_mapped = not variant.alignment.is_unmapped
            if is_mapped and not self.__is_before_index(variant.alignment, index):
                # break to keep alignments in order
                # otherwise shorter alignment could be written before longer alignment
                break

            self.__write_alignment(out_bam_file, variant.alignment, variant.is_mutated)
            # remove written alignment
            alignment_queue.remove(variant)
