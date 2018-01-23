import pysam
import numpy as np
from array import array

# TODO ? tests
from varlock_src.cigar import Cigar


def qual_str2array(qual_str: str):
    return array('B', map(lambda x: ord(x) - 33, qual_str))


def qual_array2str(qual_array: array):
    if qual_array is None:
        return ''
    return ''.join(map(lambda x: chr(x + 33), qual_array))


class AlignedVariant:
    def __init__(
            self,
            alignment: pysam.AlignedSegment,
            pos: int = None,
            end_pos: int = None,
            ref_seq: str = None
    ):
        """
        :param alignment:
        :param pos: position of variation in AlignedSegment.query_sequence in alternative sequence
        :param end_pos: position one base after INDEL variation in AlignedSegment.query_sequence in alternative sequence
        :param ref_seq: reference sequence of the variant
        """
        self.alignment = alignment
        self._pos = pos
        self._end_pos = end_pos
        self._is_snv = False
        self._is_indel = False
        self._ref_seq = ref_seq

        if self._pos is not None:
            if self._end_pos is not None:
                assert self._pos < self._end_pos
                assert ref_seq is not None
                self._is_indel = True
            else:
                self._is_snv = True
                self._end_pos = self._pos + 1

    def is_present(self):
        return self._is_snv or self._is_indel

    # TODO try to remove, pos and end_pos should be immutable (so as variant presence)
    def clear(self):
        self._is_snv = False
        self._is_indel = False
        self._pos = None
        self._end_pos = None

    @property
    def seq(self):
        """
        Get variant sequence.
        """
        if self.is_present():
            return self.alignment.query_sequence[self._pos:self._end_pos]
        else:
            return None

    @seq.setter
    def seq(self, seq):
        """
        Set variant sequence.
        :param seq: new sequence
        """
        # save quality

        from copy import copy
        quality_seq = copy(self.alignment.query_qualities)
        if self._is_indel:
            mean_qual = 34
            if self.alignment.query_qualities is not None:
                if len(self.alignment.query_qualities) > 0:
                    mean_qual = int(np.round(np.mean(self.alignment.query_qualities) + 0.5))
                len_seq = len(seq)
                len_var = self._end_pos - self._pos
                len_overlap = min(len_var, len_seq)
                quality_seq = self.alignment.query_qualities[:self._pos + len_overlap] + array('B', [mean_qual] * (len_seq - len_var)) + self.alignment.query_qualities[self._end_pos:]
                '''print(array('B', [mean_qual] * (len(seq) - self._end_pos + self._pos)), mean_qual, seq, len(seq), self.alignment.query_sequence, self._pos, self._end_pos,
                      self.alignment.query_sequence[:self._pos],
                      self.alignment.query_sequence[self._end_pos:], qual_array2str(self.alignment.query_qualities), qual_array2str(self.alignment.query_qualities[:self._pos]),
                      qual_array2str(self.alignment.query_qualities[self._end_pos:]), qual_array2str(quality_seq))'''
                # print(len(quality_seq), len(self.alignment.query_qualities[:self._pos + len_overlap]), len(self.alignment.query_qualities[:self._end_pos]))

        # fix sequence
        mut_seq = self.alignment.query_sequence[:self._pos]
        mut_seq += seq
        mut_seq += self.alignment.query_sequence[self._end_pos:]
        self.alignment.query_sequence = mut_seq

        # print('\n' + str(self.alignment.query_qualities), qual_array2str(quality_seq), self._is_snv)

        # fix vanishing quality:
        self.alignment.query_qualities = quality_seq

        # print(self.alignment.query_alignment_sequence)
        # print(self._pos)
        # print(self._ref_seq)
        # print()

        # TODO do something about MD string if present
        # TODO update quality string
        assert len(seq) > 0
        if self._is_snv:
            # TODO treat [X, =] OP cases
            # TODO at least assert that corresponding cigar letter is M
            # expecting only M OP now - it does not change with different SNV
            pass
        elif self._is_indel:
            # if self._first_bases_match(self._seqs):
            #     # VCF indels should always have single matched base prefix
            #     tmp_end_pos = self._end_pos + 1
            # else:
            #     tmp_end_pos = self._end_pos

            # remove point of mutation
            tmp_cigar = Cigar.del_subrange(
                self.alignment.cigartuples,
                self._pos,
                self._end_pos
            )

            # generate cigar for change
            variant_cigar = Cigar.variant(self._ref_seq, seq)
            pos = self._pos
            for tpl in variant_cigar:
                tmp_cigar = Cigar.place_op(tmp_cigar, pos, tpl[0], tpl[1])
                # consumes query?
                if tpl[0] not in [Cigar.OP_DEL, Cigar.OP_REF_SKIP]:
                    pos += tpl[1]

            self.alignment.cigartuples = tmp_cigar

    @staticmethod
    def _first_bases_match(seqs: list):
        assert len(seqs) > 0
        for i in range(len(seqs) - 1):
            if seqs[i] != seqs[i + 1]:
                return False

        return True
