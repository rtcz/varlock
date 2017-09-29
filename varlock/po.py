import pysam
import numpy as np

class GenomicPosition(object):
    def __init__(self, index, ref_name, ref_pos):
        self.index = index
        self.ref_name = ref_name
        self.ref_pos = ref_pos
    
    def __str__(self):
        return '#%d %s:%d' % (self.index, self.ref_name, self.ref_pos)


class DiffSnvRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, mut_map: dict):
        super().__init__(index, ref_name, ref_pos)
        self.mut_map = mut_map


class DiffIndelRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, mut_map: dict):
        super().__init__(index, ref_name, ref_pos)
        self.mut_map = mut_map
        
    def ref_seq(self):
        return


class VacSnvRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, freqs: tuple):
        super().__init__(index, ref_name, ref_pos)
        self.freqs = freqs


class VacIndelRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, freqs: list, seqs: list):
        super().__init__(index, ref_name, ref_pos)
        self.freqs = freqs
        self.seqs = seqs
        
    def ref_seq(self):
        # TODO embed reference id into VAC file
        # noinspection PyTypeChecker
        return self.seqs[np.argmax(self.freqs)]


# TODO separate file, tests ?
# noinspection PyUnresolvedReferences
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
        :param pos: position of variation in AlignedSegment.query_sequence
        :param end_pos: position one base after INDEL variation in AlignedSegment.query_sequence
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
        :param value: new sequence
        """
        # TODO do something about MD string if present
        # TODO update quality string
        assert len(seq) > 0
        if self._is_snv:
            # TODO treat [X, =] OP cases
            # expecting only M OP now
            mut_seq = self.alignment.query_sequence[:self._pos]
            mut_seq += seq
            mut_seq += self.alignment.query_sequence[self._end_pos:]
        elif self._is_indel:
            if self._first_bases_match(self._seqs):
                # VCF indels should always have single matched base preffix
                tmp_end_pos = self._end_pos + 1
            else:
                tmp_end_pos = self._end_pos
            
            tmp_cigar = Cigar.del_subrange(
                self.alignment.cigartuples,
                self._pos,
                tmp_end_pos
            )
            # assuming that first seq is the reference seq
            variant_cigar = Cigar.variant(self._ref_seq, seq)
            for tpl in variant_cigar:
                tmp_cigar = Cigar.place_op(
                    tmp_cigar,
                    self._pos,
                    tpl[0],
                    tpl[1]
                )
            
            self.alignment.cigartuples = tmp_cigar
    
    @staticmethod
    def _first_bases_match(seqs: list):
        assert len(seqs) > 0
        for i in range(len(seqs) - 1):
            if seq[i] != seq[i + 1]:
                return False
        
        return True


class FaiRecord:
    def __init__(self, index: int, name: str, start: int, length: int):
        """
        :param index: reference id
        :param name: reference name
        :param start: 0-based first position
        :param length: reference length - number of bases
        """
        self.id = index
        self.name = name
        self.start = start
        self.length = length
    
    def __str__(self):
        return '#%d %s %d %d' % (self.id, self.name, self.start, self.length)
