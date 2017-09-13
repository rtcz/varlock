import pysam


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


class VacSnvRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, freqs: tuple):
        super().__init__(index, ref_name, ref_pos)
        self.freqs = freqs


class VacIndelRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, freqs: list, seqs: list):
        super().__init__(index, ref_name, ref_pos)
        self.freqs = freqs
        self.seqs = seqs


# TODO separate file, tests ?
# noinspection PyUnresolvedReferences
class AlignedVariant:
    def __init__(self, alignment: pysam.AlignedSegment, pos: int = None, end_pos: int = None):
        """
        :param alignment:
        :param pos: position of variation in alignment sequence
        :param end_pos: position after variation in alignment sequence
        """
        self.alignment = alignment
        self._pos = pos
        self._end_pos = end_pos
        
        if self._pos is not None and self._end_pos is None:
            # SNV
            self._end_pos = self._pos + 1
    
    def is_present(self):
        return self._pos is not None and self._end_pos is not None
    
    def clear(self):
        self._pos = None
        self._end_pos = None
    
    @property
    def seq(self):
        """
        Get variant sequence.
        """
        assert self._pos < self._end_pos
        return self.alignment.query_sequence[self._pos:self._end_pos]
    
    @seq.setter
    def seq(self, seq):
        """
        Set variant sequence.
        :param value: new sequence
        """
        assert len(seq) > 0
        mut_seq = self.alignment.query_sequence[:self._pos]
        mut_seq += seq
        mut_seq += self.alignment.query_sequence[self._end_pos:]
        self.alignment.query_sequence = mut_seq


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
