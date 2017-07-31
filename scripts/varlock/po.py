class GenomicPosition(object):
    def __init__(self, index, ref_name, ref_pos):
        self.index = index
        self.ref_name = ref_name
        self.ref_pos = ref_pos
    
    def __str__(self):
        return '#%d %s:%d' % (self.index, self.ref_name, self.ref_pos)


class DiffRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, mut_map: dict):
        super().__init__(index, ref_name, ref_pos)
        self.mut_map = mut_map


class VacSnvRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, ac: list):
        super().__init__(index, ref_name, ref_pos)
        self.ac = ac


class VacIndelRecord(GenomicPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, indel_map: list):
        super().__init__(index, ref_name, ref_pos)
        self.indel_map = indel_map


class SnvAlignment:
    def __init__(self, alignment, pos):
        """
        :param alignment: pysam.AlignedSegment
        :param pos: reference position of SNV
        """
        self.alignment = alignment
        self.pos = pos


class FaiRecord:
    def __init__(self, index, name, start, length):
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
