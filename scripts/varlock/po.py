class SnvPosition(object):
    def __init__(self, index, ref_name, ref_pos):
        self.index = index
        self.ref_name = ref_name
        self.ref_pos = ref_pos


class DiffRecord(SnvPosition):
    def __init__(self, index, ref_name, ref_pos, mut_map):
        super().__init__(index, ref_name, ref_pos)
        self.mut_map = mut_map


class VacRecord(SnvPosition):
    def __init__(self, index, ref_name, ref_pos, ac):
        super().__init__(index, ref_name, ref_pos)
        self.ac = ac


class SnvAlignment:
    def __init__(self, alignment, pos):
        """
        :param alignment: pysam.AlignedSegment
        :param pos: reference position of SNV
        """
        self.alignment = alignment
        self.pos = pos


class FaiRecord:
    def __init__(self, name, start, length):
        """
        :param name: reference name
        :param start: 0-based first position
        :param length: reference length - number of bases
        """
        self.name = name
        self.start = start
        self.length = length
