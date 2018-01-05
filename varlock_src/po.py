class VariantPosition(object):
    def __init__(self, index: int, ref_name: str, ref_pos: int):
        self.index = index
        self.ref_name = ref_name
        self.ref_pos = ref_pos
    
    def __str__(self):
        return '#%d %s:%d' % (self.index, self.ref_name, self.ref_pos)


class VariantDiff(VariantPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, mut_map: dict, ref_seq: str):
        super().__init__(index, ref_name, ref_pos)
        self.mut_map = mut_map
        self.ref_seq = ref_seq


class SnvDiff(VariantDiff):
    pass


class IndelDiff(VariantDiff):
    pass


class VariantOccurrence(VariantPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, freqs: list, seqs: list, ref_id: int):
        super().__init__(index, ref_name, ref_pos)
        self.freqs = freqs
        self.seqs = seqs
        self.ref_id = ref_id
    
    @property
    def ref_freq(self):
        return self.freqs[self.ref_id]
    
    @property
    def ref_seq(self):
        return self.seqs[self.ref_id]

class SnvOccurrence(VariantOccurrence):

    def __str__(self):
        return '#%d %s:%d SNP' % (self.index, self.ref_name, self.ref_pos)


class IndelOccurrence(VariantOccurrence):

    def __str__(self):
        return '#%d %s:%d INDEL' % (self.index, self.ref_name, self.ref_pos)


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
