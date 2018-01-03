import pysam

from varlock.cigar import Cigar


class VariantPosition(object):
    def __init__(self, index: int, ref_name: str, ref_pos: int):
        self.index = index
        self.ref_name = ref_name
        self.ref_pos = ref_pos
    
    def __str__(self):
        return '#%d %s:%d' % (self.index, self.ref_name, self.ref_pos)


# TODO rename VariantDiff
class DiffRecord(VariantPosition):
    def __init__(self, index: int, ref_name: str, ref_pos: int, mut_map: dict, ref_seq: str):
        super().__init__(index, ref_name, ref_pos)
        self.mut_map = mut_map
        self.ref_seq = ref_seq


# TODO rename SnvDiff
class DiffSnvRecord(DiffRecord):
    # def __init__(self, index: int, ref_name: str, ref_pos: int, mut_map: dict, ref_base: str):
    #     super().__init__(index, ref_name, ref_pos)
    #     self.mut_map = mut_map
    #     self._ref_base = ref_base
    #
    # def ref_seq(self):
    #     return self._ref_base
    pass


# TODO rename IndelDiff
class DiffIndelRecord(DiffRecord):
    # def __init__(self, index: int, ref_name: str, ref_pos: int, mut_map: dict, ref_seq: str):
    #     super().__init__(index, ref_name, ref_pos)
    #     self.mut_map = mut_map
    #     self._ref_seq = ref_seq
    #
    # def ref_seq(self):
    #     return self._ref_seq
    pass


# TODO rename VariantOccurence
class VacRecord(VariantPosition):
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


# TODO rename SnvOccurence
class VacSnvRecord(VacRecord):
    # def __init__(self, index: int, ref_name: str, ref_pos: int, freqs: tuple, ref_id: int):
    #     super().__init__(index, ref_name, ref_pos)
    #     self.freqs = freqs
    #     self.ref_id = ref_id
    #
    # def ref_seq(self):
    #     return self.freqs[self.ref_id]
    pass


# TODO rename IndelOccurence
class VacIndelRecord(VacRecord):
    # """
    #     :param index: genomic index
    #     :param ref_name:
    #     :param ref_pos:
    #     :param freqs: alleles frequencies, the first item is reference frequency
    #     :param seqs: alleles sequences, the first item is reference sequence
    #     """
    #     super().__init__(index, ref_name, ref_pos)
    #     self.freqs = freqs
    #     self.seqs = seqs
    #     self._ref_id = ref_id
    #
    # def ref_freq(self):
    #     return self.freqs[self._ref_id]
    #
    # def ref_seq(self):
    #     return self.seqs[self._ref_id]
    pass


# TODO separate file, tests ?
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
        :param seq: new sequence
        """
        mut_seq = self.alignment.query_sequence[:self._pos]
        mut_seq += seq
        mut_seq += self.alignment.query_sequence[self._end_pos:]
        self.alignment.query_sequence = mut_seq
        
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
            #     # VCF indels should always have single matched base preffix
            #     tmp_end_pos = self._end_pos + 1
            # else:
            #     tmp_end_pos = self._end_pos
            
            # print()
            # print(self.alignment.reference_start + 1)
            # print(self.alignment.query_alignment_sequence)
            # print(self.alignment.cigarstring)
            # print(self.alignment.cigartuples)
            # print(self._pos, self._end_pos)
            
            tmp_cigar = Cigar.del_subrange(
                self.alignment.cigartuples,
                self._pos,
                # tmp_end_pos
                self._end_pos
            )
            
            # print(tmp_cigar)
            # print(self._ref_seq)
            # print(seq)
            
            variant_cigar = Cigar.variant(self._ref_seq, seq)
            
            # print(variant_cigar)
            
            pos = self._pos
            for tpl in variant_cigar:
                tmp_cigar = Cigar.place_op(
                    tmp_cigar,
                    pos,
                    tpl[0],
                    tpl[1]
                )
                if tpl[0] not in [Cigar.OP_DEL, Cigar.OP_REF_SKIP]:
                    pos += tpl[1]
            
            # print(tmp_cigar)
            
            self.alignment.cigartuples = tmp_cigar
    
    @staticmethod
    def _first_bases_match(seqs: list):
        assert len(seqs) > 0
        for i in range(len(seqs) - 1):
            if seqs[i] != seqs[i + 1]:
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
