# TODO refactor, use @property, etc
from enum import Enum


class VariantType(Enum):
    SNV = 1
    INDEL = 2


class GenomicPosition(object):
    def __init__(self, index: int, ref_name: str, ref_pos: int):
        self._index = index
        self._ref_name = ref_name
        self._ref_pos = ref_pos
    
    @property
    def index(self) -> int:
        return self._index
    
    @property
    def ref_name(self) -> str:
        return self._ref_name
    
    @property
    def ref_pos(self) -> int:
        return self._ref_pos
    
    def __str__(self):
        return '#%d %s:%d' % (self.index, self.ref_name, self.ref_pos + 1)


class Variant:
    def __init__(
            self,
            position: GenomicPosition,
            vtype: VariantType,
            alleles: list,
            ref_allele: str
    ):
        # assert ref_allele in alleles
        self._position = position
        self._alleles = alleles
        self._ref_allele = ref_allele
        self._vtype = vtype
    
    @property
    def pos(self):
        return self._position
    
    def is_type(self, vtype: VariantType) -> bool:
        return self._vtype == vtype
    
    @property
    def alleles(self) -> list:
        return self._alleles
    
    @property
    def ref_allele(self):
        return self._ref_allele
    
    def __str__(self):
        return '%s %s' % (self._position, self._vtype)


class VariantDiff(Variant):
    def __init__(
            self,
            position: GenomicPosition,
            vtype: VariantType,
            mut_map: dict,
            ref_allele: str
    ):
        super().__init__(
            position,
            vtype,
            list(mut_map.values()),
            ref_allele
        )
        # assert ref_allele in mut_map.values()
        self._mut_map = mut_map
    
    @property
    def mut_map(self) -> dict:
        return self._mut_map


class VariantOccurrence(Variant):
    def __init__(
            self,
            position: GenomicPosition,
            vtype: VariantType,
            freqs: list,
            alleles: list,
            ref_allele: str
    ):
        super().__init__(
            position,
            vtype,
            alleles,
            ref_allele
        )
        assert len(freqs) == len(alleles)
        self._freqs = freqs
    
    @property
    def freqs(self) -> list:
        return self._freqs


class FastaSequence:
    def __init__(self, index: int, name: str, start: int, length: int):
        """
        :param index: reference id
        :param name: reference name
        :param start: 0-based first position
        :param length: number of bases
        """
        self.id = index
        self.name = name
        self.start = start
        self.length = length
    
    def __str__(self):
        return '#%d %s %d %d' % (self.id, self.name, self.start, self.length)
