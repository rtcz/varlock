import numpy as np

import varlock_src.po as po
from varlock_src.common import BASES
from varlock_src.fasta_index import FastaIndex
from varlock_src.random import VeryRandom
from varlock_src.vac import Vac


# TODO do not use python iterator "interface", use custom next() and close() methods

class VariantIterator:
    def __init__(self, vac_filename: str, fai: FastaIndex, mut_p: float, rnd: VeryRandom):
        """
        Composite iterator using content of vac_file and randomness to return next VariantOccurence.
        Returned VariantOccurence is either read from vac_file or randomly generated.
        :param vac_filename:
        :param fai:
        :param mut_p: random variant (mutation) probability per genome base
        """
        self._iter1 = VacFileIterator(vac_filename, fai)
        self._iter2 = RandomSnvIterator(fai, mut_p, rnd)
        
        # initialize current values
        self._curr1 = next(self._iter1)
        self._curr2 = next(self._iter2)
    
    @property
    def counter(self):
        return self._iter1.counter + self._iter2.counter
    
    def __iter__(self):
        return self
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._iter1.__exit__(exc_type, exc_val, exc_tb)
    
    def _next1(self) -> po.VariantOccurrence:
        curr = self._curr1
        self._curr1 = next(self._iter1)
        return curr
    
    def _next2(self) -> po.VariantOccurrence:
        curr = self._curr2
        self._curr2 = next(self._iter2)
        return curr
    
    def __next__(self) -> po.VariantOccurrence:
        if self._curr1 is None and self._curr2 is None:
            # noinspection PyTypeChecker
            return None
        elif self._curr2 is None:
            return self._next1()
        elif self._curr1 is None:
            return self._next2()
        elif self._curr1.index <= self._curr2.index:
            # vac_file originated variants have precedence
            return self._next1()
        else:
            return self._next2()


class VacFileIterator:
    """
    Iterates over VAC file alternating between both SNV and INDEL records sorted by genomic index.
    """
    
    def __init__(self, vac_filename: str, fai: FastaIndex):
        """
        :param vac_filename:
        :param fai:
        does not affect bases listed in VAC file
        """
        with open(vac_filename, 'rb') as vac_file:
            # read header
            self._snv_count, self._indel_count = Vac.read_header(vac_file)
        
        self._snv_file = open(vac_filename, 'rb')
        self._snv_file.seek(Vac.HEADER_SIZE)
        
        self._indel_file = open(vac_filename, 'rb')
        self._indel_file.seek(Vac.HEADER_SIZE + self._snv_count * Vac.SNV_RECORD_SIZE)
        
        self._fai = fai
        self._counter = 0
        
        # init iteration
        self._read_snv()
        self._read_indel()
    
    @property
    def counter(self):
        return self._counter
    
    def _read_snv(self):
        if self._snv_count > 0:
            self._snv_count -= 1
            index, ref_id, freqs = Vac.read_snv_record(self._snv_file)
            ref_name, ref_pos = self._fai.index2pos(index)
            self._snv_record = po.SnvOccurrence(
                index=index,
                ref_name=ref_name,
                ref_pos=ref_pos,
                freqs=freqs,
                seqs=BASES,
                ref_id=ref_id
            )
        else:
            self._snv_record = None
            self._snv_file.close()
    
    def _read_indel(self):
        if self._indel_count > 0:
            self._indel_count -= 1
            index, counts, seqs = Vac.read_indel_record(self._indel_file)
            ref_name, ref_pos = self._fai.index2pos(index)
            self._indel_record = po.IndelOccurrence(
                index=index,
                ref_name=ref_name,
                ref_pos=ref_pos,
                freqs=counts,
                seqs=seqs,
                ref_id=0
            )
        else:
            self._indel_record = None
            self._indel_file.close()
    
    def next_snv(self):
        self._counter += 1
        snv_record = self._snv_record
        self._read_snv()
        return snv_record
    
    def next_indel(self):
        self._counter += 1
        indel_list = self._indel_record
        self._read_indel()
        return indel_list
    
    def __iter__(self):
        return self
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._snv_file.close()
        self._indel_file.close()
    
    def __next__(self) -> po.VariantOccurrence:
        if self._snv_record is not None and self._indel_record is not None:
            if self._snv_record.index == self._indel_record.index:
                raise ValueError("SNV and INDEL have the same position (%s)" % str(self._snv_record.index))
            elif self._snv_record.index < self._indel_record.index:
                return self.next_snv()
            else:
                return self.next_indel()
        
        elif self._snv_record is not None:
            # INDELs depleted
            self._indel_file.close()
            return self.next_snv()
        elif self._indel_record is not None:
            # SNVs depleted
            self._snv_file.close()
            return self.next_indel()
        else:
            # EOF
            self._snv_file.close()
            self._indel_file.close()
            # noinspection PyTypeChecker
            return None


class RandomSnvIterator:
    """
    Random rare SNV iterator. SNVs positions are randomly generated across whole genome.
    """
    
    def __init__(self, fai: FastaIndex, mut_p: float, rnd: VeryRandom):
        """
        :param fai:
        :param mut_p: random variant (mutation) probability per genome base
        :param rnd:
        """
        assert 0 <= mut_p <= 0.001
        self._fai = fai
        self._rnd = rnd
        length = fai.last_index() - fai.first_index()
        assert length >= 0
        mut_count = int(mut_p * length)
        
        # sample random genomic indices from uniform distribution
        self._indices = np.sort(rnd.rand_ints(fai.first_index(), fai.last_index(), mut_count))
        self._counter = 0
    
    @property
    def counter(self):
        return self._counter
    
    def __iter__(self):
        return self
    
    def __next__(self) -> po.VariantOccurrence:
        if self._counter >= len(self._indices):
            # noinspection PyTypeChecker
            return None
        
        index = self._indices[self._counter]
        self._counter += 1
        
        ref_name, ref_pos = self._fai.index2pos(index)
        return po.SnvOccurrence(
            index=index,
            ref_name=ref_name,
            ref_pos=ref_pos,
            freqs=[3, 3, 2, 2],  # approximate GC content in human genome
            seqs=BASES,
            # TODO use real value based on genome fasta via pyfaidx
            # this could be later utilized in CIGAR edit
            ref_id=0
        )
