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
        self._iter = VacFileIterator(vac_filename, fai)
        self._iter_rnd = RandomSnvIterator(fai, mut_p, rnd)
        
        self._last_record = None
        
        self._next_record = next(self._iter)  # type: po.VariantOccurrence
        self._next_record_rnd = next(self._iter_rnd)  # type: po.VariantOccurrence
    
    @property
    def counter(self):
        return self._iter.counter + self._iter_rnd.counter
    
    def __iter__(self):
        return self
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._iter.__exit__(exc_type, exc_val, exc_tb)
    
    def _next(self) -> po.VariantOccurrence:
        curr_record = self._next_record
        self._next_record = next(self._iter)
        return curr_record
    
    def _next_rnd(self) -> po.VariantOccurrence:
        curr_record = self._next_record_rnd
        self._next_record_rnd = next(self._iter_rnd)
        return curr_record
    
    def __next__(self) -> po.VariantOccurrence:
        if self._next_record is None and self._next_record_rnd is None:
            # noinspection PyTypeChecker
            return None
        elif self._next_record_rnd is None:
            return self._next()
        elif self._next_record is None:
            return self._next_rnd()
        elif self._next_record.pos.index <= self._next_record_rnd.pos.index:
            # known variants have precedence
            return self._next()
        else:
            return self._next_rnd()
    
    # TODO this method is used as temporary solution for overlapping variants
    def next_after(self) -> po.VariantOccurrence:
        """
        :return: variant after reference span of previous variant
        """
        while True:
            record = self.__next__()
            if record is None or self._last_record is None:
                break
            
            if self._last_record.pos.ref_name != record.pos.ref_name:
                break
            
            ref_end = self._last_record.pos.ref_pos + len(self._last_record.ref_allele)
            if record.pos.ref_pos >= ref_end:
                break
        
        self._last_record = record
        return record


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
        self._snv_variant = self._read_snv()
        self._indel_variant = self._read_indel()
    
    @property
    def counter(self):
        return self._counter
    
    def _read_snv(self) -> po.VariantOccurrence:
        if self._snv_count > 0:
            self._snv_count -= 1
            index, ref_id, freqs = Vac.read_snv_record(self._snv_file)
            ref_name, ref_pos = self._fai.index2pos(index)
            
            variant = po.VariantOccurrence(
                position=po.GenomicPosition(index, ref_name, ref_pos),
                vtype=po.VariantType.SNV,
                freqs=freqs,
                alleles=BASES,
                ref_allele=BASES[ref_id]
            )
        else:
            variant = None
            self._snv_file.close()
        
        return variant
    
    def _read_indel(self) -> po.VariantOccurrence:
        if self._indel_count > 0:
            self._indel_count -= 1
            index, counts, seqs = Vac.read_indel_record(self._indel_file)
            ref_name, ref_pos = self._fai.index2pos(index)
            
            variant = po.VariantOccurrence(
                position=po.GenomicPosition(index, ref_name, ref_pos),
                vtype=po.VariantType.INDEL,
                freqs=counts,
                alleles=seqs,
                ref_allele=seqs[0]
            )
        
        else:
            variant = None
            self._indel_file.close()
        
        return variant
    
    def next_snv(self):
        self._counter += 1
        snv_variant = self._snv_variant
        self._snv_variant = self._read_snv()
        return snv_variant
    
    def next_indel(self):
        self._counter += 1
        indel_variant = self._indel_variant
        self._indel_variant = self._read_indel()
        return indel_variant
    
    def __iter__(self):
        return self
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._snv_file.close()
        self._indel_file.close()
    
    def __next__(self) -> po.VariantOccurrence:
        if self._snv_variant is not None and self._indel_variant is not None:
            if self._snv_variant.pos.index == self._indel_variant.pos.index:
                raise ValueError("SNV and INDEL have the same position (%s)" % str(self._snv_variant.pos.index))
            elif self._snv_variant.pos.index < self._indel_variant.pos.index:
                return self.next_snv()
            else:
                return self.next_indel()
        
        elif self._snv_variant is not None:
            # INDELs depleted
            self._indel_file.close()
            return self.next_snv()
        elif self._indel_variant is not None:
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
        
        return po.VariantOccurrence(
            position=po.GenomicPosition(index, ref_name, ref_pos),
            vtype=po.VariantType.SNV,
            freqs=[3, 3, 2, 2],  # approximate GC content in human genome,
            alleles=BASES,
            ref_allele=BASES[0]  # TODO use real value based on genome fasta via pyfaidx
        )
