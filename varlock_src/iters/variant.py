import numpy as np

import varlock_src.po as po
from varlock_src.common import BASES
from varlock_src.fasta_index import FastaIndex
from varlock_src.vac import Vac
from varlock_src.random import VeryRandom


class VariantIterator:
    def __init__(self, vac_filename: str, fai: FastaIndex, mut_p: int):
        """
        :param vac_filename:
        :param fai:
        :param mut_p:
        """
        pass
    
    def __next__(self) -> po.VariantOccurrence:
        # TODO
        pass


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
    
    def __init__(self, fai: FastaIndex, mut_p: int, rnd: VeryRandom):
        """
        :param mut_p: random variant (mutation) probability per genome base
        :param fai:
        """
        assert mut_p > 0
        self._fai = fai
        self._rnd = rnd
        length = fai.last_index() - fai.first_index()
        mut_count = int(mut_p * length)
        
        # sample random genomic indices from uniform distribution
        self._indices = np.sort(rnd.rand_ints(fai.first_index(), fai.last_index(), mut_count))
        self._counter = 0
    
    def __next__(self) -> po.VariantOccurrence:
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
