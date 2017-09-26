import varlock.po as po

from varlock.fasta_index import FastaIndex
from varlock.vac import Vac


class VacIterator:
    """
    Iterates over VAC file alternating between both SNV and INDEL records in order of genomic index.
    """
    
    def __init__(self, vac_filename: str, fai: FastaIndex):
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
            index, freqs = Vac.read_snv_record(self._snv_file)
            ref_name, ref_pos = self._fai.index2pos(index)
            self._snv_record = po.VacSnvRecord(
                index=index,
                ref_name=ref_name,
                ref_pos=ref_pos,
                freqs=freqs
            )
        else:
            self._snv_record = None
            self._snv_file.close()
    
    def _read_indel(self):
        if self._indel_count > 0:
            self._indel_count -= 1
            index, counts, seqs = Vac.read_indel_record(self._indel_file)
            ref_name, ref_pos = self._fai.index2pos(index)
            self._indel_record = po.VacIndelRecord(
                index=index,
                ref_name=ref_name,
                ref_pos=ref_pos,
                freqs=counts,
                seqs=seqs
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
    
    def __next__(self):
        if self._snv_record is not None and self._indel_record is not None:
            if self._snv_record.index == self._indel_record.index:
                raise ValueError("SNV and INDEL have the same position")
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
            return None
