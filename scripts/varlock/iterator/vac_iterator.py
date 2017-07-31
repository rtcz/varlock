from varlock.fasta_index import FastaIndex
from varlock.po import VacSnvRecord, VacIndelRecord
from varlock.vac import Vac


class VacIterator:
    """
    Iterates over VAC file alternating between both SNV and INDEL records in order of genomic index.
    """
    
    def __init__(self, vac_filename: str, fai: FastaIndex):
        with open(vac_filename, 'rb') as vac_file:
            # read header
            self.snv_count, self.indel_count = Vac.read_header(vac_file)
        
        self.snv_file = open(vac_filename, 'rb')
        self.snv_file.seek(Vac.HEADER_SIZE)
        
        self.indel_file = open(vac_filename, 'rb')
        self.indel_file.seek(Vac.HEADER_SIZE)
        self.indel_file.seek(self.snv_count * Vac.SNV_RECORD_SIZE)
        
        self.fai = fai
        self.counter = 0
        
        # init iteration
        self.__read_snv()
        self.__read_indel()
    
    def __read_snv(self):
        if self.snv_count > 0:
            self.snv_count -= 1
            snv_index, snv_list = Vac.read_snv_record(self.snv_file)
            ref_name, ref_pos = self.fai.index2pos(snv_index)
            self.snv_record = VacSnvRecord(
                index=snv_index,
                ref_name=ref_name,
                ref_pos=ref_pos,
                ac=snv_list
            )
        else:
            self.snv_record = None
            self.snv_file.close()
    
    def __read_indel(self):
        if self.indel_count > 0:
            self.indel_count -= 1
            indel_index, indel_map = Vac.read_snv_record(self.snv_file)
            ref_name, ref_pos = self.fai.index2pos(indel_index)
            self.indel_record = VacIndelRecord(
                index=indel_index,
                ref_name=ref_name,
                ref_pos=ref_pos,
                indel_map=indel_map
            )
        else:
            self.indel_record = None
            self.indel_file.close()
    
    def next_snv(self):
        self.counter += 1
        snv_record = self.snv_record
        self.__read_snv()
        return snv_record
    
    def next_indel(self):
        self.counter += 1
        indel_list = self.indel_record
        self.__read_indel()
        return indel_list
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.snv_record is not None and self.indel_record is not None:
            if self.snv_record.index == self.indel_record.index:
                raise ValueError("SNV and INDEL have the same position")
            elif self.snv_record.index < self.indel_record.index:
                return self.next_snv()
            else:
                return self.next_indel()
        
        elif self.snv_record is not None:
            # INDELs depleted
            return self.next_snv()
        elif self.indel_record is not None:
            # SNVs depleted
            return self.next_indel()
        else:
            # EOF
            return None
