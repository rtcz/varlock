from varlock.po import VacRecord
from varlock.vac import Vac


class VacIterator:
    def __init__(self, vac_file, fai):
        self.vac_file = vac_file
        self.fai = fai
        self.counter = 0
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            index, ac_tuple = Vac.read_snv_record(self.vac_file)
        except EOFError:
            return None
        
        ref_name, ref_pos = self.fai.index2pos(index)
        self.counter += 1
        
        return VacRecord(
            index=index,
            ref_name=ref_name,
            ref_pos=ref_pos,
            ac=ac_tuple
        )
