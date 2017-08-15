from varlock import FastaIndex
from varlock.bdiff import BdiffIO
from varlock.po import DiffSnvRecord, DiffIndelRecord



class BdiffIterator:
    def __init__(self, filename, fai: FastaIndex, start_index, end_index):
        self._fai = fai
        self._counter = 0
        self._bdiff = BdiffIO(filename, 'r')
        
        self._read_snv()
        self._read_indel()
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._bdiff.close()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        """
        :return: next DIFF record by genomic index order
        """
        if self._snv_record is None and self._indel_record is None:
            raise StopIteration
        
        if self._indel_record is None:
            # INDELs are depleted
            return self._next_snv()
        
        if self._snv_record is None:
            # SNVs are depleted
            return self._next_indel()
        
        # both INDEL and SNV are present
        if self._snv_record.index == self._indel_record.index:
            raise ValueError("SNV and INDEL have the same position")
        elif self._snv_record.index < self._indel_record.index:
            return self._next_snv()
        else:
            return self._next_indel()
    
    def _read_snv(self):
        """
        Internally read SNV record
        """
        try:
            index, perm = self._bdiff.read_snv()
            ref_name, ref_pos = self._fai.index2pos(index)
            self._snv_record = DiffSnvRecord(
                index,
                ref_name,
                ref_pos,
                perm
            )
        except EOFError:
            self._snv_record = None
    
    def _next_snv(self):
        """
        Return current SNV and read the next one.
        :return: current SNV
        """
        snv_record = self._snv_record
        self._read_snv()
        return snv_record
    
    def _read_indel(self):
        """
        Internally read INDEL record.
        """
        try:
            index, alt_list = self._bdiff.read_indel()
            ref_name, ref_pos = self._fai.index2pos(index)
            self._indel_record = DiffIndelRecord(
                index,
                ref_name,
                ref_pos,
                alt_list
            )
        except EOFError:
            self._indel_record = None
    
    def _next_indel(self):
        """
        Return current INDEL and read the next one.
        :return: current INDEL
        """
        indel_record = self._indel_record
        self._read_indel()
        return indel_record
