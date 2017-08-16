from varlock.common import BASES
from varlock.po import DiffSnvRecord, DiffIndelRecord

from varlock.fasta_index import FastaIndex
import varlock.bdiff as bdiff


class BdiffIterator:
    def __init__(
            self,
            bdiff_io: bdiff.BdiffIO,
            fai: FastaIndex,
            start_index: int,
            end_index: int
    ):
        assert bdiff_io.is_read_mode
        self._fai = fai
        self._counter = 0
        self._bdiff_io = bdiff_io
        
        self._start_pos, self._end_pos = self._bdiff_io.tell_range(start_index, end_index)
        self._bdiff_io.seek(self._start_pos)
    
    @property
    def counter(self):
        return self._counter
    
    def __iter__(self):
        return self
    
    def __next__(self):
        """
        :return: next DIFF record (SNV or INDEL) by genomic index order
        SNV diff record
            mut_map.key: mutated base
            mut_map.value: original base
        """
        record = None
        if self._bdiff_io.tell() <= self._end_pos:
            index, alts = self._bdiff_io.read_record()
            ref_name, ref_pos = self._fai.index2pos(index)
            if isinstance(alts, tuple):
                # is SNV
                record = DiffSnvRecord(
                    index,
                    ref_name,
                    ref_pos,
                    mut_map=dict(zip(alts, BASES))
                )
            else:
                # is INDEL
                record = DiffIndelRecord(
                    index,
                    ref_name,
                    ref_pos,
                    alts
                )
        
        self._counter += 1
        return record