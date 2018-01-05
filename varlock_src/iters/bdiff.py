from varlock_src.po import SnvDiff, IndelDiff, VariantDiff

from varlock_src.fasta_index import FastaIndex
import varlock_src.bdiff as bdiff


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
    
    def __next__(self) -> VariantDiff:
        """
        :return: next DIFF record (SNV or INDEL) by genomic index order
        BDIFF reversed mut_map:
            key: mutated variation
            value: original variation
        """
        record = None
        if self._bdiff_io.tell() <= self._end_pos:
            # index, alts = self._bdiff_io.read_record()
            index, is_indel, ref_seq, mut_map = self._bdiff_io.read_record()
            ref_name, ref_pos = self._fai.index2pos(index)
            
            if is_indel:
                # it is assumed that INDEL record was saved as:
                # sorted(alts) -> permuted(alts)
                record = IndelDiff(index, ref_name, ref_pos, mut_map, ref_seq)
            else:
                # is SNV
                # it is assumed that SNV record was saved as:
                # BASES -> permuted(BASES)
                record = SnvDiff(index, ref_name, ref_pos, mut_map, ref_seq)
        
        self._counter += 1
        return record