from ..common import BASES
from ..diff import Diff
from ..po import DiffRecord


class DiffIterator:
    def __init__(self, diff_file, fai, start_index=None, end_index=None):
        assert start_index <= end_index
        self.diff_file = diff_file
        self.fai = fai
        self.start_index = start_index
        self.end_index = end_index
        self.counter = 0
        
        # TODO
        
        if start_index is not None:
            Diff.seek_pos(diff_file, start_index)
            
        if end_index is not None:
            end_pos = Diff.seek_pos(diff_file, end_index)
            diff_file.truncate(end_pos)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            index, mut_tuple = Diff.read_record(self.diff_file)
        except EOFError:
            return None
        
        ref_name, ref_pos = self.fai.index2pos(index)
        self.counter += 1
        
        return DiffRecord(
            index=index,
            ref_name=ref_name,
            ref_pos=ref_pos,
            mut_map=dict(zip(mut_tuple, BASES))
        )
    
    def counter(self):
        return self.counter
