from ..common import BASES
from ..diff import Diff
from ..po import DiffRecord


class DiffIterator:
    def __init__(self, diff_file, fai, start_index, end_index):
        assert start_index < end_index
        
        self.diff_file = Diff.slice(diff_file, start_index, end_index)
        self.diff_file.seek(Diff.HEADER_LENGTH)
        self.fai = fai
        self.start_index = start_index
        self.end_index = end_index
        self.counter = 0
    
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
