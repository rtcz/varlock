from ..diff import Diff
from ..po import DiffSnvRecord
from ..common import BASES


# TODO refactor - use _filename as parameter

class DiffIterator:
    def __init__(self, diff_file, fai, start_index, end_index):
        assert start_index < end_index
        
        self.diff_file = Diff.slice(diff_file, start_index, end_index)
        self.diff_file.seek(Diff.HEADER_SIZE)
        self.fai = fai
        self.start_index = start_index
        self.end_index = end_index
        self.counter = 0
    
    def __iter__(self):
        return self
    
    def __next__(self):
        """
        :return: diff record with reversed mutation map
        mut_map.key: mutated base
        mut_map.value: original base
        """
        try:
            index, mut_tuple = Diff.read_record(self.diff_file)
        except EOFError:
            self.diff_file.close()
            return None
        
        ref_name, ref_pos = self.fai.index2pos(index)
        self.counter += 1
        
        return DiffSnvRecord(
            index=index,
            ref_name=ref_name,
            ref_pos=ref_pos,
            mut_map=dict(zip(mut_tuple, BASES))
        )
    
    def __enter__(self):
        return self
    
    def __exit__(self):
        self.diff_file.close()
