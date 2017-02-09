from ..fasta_index import FastaIndex


class BamIterator:
    def __init__(self, bam_file, start_index=None, end_index=None):
        """
        :param bam_file: pysam.AlignmentFile
        :param start_index:
        :param end_index:
        """
        assert start_index <= end_index
        
        self.bam_file = bam_file
        self.fai = FastaIndex(bam_file)
        # empty iterator
        self.iterator = iter(())
        
        if end_index is not None:
            end_ref_name, end_ref_pos = self.fai.index2pos(end_index)
        
        self.start_pos = start_ref_pos
        self.end_pos = end_ref_pos
        
        self.start_ref_id = None
        self.curr_ref_id = None
        self.end_ref_id = None
        
        if start_ref_name == end_ref_name:
            # single iterator
            # fetch all if None or fetch single reference
            self.iterator = bam_file.fetch(start_ref_name, start_ref_pos, end_ref_pos)
        else:
            # multiple iterators
            self.start_ref_id = self.__ref2id(start_ref_name)
            self.curr_ref_id = self.start_ref_id
            self.end_ref_id = self.__ref2id(end_ref_name)
            
            if self.curr_ref_id > self.end_ref_id:
                raise ValueError("Start reference has position after end reference.")
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            return next(self.iterator)
        except StopIteration:
            if self.curr_ref_id is None or self.curr_ref_id > self.end_ref_id:
                # end iteration
                return None
            elif self.curr_ref_id == self.start_ref_id:
                # first iterator
                start_ref_name = self.fai_list[self.start_ref_id].name
                self.iterator = self.bam_file.fetch(start_ref_name, self.start_pos, None)
            elif self.curr_ref_id == self.end_ref_id:
                # last iterator
                end_ref_name = self.fai_list[self.end_ref_id].name
                self.iterator = self.bam_file.fetch(end_ref_name, None, self.end_pos)
            else:
                # intermediate iterator
                self.iterator = self.bam_file.fetch(self.fai_list[self.curr_ref_id].name)
            
            self.curr_ref_id += 1
            return self.__next__()
