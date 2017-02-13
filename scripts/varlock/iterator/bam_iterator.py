from ..fasta_index import FastaIndex


class BamIterator:
    def __init__(self, bam_file, start_index=None, end_index=None):
        """
        :param bam_file: pysam.AlignmentFile
        :param start_index:
        :param end_index:
        """
        if start_index is end_index is not None:
            assert start_index <= end_index
        
        self.bam_file = bam_file
        self.fai = FastaIndex(bam_file)
        # empty iterator
        self.iterator = iter(())
        
        self.start_ref_id = None
        self.curr_ref_id = None
        self.end_ref_id = None
        
        if start_index is end_index is None:
            # fetch all
            self.iterator = bam_file.fetch()
        else:
            self.start_ref_name, self.start_ref_pos = self.fai.resolve_start_pos(start_index)
            self.end_ref_name, self.end_ref_pos = self.fai.resolve_end_pos(end_index)
            
            if self.start_ref_name == self.end_ref_name:
                # single iterator
                self.iterator = bam_file.fetch(
                    reference=self.start_ref_name,
                    start=self.start_ref_pos,
                    end=self.end_ref_pos
                )
            else:
                # multiple iterators
                self.start_ref_id = self.fai.ref_id(self.start_ref_name)
                self.curr_ref_id = self.start_ref_id
                self.end_ref_id = self.fai.ref_id(self.end_ref_name)
                
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
                self.iterator = self.bam_file.fetch(
                    reference=self.start_ref_name,
                    start=self.start_ref_pos,
                    end=None
                )
            elif self.curr_ref_id == self.end_ref_id:
                # last iterator
                self.iterator = self.bam_file.fetch(
                    refenrece=self.end_ref_name,
                    start=None,
                    end=self.end_ref_pos
                )
            else:
                # intermediate iterator
                self.iterator = self.bam_file.fetch(
                    reference=self.fai.ref_name(self.curr_ref_id)
                )
                self.curr_ref_id += 1
        return next(self)
