from ..common import is_placed_alignment
from ..fasta_index import FastaIndex


class BaiBamIterator:
    def __init__(self, bam_file, start_index=None, end_index=None):
        """
        :param bam_file: pysam.AlignmentFile
        :param start_index: iterate from 0-based index inclusive
        :param end_index: iterate to 0-based index inclusive
        
        Iterates over indexed reads only. Placed unmapped reads are included.
        """
        if start_index is end_index is not None:
            assert start_index <= end_index
        
        if not bam_file.has_index():
            raise IndexError('BAM has no index')
        
        self.bam_file = bam_file
        self.fai = FastaIndex(bam_file)
        # empty iterator
        self.iterator = iter(())
        
        self.start_ref_id = None
        self.curr_ref_id = None
        self.end_ref_id = None
        self.counter = 0
        
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
                    end=self.end_ref_pos + 1  # add one for inclusive range
                )
            else:
                # multiple iterators
                self.start_ref_id = self.fai.ref_id(self.start_ref_name)
                self.curr_ref_id = self.start_ref_id
                self.end_ref_id = self.fai.ref_id(self.end_ref_name)
                
                if self.curr_ref_id > self.end_ref_id:
                    raise ValueError("Start reference has position after end reference.")
    
    def __next_iterator(self):
        if self.curr_ref_id is None or self.curr_ref_id > self.end_ref_id:
            # end iteration
            iterator = None
        else:
            if self.curr_ref_id == self.start_ref_id:
                # first iterator
                iterator = self.bam_file.fetch(
                    reference=self.start_ref_name,
                    start=self.start_ref_pos,
                    end=None
                )
            elif self.curr_ref_id == self.end_ref_id:
                # last iterator
                iterator = self.bam_file.fetch(
                    reference=self.end_ref_name,
                    start=None,
                    end=self.end_ref_pos + 1  # half open interval
                )
            else:
                # intermediate iterator
                iterator = self.bam_file.fetch(
                    reference=self.fai.ref_name(self.curr_ref_id)
                )
            
            self.curr_ref_id += 1
        
        return iterator
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            alignment = next(self.iterator)
            self.counter += 1
            return alignment
        except StopIteration:
            self.iterator = self.__next_iterator()
            if self.iterator is None:
                return None
            else:
                return next(self)


class UnmappedOnlyBamIterator:
    def __init__(self, bam_file):
        """
        :param bam_file: pysam.AlignmentFile
        
        Iterates only on both placed and unplaced unmapped reads.
        """
        self.iterator = bam_file.fetch(until_eof=True)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            alignment = next(self.iterator)
            # proceed to (placed or unplaced) unmapped alignment if mapped
            while not alignment.is_unmapped:
                alignment = next(self.iterator)
            
            return alignment
        except StopIteration:
            return None


class FullBamIterator:
    def __init__(self, bam_file):
        """
        Iterates over all reads in BAM file.
        """
        self.iterator = bam_file.fetch(until_eof=True)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            return next(self.iterator)
        except StopIteration:
            return None


class UnmappedBamIterator:
    def __init__(self, bam_file, start_index, end_index):
        """
        :param bam_file: pysam.AlignmentFile
        :param start_index: iterate from 0-based index inclusive
        :param end_index: iterate to 0-based index inclusive
        
        Iterator includes all unplaced unmapped reads. BAM index ensures that BAM is sorted
        and is needed to resolve range.
        Iterator assumes that unplaced alignment can be anywhere in BAM file.
        """
        assert start_index <= end_index
        
        if not bam_file.has_index():
            raise IndexError('BAM has no index')
        
        self.fai = FastaIndex(bam_file)
        self.start_ref_name, self.start_ref_pos = self.fai.resolve_start_pos(start_index)
        self.end_ref_name, self.end_ref_pos = self.fai.resolve_end_pos(end_index)
        
        self.start_ref_id = self.fai.ref_id(self.start_ref_name)
        self.end_ref_id = self.fai.ref_id(self.end_ref_name)
        
        self.iterator = bam_file.fetch(until_eof=True)
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            while True:
                alignment = next(self.iterator)
                if alignment.is_unmapped:
                    return alignment
                
                curr_ref_id = self.fai.ref_id(alignment.reference_name)
                if curr_ref_id < self.start_ref_id or (
                                curr_ref_id == self.start_ref_id and alignment.reference_end <= self.start_ref_pos):
                    # before range
                    continue
                
                if curr_ref_id > self.end_ref_id or (
                                curr_ref_id == self.end_ref_id and alignment.reference_start > self.end_ref_pos):
                    # after range
                    continue
                
                return alignment
        
        except StopIteration:
            return None
