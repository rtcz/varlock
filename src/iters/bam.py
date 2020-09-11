import pysam

from src.bam import open_bam
from src.fasta_index import FastaIndex


class BamIterator:
    def __iter__(self):
        return self
    
    def __next__(self) -> pysam.AlignedSegment:
        # noinspection PyTypeChecker
        return None


def bam_iterator(
        filename: iter,
        start_index: int,
        end_index: int,
        unmapped_only: bool,
        include_unmapped: bool
):
    """
    Factory for selection of proper iterator. FullBamIterator is excluded.
    :param filename: BAM filename
    :param start_index:
    :param end_index:
    :param unmapped_only:
    :param include_unmapped:
    :return:
    """
    if unmapped_only:
        bam_iter = UnmappedBamIterator(filename)
    elif include_unmapped:
        bam_iter = RangedBamIterator(filename, start_index, end_index)
    else:
        # mapped only
        bam_iter = MappedBamIterator(filename, start_index, end_index)
    return bam_iter


class MappedBamIterator(BamIterator):
    def __init__(self, bam_filename: str, start_index: int = None, end_index: int = None):
        """
        :param bam_filename:
        :param start_index: iterate from 0-based index inclusive
        :param end_index: iterate to 0-based index inclusive
        Iterates over mapped reads only.
        """
        if start_index is not None and end_index is not None:
            assert start_index <= end_index
        
        self._bam_file = open_bam(bam_filename, 'rb')
        
        if not self._bam_file.has_index():
            raise IndexError('BAM has no index')
        
        self._fai = FastaIndex.from_bam(self._bam_file)
        # empty iterator
        self._iterator = iter(())
        
        self.start_ref_id = None
        self.curr_ref_id = None
        self.end_ref_id = None
        self.counter = 0

        if start_index is end_index is None:
            # fetch all
            self._iterator = self._bam_file.fetch()
        else:
            self.start_ref_name, self.start_ref_pos = self._fai.resolve_start_pos(start_index)
            self.end_ref_name, self.end_ref_pos = self._fai.resolve_end_pos(end_index)
            
            if self.start_ref_name == self.end_ref_name:
                # single iterator
                self._iterator = self._bam_file.fetch(
                    reference=self.start_ref_name,
                    start=self.start_ref_pos,
                    end=self.end_ref_pos
                )
            else:
                # multiple iterators
                self.start_ref_id = self._fai.ref_id(self.start_ref_name)
                self.curr_ref_id = self.start_ref_id
                self.end_ref_id = self._fai.ref_id(self.end_ref_name)
                
                if self.curr_ref_id > self.end_ref_id:
                    raise ValueError("Start reference has position after end reference.")
    
    def __next_iterator(self):
        if self.curr_ref_id is None or self.curr_ref_id > self.end_ref_id:
            # end iteration
            iterator = None
        else:
            if self.curr_ref_id == self.start_ref_id:
                # first iterator
                iterator = self._bam_file.fetch(
                    reference=self.start_ref_name,
                    start=self.start_ref_pos,
                    end=None
                )
            elif self.curr_ref_id == self.end_ref_id:
                # last iterator
                iterator = self._bam_file.fetch(
                    reference=self.end_ref_name,
                    start=None,
                    end=self.end_ref_pos + 1  # half open interval
                )
            else:
                # intermediate iterator
                iterator = self._bam_file.fetch(
                    reference=self._fai.ref_name(self.curr_ref_id)
                )
            
            self.curr_ref_id += 1
        
        return iterator
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._bam_file.close()
    
    def __next__(self) -> pysam.AlignedSegment:
        try:
            alignment = next(self._iterator)
            while alignment.is_unmapped:
                alignment = next(self._iterator)
            self.counter += 1
            return alignment
        except StopIteration:
            self._iterator = self.__next_iterator()
            if self._iterator is None:
                self._bam_file.close()
                # noinspection PyTypeChecker
                return None
            else:
                return next(self)


class UnmappedBamIterator(BamIterator):
    def __init__(self, bam_filename: str):
        """
        :param bam_filename:
        Iterates only on both placed and unplaced unmapped reads.
        
        Iterator does not use BAI.
        """
        self._bam_file = open_bam(bam_filename, 'rb')
        self._iterator = self._bam_file.fetch(until_eof=True)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._bam_file.close()
    
    def __iter__(self):
        return self
    
    def __next__(self) -> pysam.AlignedSegment:
        try:
            alignment = next(self._iterator)
            # proceed to (placed or unplaced) unmapped alignment if mapped
            while not alignment.is_unmapped:
                alignment = next(self._iterator)
            
            return alignment
        except StopIteration:
            self._bam_file.close()
            # noinspection PyTypeChecker
            return None


class FullBamIterator(BamIterator):
    def __init__(self, bam_filename: str):
        """
        :param bam_filename:
        Iterates over all reads within BAM file.
        
        Iterator does not use BAI.
        """
        self._bam_file = open_bam(bam_filename, 'rb')
        self._iterator = self._bam_file.fetch(until_eof=True)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._bam_file.close()
    
    def __next__(self) -> pysam.AlignedSegment:
        try:
            return next(self._iterator)
        except StopIteration:
            self._bam_file.close()
            # noinspection PyTypeChecker
            return None


class RangedBamIterator(BamIterator):
    def __init__(self, bam_filename: str, start_index: int, end_index: int):
        """
        :param bam_filename:
        :param start_index: iterate from 0-based index inclusive
        :param end_index: iterate to 0-based index inclusive
        
        Iterator includes mapped reads within the range and all unmapped reads.
        BAM index ensures that BAM is sorted and is needed to resolve range.
        Iterator assumes that unplaced alignment can be anywhere in BAM file.
        """
        assert start_index <= end_index
        
        self._bam_file = open_bam(bam_filename, 'rb')
        
        if not self._bam_file.has_index():
            raise IndexError('BAM has no index')
        
        self._fai = FastaIndex.from_bam(self._bam_file)
        self._start_ref_name, self._start_ref_pos = self._fai.resolve_start_pos(start_index)
        self._end_ref_name, self._end_ref_pos = self._fai.resolve_end_pos(end_index)
        
        self._start_ref_id = self._fai.ref_id(self._start_ref_name)
        self._end_ref_id = self._fai.ref_id(self._end_ref_name)
        
        self._iterator = self._bam_file.fetch(until_eof=True)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._bam_file.close()
    
    def __next__(self) -> pysam.AlignedSegment:
        try:
            while True:
                alignment = next(self._iterator)
                if alignment.is_unmapped:
                    return alignment
                
                curr_ref_id = self._fai.ref_id(alignment.reference_name)
                if curr_ref_id < self._start_ref_id or (
                                curr_ref_id == self._start_ref_id and alignment.reference_end <= self._start_ref_pos):
                    # before range
                    continue
                
                if curr_ref_id > self._end_ref_id or (
                                curr_ref_id == self._end_ref_id and alignment.reference_start >= self._end_ref_pos):
                    # after range
                    continue

                return alignment
        
        except StopIteration:
            self._bam_file.close()
            # noinspection PyTypeChecker
            return None
