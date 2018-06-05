import varlock_src.bdiff as bdiff
from varlock_src import po
from varlock_src.fasta_index import FastaIndex


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
        
        self._next_record = self._next()  # type: po.VariantDiff
    
    @property
    def counter(self):
        return self._counter
    
    def __iter__(self):
        return self
    
    def __next__(self) -> po.VariantDiff:
        """
        :return: next DIFF record (SNV or INDEL) by genomic index order
        BDIFF reversed mut_map:
            key: mutated variation
            value: original variation
        """
        # record = None
        # if self._bdiff_io.tell() <= self._end_pos:
        #     # index, alts = self._bdiff_io.read_record()
        #     index, is_indel, ref_seq, mut_map = self._bdiff_io.read_record()
        #     ref_name, ref_pos = self._fai.index2pos(index)
        #
        #     record = po.VariantDiff(
        #         position=po.GenomicPosition(index, ref_name, ref_pos),
        #         vtype=po.VariantType.INDEL if is_indel else po.VariantType.SNV,
        #         mut_map=mut_map,
        #         ref_allele=ref_seq
        #     )
        #
        #     self._counter += 1
        #
        # return record
        
        record = self._next_record
        if record is not None:
            self._next_record = self._next()
            self._counter += 1
        
        return record
    
    def _next(self) -> po.VariantDiff:
        record = None
        if self._bdiff_io.tell() <= self._end_pos:
            index, is_snv, ref_seq, mut_map = self._bdiff_io.read_record()
            ref_name, ref_pos = self._fai.index2pos(index)
            
            record = po.VariantDiff(
                position=po.GenomicPosition(index, ref_name, ref_pos),
                vtype=po.VariantType.SNV if is_snv else po.VariantType.INDEL,
                mut_map=mut_map,
                ref_allele=ref_seq
            )
        
        return record
    
    # TODO use this method to unmask overlapping variants (in reverse orded)
    def next_cluster(self) -> list:
        """
        Experimental method intended for unmasking overlapping variants
        :return: list of overlapping variants
        """
        cluster = []
        while True:
            record = self._next_record
            if record is None:
                break
            
            cluster.append(record)
            self._next_record = self._next()
            self._counter += 1
            
            if record.pos.ref_name != self._next_record.pos.ref_name:
                break
            
            ref_end = record.pos.ref_pos + len(record.ref_allele)
            if self._next_record.pos.ref_pos >= ref_end:
                break
        
        return cluster
