from typing import Optional

import pysam

from varlock_src.cigar import Cigar, NotFoundError
from varlock_src.po import Variant, VariantType


class AlleleAlignment:
    def __init__(
            self,
            alignment: pysam.AlignedSegment,
            variant: Variant = None,
            is_mutated: bool = False
    ):
        """
        :param alignment:
        :param variant:
        :param is_mutated: True if alignment is mutated already
        """
        self._alignment = alignment
        self._variant = variant
        self._is_mutated = is_mutated
        
        if variant is not None:
            assert alignment.reference_name == variant.pos.ref_name
            if not alignment.is_unmapped and alignment.reference_start <= variant.pos.ref_pos < alignment.reference_end:
                # alignment is mapped and is covering variant
                self._seq_pos = self._ref_pos2seq_pos(alignment, variant.pos.ref_pos)
                if self._seq_pos is not None:
                    # find matching alleles with respect to CIGAR string
                    self._exp_cigar = Cigar.tuples2exp_str(alignment.cigartuples)
                    self._cigar_pos = Cigar.seq_pos2cigar_pos(self._exp_cigar, self._seq_pos)
                    assert self._cigar_pos is not None
                    
                    self._allele, self._allele_cigar = self._find_allele()
    
    def _find_allele(self) -> (str, str):
        allele = None
        allele_cigar = None
        
        if self.is_snv:
            allele = self._alignment.query_sequence[self._seq_pos]
            allele_cigar = self._exp_cigar[self._cigar_pos]
            # TODO consider present CIGAR and =/X operations
            assert allele_cigar == Cigar.OP_MATCH
        
        elif self.is_indel:
            try:
                # if self._alignment.query_name == 'ERR015528.7173436' and self._alignment.reference_start == 1664944:
                #     print()
                #     print(self._alignment.query_sequence)
                #     print('alleles %s' % str(self._variant.alleles))
                #     # print('actual allele %s' % self._allele)
                #     print('reference allele %s' % self._variant.ref_allele)
                #     print(self._exp_cigar)
                #     print(self._seq_pos)
                #     print(self._cigar_pos)
                #     print(self._variant)
                #     print()
                
                allele, allele_cigar = Cigar.matching_allele(
                    seq=self._alignment.query_sequence,
                    exp_cigar=self._exp_cigar,
                    alleles=self._variant.alleles,
                    ref_allele=self._variant.ref_allele,
                    seq_pos=self._seq_pos,
                    cigar_pos=self._cigar_pos
                )
                if self._cigar_pos + len(allele_cigar) >= len(self._exp_cigar):
                    # Only with CIGAR after an allele it is guranteed
                    # that replacement allele can be found too.
                    allele, allele_cigar = None, None
            
            except NotFoundError:
                pass
        
        return allele, allele_cigar
    
    @staticmethod
    def _ref_pos2seq_pos(alignment: pysam.AlignedSegment, ref_pos: int) -> int:
        """
        Retrieve base position in sequence string at refence position.
        Alignment and ref_pos are assumed to be of the same reference.
        :param alignment: pysam.AlignedSegment
        :param ref_pos: reference position of base
        :return: AlignedSegment.query_sequence position matched to ref_pos.
        None is returned if matching position is not found.
        """
        # TODO optimalize: (try matches_only=True)
        # TODO optimalize: case when alignment is full matched based on CIGAR (e.g. 30M)
        
        seq_pos = None
        for current_seq_pos, current_ref_pos in alignment.get_aligned_pairs(matches_only=False, with_seq=False):
            # search for base in snv position
            if current_ref_pos == ref_pos:
                seq_pos = current_seq_pos
                break
        
        return seq_pos
    
    @property
    def variant(self) -> Variant:
        return self._variant
    
    @property
    def alignment(self) -> pysam.AlignedSegment:
        return self._alignment
    
    @property
    def is_known(self) -> bool:
        try:
            return self._allele is not None
        except AttributeError:
            return False
    
    @property
    def is_mutated(self) -> bool:
        return self._is_mutated
    
    @property
    def is_snv(self) -> bool:
        return self._variant.is_type(VariantType.SNV)
    
    @property
    def is_indel(self) -> bool:
        return self._variant.is_type(VariantType.INDEL)
    
    @property
    def allele(self) -> Optional[str]:
        """
        Get variant sequence.
        """
        if self.is_known:
            return self._allele
        else:
            return None
    
    @allele.setter
    def allele(self, allele: str):
        """
        Set variant sequence.
        :param allele: new allele sequence
        """
        assert len(allele) > 0
        if allele not in self._variant.alleles:
            ValueError('allele %s is not in available list %s' % (allele, self._variant.alleles))
        
        if not self.is_known:
            return
        
        if allele == self.allele:
            return
        
        self._is_mutated = True
        
        # if self.alignment.query_name == 'ERR015528.25489857' and self.alignment.reference_start == 1272714:
        #     print()
        #     print(self.alignment.query_name)
        #     print(self.alignment.reference_name)
        #     print(self.alignment.reference_start)
        #     print(self.alignment.query_sequence)
        #     print(self._exp_cigar)
        #     print('seq_pos %d' % self._seq_pos)
        #     print('reference allele %s' % self._variant.ref_allele)
        #     print('actual allele %s' % self._allele)
        #     print('masking allele %s' % allele)
        #     print('alleles %s' % self._variant.alleles)
        #     print('variant %s' % self._variant)
        #     print()
        
        # TODO do something about MD string if present
        
        if self.is_snv:
            assert len(allele) == 1
            # TODO treat extended CIGAR opearations (=, X)
            # expecting only M operation in place of SNV
            query_qualities = self._alignment.query_qualities
            pre_seq = self._alignment.query_sequence[:self._seq_pos]
            post_seq = self._alignment.query_sequence[self._seq_pos + 1:]
            self._alignment.query_sequence = pre_seq + allele + post_seq
            self._alignment.query_qualities = query_qualities
        
        elif self.is_indel:
            # save quality
            # TODO update quality string, save deleted qualities
            # assigning to query_sequence removes query_qualities
            
            seq, exp_cigar = self._replace_indel(allele, self._allele, self._allele_cigar)
            self._alignment.cigartuples = Cigar.exp_str2tuples(exp_cigar)
            self._alignment.query_sequence = seq
            
            assert len(self._alignment.query_sequence) == self._alignment.infer_query_length()
            
            # if self.alignment.query_name == 'ERR015528.25489857' and self.alignment.reference_start == 1272714:
            #     print()
            #     print(seq)
            #     print(exp_cigar)
            #     print()
    
    def _replace_indel(self, allele: str, matched_allele: str, matched_allele_cigar: str) -> (str, str):
        # TODO extended CIGAR operations (=,X)
        
        # update sequence
        pre_seq = self._alignment.query_sequence[:self._seq_pos]
        post_seq = self._alignment.query_sequence[self._seq_pos + len(matched_allele):]
        result_seq = pre_seq + allele + post_seq
        
        # update CIGAR string
        allele_cigar = Cigar.allele(
            allele,
            self._variant.ref_allele
        )
        pre_cigar = self._exp_cigar[:self._cigar_pos]
        post_cigar = self._exp_cigar[self._cigar_pos + len(matched_allele_cigar):]
        result_cigar = pre_cigar + allele_cigar + post_cigar
        
        return result_seq, result_cigar
