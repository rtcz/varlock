import typing
from array import array
from typing import Optional

import numpy as np
import pysam

from src.cigar import Cigar, NotFoundError
from src.common import BASES
from src.po import Variant, VariantType


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
        self._allele = None
        self._allele_cigar = None
        
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
                    
                    if self.is_snv:
                        self._allele, self._allele_cigar = self._find_snv_allele()
                    elif self.is_indel:
                        self._allele, self._allele_cigar = self._find_indel_allele()
    
    def _find_indel_allele(self) -> (str, str):
        # if self._alignment.query_name == 'ERR013136.13215553' and self._alignment.reference_start == 10365830:
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
        try:
            allele, allele_cigar = Cigar.matching_allele(
                seq=self._alignment.query_sequence,
                exp_cigar=self._exp_cigar,
                alleles=self._variant.alleles,
                ref_allele=self._variant.ref_allele,
                seq_pos=self._seq_pos,
                cigar_pos=self._cigar_pos
            )
        except NotFoundError:
            return None, None
        
        if len(self._exp_cigar) > self._cigar_pos + len(allele_cigar):
            # Only with CIGAR operation after an allele it is guranteed
            # that replacement allele can be found too.
            return allele, allele_cigar
        else:
            return None, None
    
    def _find_snv_allele(self) -> (str, str):
        allele = self._alignment.query_sequence[self._seq_pos]
        allele_cigar = self._exp_cigar[self._cigar_pos]
        
        # TODO consider present CIGAR and =/X operations
        assert allele_cigar == Cigar.OP_MATCH
        
        if allele in BASES:
            return allele, allele_cigar
        else:
            return None, None
    
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
        return self._allele is not None
    
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
        return self._allele
    
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
        
        # if self._alignment.query_name == 'ERR015528.5973067' and self._alignment.reference_start == 7437170:
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
            # TODO do something about deleted qualities
            query_qualities = self._replace_qualities(
                allele_len=len(allele),
                matched_allele_len=len(self._allele),
                seq_pos=self._seq_pos,
                qualities=self._alignment.query_qualities
            )
            
            seq, exp_cigar = self._replace_indel(allele, self._allele, self._allele_cigar)
            self._alignment.cigartuples = Cigar.exp_str2tuples(exp_cigar)
            self._alignment.query_sequence = seq
            self._alignment.query_qualities = query_qualities
            assert len(self._alignment.query_sequence) == self._alignment.infer_query_length()
            
            # if self._alignment.query_name == 'ERR015528.5973067' and self._alignment.reference_start == 7437170:
            #     print()
            #     print(seq)
            #     print(exp_cigar)
            #     print()
        
        self._allele = allele
    
    @staticmethod
    def _replace_qualities(allele_len: int, matched_allele_len: int, seq_pos: int, qualities: array) -> array:
        result = array('B')
        if qualities is not None and len(qualities) > 0:
            len_diff = allele_len - matched_allele_len
            if len_diff > 0:
                # new allele is longer
                mean_qual = int(np.ceil(np.mean(qualities)))
                result = qualities[:seq_pos + matched_allele_len] + \
                         array('B', [mean_qual]) * len_diff + \
                         qualities[seq_pos + matched_allele_len:]
            elif len_diff < 0:
                # new allele is shorter
                result = qualities[:seq_pos + matched_allele_len + len_diff] + \
                         qualities[seq_pos + matched_allele_len:]
            else:
                # new allele has the same length
                result = qualities
        
        return result
    
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


def pileup_alleles(variants: typing.List[AlleleAlignment]):
    """
    :param variants: list of po.AlignedVariant
    :return: list of bases at specific position
    """
    result = []
    for variant in variants:
        if variant.is_known:
            result.append(variant.allele)
    
    return result
