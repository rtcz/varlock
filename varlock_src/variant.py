from array import array

import pysam

from varlock_src.cigar import Cigar


# TODO WTF1
def qual_str2array(qual_str: str):
    return array('B', map(lambda x: ord(x) - 33, qual_str))


# TODO WTF2
def qual_array2str(qual_array: array):
    if qual_array is None:
        return ''
    return ''.join(map(lambda x: chr(x + 33), qual_array))


class AlignedVariant:
    def __init__(
            self,
            alignment: pysam.AlignedSegment,
            pos: int = None,
            end_pos: int = None,
            ref_seq: str = None,
            is_mutated: bool = False
    ):
        """
        :param alignment:
        :param pos: position of variation in AlignedSegment.query_sequence in alternative sequence
        :param end_pos: position one base after INDEL variation in AlignedSegment.query_sequence in alternative sequence
        :param ref_seq: reference sequence of the variant
        :param is_mutated: True if alignment is mutated already
        """
        # TODO check if pos is a cigar match
        
        self.alignment = alignment
        self._pos = pos
        self._end_pos = end_pos
        self._is_snv = False
        self._is_indel = False
        self._ref_seq = ref_seq
        self._is_mutated = is_mutated
        
        if self._pos is not None:
            if self._end_pos is not None:
                assert self._pos < self._end_pos
                assert ref_seq is not None
                self._is_indel = True
            else:
                self._is_snv = True
                self._end_pos = self._pos + 1
    
    def is_present(self):
        return self._is_snv or self._is_indel
    
    @property
    def is_mutated(self):
        return self._is_mutated
    
    # TODO try to remove, pos and end_pos should be immutable (so as variant presence)
    def clear(self):
        self._is_snv = False
        self._is_indel = False
        self._pos = None
        self._end_pos = None
    
    @property
    def seq(self):
        """
        Get variant sequence.
        """
        if self.is_present():
            return self.alignment.query_sequence[self._pos:self._end_pos]
        else:
            return None
    
    @seq.setter
    def seq(self, seq):
        """
        Set variant sequence.
        :param seq: new sequence
        """
        assert len(seq) > 0
        if seq == self.seq:
            return
        
        # if self.alignment.query_name == 'ERR015528.13775373' and self.alignment.reference_start == 1421565:
        #     print('xxx')
        #     print(self.alignment.query_name)
        #     print(self.alignment.reference_name)
        #     print(self.alignment.reference_start)
        #     print(self.alignment.query_sequence)
        #     print(self.alignment.cigarstring)
        #     print('pos %d' % self._pos)
        #     print('end_pos %d' % self._end_pos)
        #     print('reference allele %s' % self._ref_seq)
        #     print('alternative allele %s' % self.seq)
        #     print('masking allele %s' % seq)
        #     print('xxx')
        
        self._is_mutated = True
        
        # TODO
        old_sequence = self.alignment.query_sequence
        old_sequence_len = len(self.alignment.query_sequence)
        old_cigar_len = self.alignment.infer_query_length()
        old_cigar_str = self.alignment.cigarstring
        
        # save quality
        # TODO update quality string, save deleted qualities
        # assigning to query_sequence removes query_qualities
        
        # TODO do something about MD string if present
        
        if self._is_snv:
            # TODO treat extended CIGAR opearations (=, X)
            # TODO at least assert that corresponding cigar letter is M
            # expecting only M OP now - it does not change with different SNV
            query_qualities = self.alignment.query_qualities
            mut_seq = self.alignment.query_sequence[:self._pos] + seq + self.alignment.query_sequence[self._end_pos:]
            self.alignment.query_sequence = mut_seq
            self.alignment.query_qualities = query_qualities
        
        elif self._is_indel:
            
            seq, exp_cigar = Cigar.replace_allele(
                seq=self.alignment.query_sequence,
                exp_cigar=Cigar.tuples2exp_str(self.alignment.cigartuples),
                alleles=None,  # TODO
                allele_id=0,
                ref_allele_id=0,
                seq_pos=self._pos
            )
            self.alignment.cigartuples = Cigar.exp_str2tuples(exp_cigar)
            self.alignment.query_sequence = seq
            
            # if len(self.alignment.query_sequence) == self.alignment.infer_query_length():
            #     print(self.alignment.reference_start)
            #     exit(0)
            
            # assert len(self.alignment.query_sequence) == self.alignment.infer_query_length()
            
            if len(self.alignment.query_sequence) != self.alignment.infer_query_length():
                print(self.alignment.query_name)
                print(self.alignment.reference_name)
                print(self.alignment.reference_start)
                
                print("pos %s" % self._pos)
                print("end_pos %s" % self._end_pos)
                print("masking allele %s" % seq)
                print("alternative allele %s" % old_sequence[self._pos:self._end_pos])
                print("reference allele %s" % self._ref_seq)
                
                print('OLD')
                print(old_sequence)
                print(old_sequence_len)
                print(old_cigar_str)
                print(old_cigar_len)
                
                print('NEW')
                print(self.alignment.query_sequence)
                print(len(self.alignment.query_sequence))
                print(self.alignment.cigarstring)
                print(self.alignment.infer_query_length())
                exit(0)
