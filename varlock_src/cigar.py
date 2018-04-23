import re


class Cigar:
    OP_MATCH = 'M'  # match (does not have to be equal to the refrence)
    OP_INS = 'I'  # insertion
    OP_DEL = 'D'  # deletion
    OP_REF_SKIP = 'N'  # intron - similar to D
    OP_SOFT_CLIP = 'S'  # soft clipped bases (not aligned, present in read)
    OP_HARD_CLIP = 'H'  # hard clipped bases (not aligned, not present in read)
    OP_PAD = 'P'  # padding / gap
    OP_EQUAL = '='  # equal match
    OP_DIFF = 'X'  # different match
    
    ID2OP = {
        0: OP_MATCH,
        1: OP_INS,
        2: OP_DEL,
        3: OP_REF_SKIP,
        4: OP_SOFT_CLIP,
        5: OP_HARD_CLIP,
        6: OP_PAD,
        7: OP_EQUAL,
        8: OP_DIFF
    }
    
    OP2ID = {
        OP_MATCH: 0,
        OP_INS: 1,
        OP_DEL: 2,
        OP_REF_SKIP: 3,
        OP_SOFT_CLIP: 4,
        OP_HARD_CLIP: 5,
        OP_PAD: 6,
        OP_EQUAL: 7,
        OP_DIFF: 8
    }
    
    # TODO remove
    OP_MATCH_ID = 0  # M,
    OP_INS_ID = 1  # I,
    OP_DEL_ID = 2  # D,
    OP_REF_SKIP_ID = 3  # N,
    OP_SOFT_CLIP_ID = 4  # S,
    OP_HARD_CLIP_ID = 5  # H,
    OP_PAD_ID = 6  # P,
    OP_EQUAL_ID = 7  # =, type of M
    OP_DIFF_ID = 8  # X, type of M
    
    # @staticmethod
    # def str2tuples(cigar_str: str) -> list:
    #     result = []
    #     for count, op in re.findall(re.compile('([0-9]+)([A-Z]+)'), cigar_str):
    #         result.append((Cigar.OPS.index(op), int(count)))
    #
    #     return result
    #
    # @staticmethod
    # def tuples2str(cigar_tuples: list) -> str:
    #     return ''.join([str(count) + Cigar.OPS[int(op)] for op, count in cigar_tuples])
    
    @staticmethod
    def expand(cigar_str: str) -> str:
        """
        :param cigar_str: standard CIGAR string
        :return:
        """
        result = ''
        for count, op in re.findall(re.compile('([0-9]+)([A-Z]+)'), cigar_str):
            result += op * int(count)
        
        return result
    
    @staticmethod
    def shrink(exp_cigar: str):
        """
        :param exp_cigar: expanded CIGAR string
        :return:
        """
        result = ''
        prev_op = None
        op_count = 0
        for op in exp_cigar:
            if prev_op is None or op == prev_op:
                op_count += 1
            else:
                result += str(op_count) + prev_op
                op_count = 1
            
            prev_op = op
        
        if op_count > 0:
            result += str(op_count) + prev_op
        
        return result
    
    @staticmethod
    def tuples2exp_str(cigar_tuples: list) -> str:
        """ll
        :param cigar_tuples:
        :return: expanded CIGAR string
        """
        return ''.join([Cigar.ID2OP[op] * count for op, count in cigar_tuples])
    
    @staticmethod
    def exp_str2tuples(exp_cigar: str) -> list:
        """
        :param exp_cigar: expanded CIGAR string
        :return:
        """
        result = []
        prev_op = None
        op_count = 0
        for op in exp_cigar:
            if prev_op is None or op == prev_op:
                op_count += 1
            else:
                result.append((Cigar.OP2ID[prev_op], op_count))
                op_count = 1
            
            prev_op = op
        
        if op_count > 0:
            result.append((Cigar.OP2ID[prev_op], op_count))
        
        return result
    
    @classmethod
    def _consumes_alt(cls, cigar_op: str):
        """
        :param cigar_op: CIGAR operation on single position
        :return:
        """
        return cigar_op in [cls.OP_MATCH, cls.OP_EQUAL, cls.OP_DIFF, cls.OP_INS, cls.OP_SOFT_CLIP]
    
    @classmethod
    def _consumes_ref(cls, cigar_op: str):
        """
        :param cigar_op: CIGAR operation on single position
        :return:
        """
        return cigar_op in [cls.OP_MATCH, cls.OP_EQUAL, cls.OP_DIFF, cls.OP_DEL, cls.OP_REF_SKIP]
    
    # TODO extended CIGAR operations (=,X)
    @classmethod
    def variant(cls, seq: str, ref_seq: str) -> str:
        """
        :param seq:
        :param ref_seq:
        :return: expanded CIGAR of variant
        """
        if len(seq) == len(ref_seq):
            cigar = cls.OP_MATCH * len(seq)
        else:
            max_shared_len = min(len(seq), len(ref_seq))
            cigar = cls.OP_MATCH * max_shared_len
            
            len_diff = len(seq) - len(ref_seq)
            if len_diff > 0:
                cigar += cls.OP_INS * len_diff
            elif len_diff < 0:
                cigar += cls.OP_DEL * -len_diff
        
        return cigar
    
    @classmethod
    def ref_len(cls, exp_cigar: str) -> int:
        result = 0
        for op in exp_cigar:
            result += cls._consumes_ref(op)
        
        return result
    
    @classmethod
    def alt_len(cls, exp_cigar: str) -> int:
        result = 0
        for op in exp_cigar:
            result += cls._consumes_alt(op)
        
        return result

    # TODO deprecated
    @classmethod
    def mask(cls, exp_cigar: str, seq_pos: int, alt_seq: str, ref_seq: str) -> str:
        """
        :param exp_cigar: expanded CIGAR string of alternative sequence
        :param seq_pos: masking position within or right after alternative sequence associated with cigar
        :param alt_seq: allele sequence
        :param ref_seq: reference allele sequence
        :return:
        """
        # if
        cigar_pos = cls.seq_pos2cigar_pos(exp_cigar, seq_pos)
        if cigar_pos is None:
            raise ValueError('sequence position %d not found in CIGAR %s' % (seq_pos, exp_cigar))
        
        pre_cigar = exp_cigar[:cigar_pos]
        ref_pos = cls.ref_len(pre_cigar)
        masking_cigar = cls.variant(alt_seq, ref_seq)
        ending_ref_pos = ref_pos + cls.ref_len(masking_cigar)
        
        post_cigar = ''
        for i in range(cigar_pos, len(exp_cigar)):
            if ref_pos > ending_ref_pos:
                break
            elif ref_pos == ending_ref_pos:
                # move further on cigar if ref_seq was not consumed
                post_cigar = exp_cigar[i:]
                break
            
            ref_pos += cls._consumes_ref(exp_cigar[i])
        
        return pre_cigar + masking_cigar + post_cigar
    
    @classmethod
    def ref_pos2cigar_pos(cls, exp_cigar: str, ref_pos: int):
        pass
    
    @classmethod
    def seq_pos2cigar_pos(cls, exp_cigar: str, seq_pos: int) -> int:
        """
        :param exp_cigar: expanded CIGAR string
        :param seq_pos: alternative sequence position
        :return: CIGAR position of seq_pos
        """
        result = None
        alt_pos = -1
        
        for i in range(len(exp_cigar)):
            alt_pos += cls._consumes_alt(exp_cigar[i])
            if alt_pos == seq_pos:
                result = i
                break
        
        return result
    
    @classmethod
    def _matching_allele(
            cls,
            seq: str,
            exp_cigar: str,
            alleles: list,
            ref_allele_id: int,
            seq_pos: int,
            cigar_pos: int = None,
    ) -> (str, str):
        """
        :param seq:
        :param exp_cigar:
        :param alleles:
        :param ref_allele_id:
        :param seq_pos:
        :param cigar_pos:
        :return:
        """
        if cigar_pos is None:
            cigar_pos = cls.seq_pos2cigar_pos(exp_cigar, seq_pos)
        
        for allele in reversed(sorted(alleles, key=len)):  # type: str
            if seq_pos + len(allele) > len(seq):
                # allele exceeds sequence
                # TODO treat partialy matched alleles
                break
            
            is_seq_match = seq[seq_pos: seq_pos + len(allele)] == allele
            if not is_seq_match:
                continue
            
            allele_cigar = cls.variant(allele, alleles[ref_allele_id])
            assert cigar_pos + len(allele_cigar) <= len(exp_cigar)
            
            # if cigar_pos + len(allele_cigar) > len(exp_cigar):
            #     breakNerr
            
            is_cigar_match = exp_cigar[cigar_pos:cigar_pos + len(allele_cigar)] == allele_cigar
            if is_seq_match and is_cigar_match:
                return allele, allele_cigar
        
        # TODO
        raise ValueError
    
    # TODO merge match_allele() and mask() methods
    @classmethod
    def replace_allele(
            cls,
            seq: str,
            exp_cigar: str,
            alleles: list,
            allele_id: int,
            ref_allele_id: int,
            seq_pos: int
    ) -> (str, str):
        """
        Replace allele present in the sequence with another allele.
        Allele is replaced only if some allele from the list is present in the sequence.
        :param seq: sequence
        :param exp_cigar: expanded CIGAR string of the sequence
        :param seq_pos: position of alleles within sequence
        :param alleles: list of known alleles for the position
        :param allele_id: ID of desired allele within the list
        :param ref_allele_id: ID of reference allele within the list
        :return:tuple of sequence and exp_cigar
        """
        assert 0 <= allele_id < len(alleles)
        
        cigar_pos = cls.seq_pos2cigar_pos(exp_cigar, seq_pos)
        if cigar_pos is None:
            raise ValueError('sequence position %d not found in CIGAR %s' % (seq_pos, exp_cigar))
        
        matched_allele, matched_allele_cigar = cls._matching_allele(
            seq,
            exp_cigar,
            alleles,
            ref_allele_id,
            seq_pos,
            cigar_pos
        )
        
        if alleles[allele_id] != matched_allele:
            # allele was found and is different from desired allele
            
            # update sequence
            result_seq = seq[:seq_pos] + alleles[allele_id] + seq[seq_pos + len(matched_allele):]
            
            # update CIGAR string
            allele_cigar = cls.variant(alleles[allele_id], alleles[ref_allele_id])
            result_cigar = exp_cigar[:cigar_pos] + allele_cigar + exp_cigar[cigar_pos + len(matched_allele_cigar):]
            
            return result_seq, result_cigar
        else:
            return seq, exp_cigar
    
    # TODO deprecated
    @classmethod
    def matching_alleles(cls, exp_cigar: str, allele_pos: int, alleles: list, ref_allele_id: int) -> list:
        """
        :param exp_cigar:
        :param allele_pos:
        :param alleles:
        :param ref_allele_id:
        :return: list of CIGAR matched alleles with the biggest possible length
        """
        assert 0 <= ref_allele_id < len(alleles)
        matches = []
        
        cigar_pos = cls.seq_pos2cigar_pos(exp_cigar, allele_pos)
        if cigar_pos is None:
            raise ValueError('sequence position %d not found in CIGAR %s' % (allele_pos, exp_cigar))
        
        for allele in reversed(sorted(alleles, key=len)):  # type: str
            if len(matches) and len(matches[-1]) > len(allele):
                # previously matched allele is longer
                break
            
            allele_cigar = cls.variant(allele, alleles[ref_allele_id])
            if cigar_pos + len(allele_cigar) > len(exp_cigar):
                # allele exceedes CIGAR end
                # TODO treat partialy matched alleles
                break
            
            if exp_cigar[cigar_pos:cigar_pos + len(allele_cigar)] == allele_cigar:
                matches.append(allele)
        
        return matches
