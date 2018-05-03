import re


class NotFoundError(Exception):
    pass


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
    # TODO deprecated
    @classmethod
    def allele(cls, seq: str, ref_seq: str) -> str:
        """
        :param seq: allele sequence
        :param ref_seq:  reference allele sequence
        :return: expanded CIGAR of the allele
        """
        if len(seq) == len(ref_seq):
            cigar = cls.OP_MATCH * len(seq)
        else:
            mutual_len = min(len(seq), len(ref_seq))
            cigar = cls.OP_MATCH * mutual_len
            
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
    def matching_allele(
            cls,
            seq: str,
            exp_cigar: str,
            alleles: list,
            ref_allele: str,
            seq_pos: int,
            cigar_pos: int = None,
    ) -> (str, str):
        """
        :param seq:
        :param exp_cigar:
        :param alleles:
        :param ref_allele:
        :param seq_pos:
        :param cigar_pos:
        :return:
        """
        assert ref_allele in alleles
        if cigar_pos is None:
            cigar_pos = cls.seq_pos2cigar_pos(exp_cigar, seq_pos)
        
        for allele in reversed(sorted(alleles, key=len)):  # type: str
            
            allele_cigar = cls.allele(allele, ref_allele)
            mutual_len = min(len(allele_cigar), len(exp_cigar) - cigar_pos)
            if allele_cigar[:mutual_len] != exp_cigar[cigar_pos: cigar_pos + mutual_len]:
                continue
            
            if cigar_pos + len(allele_cigar) > len(exp_cigar):
                break
            
            is_seq_match = seq[seq_pos: seq_pos + len(allele)] == allele
            if not is_seq_match:
                continue
            
            is_cigar_match = exp_cigar[cigar_pos:cigar_pos + len(allele_cigar)] == allele_cigar
            if is_seq_match and is_cigar_match:
                return allele, allele_cigar
        
        raise NotFoundError
