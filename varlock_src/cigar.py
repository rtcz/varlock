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
    
    OPS = [OP_MATCH, OP_INS, OP_DEL, OP_REF_SKIP, OP_SOFT_CLIP, OP_HARD_CLIP, OP_PAD, OP_EQUAL, OP_DIFF]
    
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
    
    @classmethod
    def _merge_append(cls, cigar: list, op_type: int, op_length: int):
        """
        Merging new OP with last OP from cigar
        when OPs are of the same type.
        Otherwise new OP is appended to cigar.
        :param cigar: [<op_type, op_length>, ...]
        """
        if len(cigar):
            prev_type, prev_length = cigar[-1]
            if prev_type == op_type:
                # merge OPs
                cigar[-1] = (op_type, op_length + prev_length)
            else:
                cigar.append((op_type, op_length))
        else:
            cigar.append((op_type, op_length))
    
    @staticmethod
    def str2tuples(cigar_str: str) -> list:
        result = []
        for count, op in re.findall(re.compile('([0-9]+)([A-Z]+)'), cigar_str):
            result.append((Cigar.OPS.index(op), int(count)))
        
        return result
    
    @staticmethod
    def tuples2str(cigar_tuples: list) -> str:
        return ''.join([str(count) + Cigar.OPS[int(op)] for op, count in cigar_tuples])
    
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
    def shrink(cigar_str: str):
        """
        :param cigar_str: expanded CIGAR string
        :return:
        """
        result = ''
        prev_op = None
        op_count = 0
        for op in cigar_str:
            if prev_op is None or op == prev_op:
                op_count += 1
            else:
                result += str(op_count) + prev_op
                op_count = 1
            
            prev_op = op
        
        if op_count > 0:
            result += str(op_count) + prev_op
        
        return result
    
    @classmethod
    def mask(cls, cigar: str, pos: int, seq: str, ref_seq: str):
        """
        :param cigar: CIGAR string of alternative sequence
        :param pos: masking position within sequence described by cigar
        :param seq: masking sequence
        :param ref_seq: reference of masking sequence
        :return:
        """
        alt_pos = 0
        ref_pos = 0
        exp_cigar = cls.expand(cigar)
        i = 0
        
        # moving on both reference and alternative sequences until pos
        while alt_pos < pos:
            op = exp_cigar[i]
            if op in [cls.OP_MATCH, cls.OP_EQUAL, cls.OP_DIFF]:
                # consumes both
                alt_pos += 1
                ref_pos += 1
            
            elif op in [cls.OP_INS, cls.OP_SOFT_CLIP]:
                # consumes alternative sequence
                alt_pos += 1
            
            elif op in [cls.OP_DEL, cls.OP_REF_SKIP]:
                # consumes reference sequence
                ref_pos += 1
            
            elif op in [cls.OP_HARD_CLIP, cls.OP_PAD]:
                # consumes neither sequence
                pass
            else:
                raise ValueError('Unknwon CIGAR operation %s' % op)
            
            i += 1
        
        assert pos == alt_pos
        
        new_cigar = exp_cigar[:i]
        
        # moving only on reference sequence
        max_shared_len = min(len(seq), len(ref_seq))
        new_cigar += cls.OP_MATCH * max_shared_len
        target_ref_pos = ref_pos + max_shared_len
        
        len_diff = len(seq) - len(ref_seq)
        if len_diff > 0:
            new_cigar += cls.OP_INS * len_diff
        elif len_diff < 0:
            new_cigar += cls.OP_DEL * -len_diff
            target_ref_pos += -len_diff
        
        while ref_pos <= target_ref_pos:
            op = exp_cigar[i]
            if op in [cls.OP_MATCH, cls.OP_EQUAL, cls.OP_DIFF]:
                # consumes both
                alt_pos += 1
                ref_pos += 1
            
            elif op in [cls.OP_INS, cls.OP_SOFT_CLIP]:
                # consumes alternative sequence
                alt_pos += 1
            
            elif op in [cls.OP_DEL, cls.OP_REF_SKIP]:
                # consumes reference sequence
                ref_pos += 1
            
            elif op in [cls.OP_HARD_CLIP, cls.OP_PAD]:
                # consumes neither sequence
                pass
            else:
                raise ValueError('Unknwon CIGAR operation %s' % op)
            
            i += 1
        
        new_cigar += exp_cigar[(i - 1):]
        
        # for i in range(max(len(seq), len(ref_seq))):
        # if i < len(seq) and i < len(ref_seq):
        #     alt_base = seq[i]
        #     ref_base = ref_seq[i]
        #     if alt_base == ref_base:
        #         if mismatch_count == 1:
        #             # single base mismatch
        #             new_cigar += cls.OP_MATCH
        #         elif mismatch_count > 1:
        #             # multiple previous mismatches
        #             new_cigar +=
        #
        #         if mismatch_count > 0:
        #             new_cigar
        #             mismatch_count = 0
        #
        #         new_cigar += cls.OP_MATCH
        #     else:
        #         mismatch_count += 1
        
        return cls.shrink(new_cigar)
    
    @classmethod
    def del_subrange(cls, cigar: list, start_pos: int, end_pos: int):
        # TODO refactor ???
        """
        Removes range from CIGAR tuples, starting with start_pos
        and ending right before end_pos.
        CIGAR tuples inside the range are deleted and lengths
        of intersecting tuples are lowered respectively.
        :param cigar: [(op_type, op_length), ...]
        :param start_pos: position in AlignedSegment.query_sequence
        :param end_pos: position in AlignedSegment.query_sequence
        :param ref_range_len: length of subrange on reference sequence
        :return: new CIGAR tuples
        """
        assert end_pos > start_pos >= 0
        new_cigar = []
        curr_pos = 0
        started = None
        ended = None
        normalize = False
        for i in range(len(cigar)):
            curr_type, curr_length = cigar[i]
            prev_pos = curr_pos
            if curr_type not in [cls.OP_DEL, cls.OP_HARD_CLIP, cls.OP_REF_SKIP]:
                curr_pos += curr_length
            
            if started is None and start_pos < curr_pos:
                started = i
            
            if ended is None and end_pos <= curr_pos:
                ended = i
            
            if started is None:
                # current OP is before the OP that contains start_pos
                new_cigar.append((curr_type, curr_length))
            elif ended is None:
                # started but not ended
                if started == i:
                    # current OP contains start_pos only
                    new_length = start_pos - prev_pos
                    if new_length:
                        new_cigar.append((curr_type, new_length))
                else:
                    # current OP is omitted
                    pass
            elif started == i and ended == i:
                # current OP contains both start_pos and end_pos
                new_length = (start_pos - prev_pos) + (curr_pos - end_pos)
                if new_length:
                    new_cigar.append((curr_type, new_length))
                else:
                    # current OP is omitted
                    normalize = True
            elif ended == i:
                # current OP contains end_pos only
                new_length = curr_pos - end_pos
                if new_length:
                    cls._merge_append(new_cigar, curr_type, new_length)
                else:
                    # current OP is omitted
                    normalize = True
            elif curr_type == cls.OP_DEL and curr_pos == end_pos:
                # deletion OP after the OP, skip it
                normalize = True
            elif normalize:
                normalize = False
                cls._merge_append(new_cigar, curr_type, curr_length)
            else:
                new_cigar.append((curr_type, curr_length))
        
        return new_cigar
    
    @classmethod
    def place_op(cls, cigar: list, op_pos: int, op_type: int, op_length: int):
        """
        Place CIGAR operation into existing CIGAR string (in form of list of tuples).
        Function does not support placing hard and soft clip operations
        since they can appear only at CIGAR ends. Moreover HARD clipped bases
        are not part of AlignedSegment.query_sequence.
        :param cigar: [<op_type, op_length>, ...]
        :param op_pos: position in AlignedSegment.query_sequence
        :param op_type: CIGAR operation type id
        :param op_length: CIGAR operation length
        :return: new CIGAR tuples
        """
        # TODO OP_PAD
        assert op_type in [cls.OP_MATCH_ID, cls.OP_EQUAL_ID, cls.OP_DIFF_ID, cls.OP_INS_ID, cls.OP_DEL_ID,
                           cls.OP_REF_SKIP_ID]
        assert op_length > 0
        assert op_pos >= 0
        curr_pos = 0
        new_cigar = []
        
        for i in range(len(cigar)):
            curr_type, curr_length = cigar[i]
            # OP_HARD_CLIP is not contained in cigar_tuples nor query_sequence
            prev_pos = curr_pos
            if curr_type in [cls.OP_DEL_ID, cls.OP_REF_SKIP_ID]:
                new_cigar.append((curr_type, curr_length))
            else:
                curr_pos += curr_length
                if prev_pos == op_pos:
                    # the OP is in place of current OP
                    if op_type == curr_type:
                        # merge OPs
                        new_cigar.append((curr_type, curr_length + op_length))
                    else:
                        # place the OP and shift current OP
                        cls._merge_append(new_cigar, op_type, op_length)
                        new_cigar.append((curr_type, curr_length))
                
                elif prev_pos < op_pos < curr_pos:
                    # place the OP inside current OP
                    if op_type == curr_type:
                        # merge the OP
                        new_cigar.append((curr_type, curr_length + op_length))
                    else:
                        new_cigar.append((curr_type, op_pos - prev_pos))
                        new_cigar.append((op_type, op_length))
                        new_cigar.append((curr_type, curr_pos - op_pos))
                else:
                    # continue
                    new_cigar.append((curr_type, curr_length))
        
        if op_pos == curr_pos:
            # place the OP at the end of CIGAR
            cls._merge_append(new_cigar, op_type, op_length)
        
        return new_cigar
    
    @classmethod
    def compute(cls, ref_seq: str, alt_seq: str):
        """
        Computes variant CIGAR with respect to the reference.
        Cigar can contain only insert, delete and match operations.
        :param ref_seq:
        :param alt_seq:
        :return: [<op_type, op_length>, ...]
        """
        assert len(alt_seq) > 0 or len(ref_seq) > 0
        cigar = []
        match_len = ins_len = del_len = 0
        
        if ref_seq == alt_seq:
            # reference variant
            cigar.append((Cigar.OP_MATCH_ID, len(alt_seq)))
        else:
            for i in range(max(len(ref_seq), len(alt_seq))):
                match_op = ins_op = del_op = False
                if i >= len(ref_seq):
                    # alt_seq is longer
                    ins_op = True
                
                elif i >= len(alt_seq):
                    # ref_seq is longer
                    del_op = True
                
                else:
                    # TODO treat [X, =] OP cases
                    if ref_seq[i] == alt_seq[i]:
                        match_op = True
                    else:
                        del_op = True
                        ins_op = True
                
                match_len += match_op
                ins_len += ins_op
                del_len += del_op
                
                if match_op:
                    if del_len == ins_len == 1:
                        # treat 1D1I as a match
                        cls._merge_append(cigar, cls.OP_MATCH_ID, 1)
                    else:
                        if ins_len:
                            cigar.append((cls.OP_INS_ID, ins_len))
                        if del_len:
                            cigar.append((cls.OP_DEL_ID, del_len))
                    
                    ins_len = del_len = 0
                
                if ins_op or del_op:
                    if match_len:
                        cls._merge_append(cigar, cls.OP_MATCH_ID, match_len)
                    
                    match_len = 0
            
            # finalize cigar
            if match_len:
                cls._merge_append(cigar, cls.OP_MATCH_ID, match_len)
            elif del_len == ins_len == 1:
                cls._merge_append(cigar, cls.OP_MATCH_ID, 1)
            else:
                if del_len:
                    cigar.append((cls.OP_DEL_ID, del_len))
                if ins_len:
                    cigar.append((cls.OP_INS_ID, ins_len))
        
        return cigar
