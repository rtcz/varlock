class Cigar:
    OP_MATCH = 0  # M, match (does not have to be equal to the refrence)
    OP_INS = 1  # I, insertion
    OP_DEL = 2  # D, deletion
    OP_REF_SKIP = 3  # N, intron - similar to D
    OP_SOFT_CLIP = 4  # S, soft masked bases (not aligned, present in read)
    OP_HARD_CLIP = 5  # H, hard masked bases (not aligned, not present in read)
    OP_PAD = 6  # P, padding / gap
    OP_EQUAL = 7  # =, type of M
    OP_DIFF = 8  # X, type of M
    
    # custom OPs, not part of CIGAR specification
    OP_DEL_INS = 9  # deletion and insertion at one place
    
    @classmethod
    def _safe_append(cls, cigar: list, op_type: int, op_length: int):
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
    
    @classmethod
    def del_subrange(cls, cigar: list, start_pos: int, end_pos: int):
        # TODO refactor ???
        """
        Removes range from CIGAR tuples, starting with start_pos
        and ending right before end_pos.
        CIGAR tuples inside the range are deleted and lengths
        of intersecting tuples are lowered respectively.
        :param cigar: [<op_type, op_length>, ...]
        :param start_pos: position in AlignedSegment.query_sequence
        :param end_pos: position in AlignedSegment.query_sequence
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
                    new_cigar.append((curr_type, start_pos - prev_pos))
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
                    cls._safe_append(new_cigar, curr_type, new_length)
                else:
                    # current OP is omitted
                    normalize = True
            elif curr_type == cls.OP_DEL and curr_pos == end_pos:
                # deletion OP after the OP, skip it
                normalize = True
            elif normalize:
                normalize = False
                cls._safe_append(new_cigar, curr_type, curr_length)
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
        assert op_type in [cls.OP_MATCH, cls.OP_EQUAL, cls.OP_DIFF, cls.OP_INS, cls.OP_DEL, cls.OP_REF_SKIP]
        assert op_length > 0
        assert op_pos >= 0
        curr_pos = 0
        new_cigar = []
        
        for i in range(len(cigar)):
            curr_type, curr_length = cigar[i]
            # OP_HARD_CLIP is not contained in cigar_tuples nor query_sequence
            prev_pos = curr_pos
            if curr_type not in [cls.OP_DEL, cls.OP_REF_SKIP]:
                curr_pos += curr_length
            
            if prev_pos == op_pos:
                # place the OP before current OP
                if op_type == curr_type:
                    # merge the OP
                    new_cigar.append((curr_type, curr_length + op_length))
                else:
                    cls._safe_append(new_cigar, op_type, op_length)
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
                new_cigar.append((curr_type, curr_length))
        
        if op_pos == curr_pos:
            # place the OP at the end of CIGAR
            cls._safe_append(new_cigar, op_type, op_length)
        
        return new_cigar
    
    @classmethod
    def variant(cls, ref_seq: str, var_seq: str):
        """
        Computes variant CIGAR with respect to the reference.
        Cigar can contain only insert, delete and match operations.
        :param ref_seq:
        :param var_seq:
        :return: [<op_type, op_length>, ...]
        """
        assert len(var_seq) > 0 or len(ref_seq) > 0
        cigar = []
        match_len = ins_len = del_len = 0
        
        if ref_seq == var_seq:
            # reference variant
            cigar.append((Cigar.OP_MATCH, len(var_seq)))
        else:
            for i in range(max(len(ref_seq), len(var_seq))):
                match_op = ins_op = del_op = False
                if i >= len(ref_seq):
                    # seq is longer
                    ins_op = True
                
                elif i >= len(var_seq):
                    # ref_seq is longer
                    del_op = True
                
                else:
                    # TODO treat [X, =] OP cases
                    if ref_seq[i] == var_seq[i]:
                        match_op = True
                    else:
                        del_op = True
                        ins_op = True
                
                match_len += match_op
                ins_len += ins_op
                del_len += del_op
                
                if match_op:
                    if del_len == ins_len == 1:
                        # treat as a match
                        cls._safe_append(cigar, cls.OP_MATCH, 1)
                    else:
                        if ins_len:
                            cigar.append((cls.OP_INS, ins_len))
                        if del_len:
                            cigar.append((cls.OP_DEL, del_len))
                    
                    ins_len = del_len = 0
                
                if ins_op or del_op:
                    if match_len:
                        cls._safe_append(cigar, cls.OP_MATCH, match_len)
                    
                    match_len = 0
            
            if match_len:
                cls._safe_append(cigar, cls.OP_MATCH, match_len)
            elif del_len == ins_len == 1:
                cls._safe_append(cigar, cls.OP_MATCH, 1)
            else:
                if del_len:
                    cigar.append((cls.OP_DEL, del_len))
                if ins_len:
                    cigar.append((cls.OP_INS, ins_len))
        
        return cigar