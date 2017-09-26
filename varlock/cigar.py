from varlock import po


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
    
    @classmethod
    def _safe_append(cls, cigar_tuples: list, op_type: int, op_length: int):
        """
        Merging new OP with last OP from cigar_tuples
        when OPs are of the same type.
        Otherwise new OP is appended to cigar_tuples.
        :param cigar_tuples:
        :return:
        """
        if len(cigar_tuples):
            prev_type, prev_length = cigar_tuples[-1]
            if prev_type == op_type:
                # merge OPs
                cigar_tuples[-1] = (op_type, op_length + prev_length)
            else:
                cigar_tuples.append((op_type, op_length))
        else:
            cigar_tuples.append((op_type, op_length))
    
    @classmethod
    def del_subrange(cls, cigar_tuples: list, start_pos: int, end_pos: int):
        """
        :param cigar_tuples: [<op_type, op_length>, ...]
        :param start_pos: position in AlignedSegment.query_sequence
        :param end_pos: position in AlignedSegment.query_sequence
        :return:
        """
        assert end_pos > start_pos >= 0
        new_cigar = []
        curr_pos = 0
        started = None
        ended = None
        normalize = False
        for i in range(len(cigar_tuples)):
            curr_type, curr_length = cigar_tuples[i]
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
            elif normalize:
                normalize = False
                cls._safe_append(new_cigar, curr_type, curr_length)
            else:
                new_cigar.append((curr_type, curr_length))
        
        return new_cigar
    
    @classmethod
    def place_op(cls, cigar_tuples: list, op_pos: int, op_type: int, op_length: int):
        """
        Place CIGAR operation into existing CIGAR string (in form of list of tuples).
        Function does not support placing hard and soft clip operations
        since they can appear only at CIGAR ends. Moreover HARD clipped bases
        are not part of AlignedSegment.query_sequence.
        :param cigar_tuples: [<op_type, op_length>, ...]
        :param op_pos: position in AlignedSegment.query_sequence
        :param op_type: CIGAR operation type id
        :param op_length: CIGAR operation length
        :return: New CIGAR in form of tuples.
        """
        assert op_type in [cls.OP_MATCH, cls.OP_EQUAL, cls.OP_DIFF, cls.OP_INS, cls.OP_PAD, cls.OP_DEL]
        assert op_length > 0
        assert op_pos >= 0
        curr_pos = 0
        new_cigar = []
        is_placed = False
        
        for i in range(len(cigar_tuples)):
            curr_type, curr_length = cigar_tuples[i]
            # OP_HARD_CLIP is not contained in cigar_tuples nor query_sequence
            if not is_placed and curr_type != cls.OP_DEL:
                # OP_MATCH, OP_EQUAL, OP_DIFF, OP_INS, OP_PAD
                if op_pos >= curr_pos + curr_length:
                    new_cigar.append((curr_type, curr_length))
                elif op_type == curr_type:
                    new_cigar.append((op_type, curr_length + op_length))
                else:
                    new_cigar.append((curr_type, curr_pos - op_length))
                    new_cigar.append((op_type, op_length))
                    new_cigar.append((curr_type, curr_pos - curr_length + op_length))
                
                if curr_pos == op_pos:
                    # place the OP before a current OP
                    if curr_type == op_type:
                        # current OP is equal to the OP, merge them
                        new_cigar.append((op_type, curr_length + op_length))
                    elif len(new_cigar) and new_cigar[-1][0] == op_type:
                        # previous OP is equal to the OP, merge them
                        new_cigar[-1] = (new_cigar[-1][0], new_cigar[-1][1] + op_length)
                    else:
                        # place the OP in begining of CIGAR
                        new_cigar.append((op_type, op_length))
                        new_cigar.append((curr_type, curr_length))
                    
                    # place a current OP after the OP
                    is_placed = True
                elif op_pos < curr_pos < curr_pos + curr_length:
                    # place the OP inside a current OP
                    new_cigar.append((curr_type, curr_pos - op_length))
                    new_cigar.append((op_type, op_length))
                    new_cigar.append((curr_type, curr_pos - curr_length + op_length))
                    
                    is_placed = True
                
                curr_pos += curr_length
            else:
                new_cigar.append((curr_type, curr_length))
        
        if not is_placed and curr_pos == op_pos:
            # place OP
            if new_cigar[-1][0] == op_type:
                new_cigar[-1] = (new_cigar[-1][0], new_cigar[-1][1] + op_length)
            else:
                new_cigar.append((op_type, op_length))
        
        return new_cigar
