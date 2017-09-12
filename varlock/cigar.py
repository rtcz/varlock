from varlock import po


class Cigar:
    CIGAR_MATCH = 0  # M
    CIGAR_INS = 1  # I
    CIGAR_DEL = 2  # D
    CIGAR_REF_SKIP = 3  # N
    CIGAR_SOFT_CLIP = 4  # S
    CIGAR_HARD_CLIP = 5  # H
    CIGAR_PAD = 6  # P
    CIGAR_EQUAL = 7  # E
    CIGAR_DIFF = 8  # X
    NM_TAG = 9
    
    @classmethod
    def validate(cls, variant: po.AlignedVariant):
        """
        Validate BAM's cigar operation at mutated variant position.
        
        :param variant: mutated variant
        """
        # TODO: update cigar and quality strings.
        pass
        # start at -1 to work with 0-based position
        # seq_pos = -1
        # for cig_op_id, cig_op_len in variant.alignment.cigartuples:
        #     if cig_op_id == cls.CIGAR_MATCH or cig_op_id == cls.CIGAR_INS:
        #         # match and insert are present in aligned sequence
        #         tmp_seq_pos = seq_pos + cig_op_len
        #     elif cig_op_id == cls.CIGAR_DEL:
        #         # deletion is not present in aligned sequence
        #         continue
        #     else:
        #         raise ValueError("Unsupported CIGAR operation %d" % cig_op_id)
        #
        #     if tmp_seq_pos < snv_pos:
        #         # curr_seq_pos has not reached snv_seq_pos
        #         seq_pos = tmp_seq_pos
        #     else:
        #         # cigar segment covers snv position, it must be a match
        #         if cig_op_id != cls.CIGAR_MATCH:
        #             # snv is not matched
        #             # TODO update mapping quality
        #             # TODO update cigar string
        #
        #             print("snv_pos %d" % snv_pos)
        #             print("seq_pos %d" % seq_pos)
        #             print(alignment.cigarstring)
        #             print(alignment.get_aligned_pairs(matches_only=False, with_seq=False))
        #
        #             raise ValueError("Unsupported CIGAR operation %d" % cig_op_id)
        #         break
