import unittest

import pysam

from varlock_src.cigar import Cigar
from varlock_src.alignment import AlleleAlignment


class TestVariant(unittest.TestCase):
    @staticmethod
    def build_alignment():
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = "ATGC" * 10
        alignment.reference_start = 1000
        alignment.cigartuples = (
            (Cigar.OP_MATCH_ID, 10),
            (Cigar.OP_DEL_ID, 1),
            (Cigar.OP_EQUAL_ID, 3),
            (Cigar.OP_DIFF_ID, 3),
            (Cigar.OP_INS_ID, 1),
            (Cigar.OP_MATCH_ID, 20)
        )
        return alignment
    
    @staticmethod
    def build_flanked_alignment():
        """
        SAM/BAM files may include extra flanking bases that are not part of the alignment.
        These bases may be the result of the Smith-Waterman or other algorithms,
        which may not require alignments that begin at the first residue or end at the last.
        In addition, extra sequencing adapters, multiplex identifiers,
        and low-quality bases that were not considered for alignment may have been retained.
        """
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = "X" * 10 + "ATGC" * 10 + "X" * 10
        alignment.reference_start = 1000
        alignment.cigartuples = (
            (Cigar.OP_SOFT_CLIP_ID, 10),
            (Cigar.OP_MATCH_ID, 10),
            (Cigar.OP_DEL_ID, 1),
            (Cigar.OP_MATCH_ID, 9),
            (Cigar.OP_INS_ID, 1),
            (Cigar.OP_MATCH_ID, 20),
            (Cigar.OP_SOFT_CLIP_ID, 10)
        )
        return alignment
    
    def test_pos_conversion(self):
        alignment = self.build_alignment()
        self.assertIsNone(AlleleAlignment._ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(9, AlleleAlignment._ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
        
        alignment = self.build_flanked_alignment()
        self.assertIsNone(AlleleAlignment._ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(19, AlleleAlignment._ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
        
        # def test_replace_allele(self):
        #     self.assertTupleEqual(('TTAAAC', 'MMMMMM'), Cigar.replace_allele('TTAAAC', 'MMMMMM', ['A', 'AAAC'], 0, 0, 2))
        #     self.assertTupleEqual(('TTA', 'MMM'), Cigar.replace_allele('TTAAAC', 'MMMIII', ['A', 'AAAC'], 0, 0, 2))
        #
        #     self.assertTupleEqual(('TTA', 'MMM'), Cigar.replace_allele('TTAAAC', 'MMMIII', ['A', 'AAAC'], 0, 0, 2))
