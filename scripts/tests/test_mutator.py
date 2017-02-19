import unittest
from random import Random

import pysam

from varlock.common import *
from varlock.mutator import Mutator
from varlock.po import SnvAlignment


class TestMutator(unittest.TestCase):
    @staticmethod
    def build_alignment():
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = "ATGC" * 10
        alignment.reference_start = 1000
        alignment.cigartuples = (
            (Mutator.CIGAR_MATCH, 10),
            (Mutator.CIGAR_DEL, 1),
            (Mutator.CIGAR_EQUAL, 3),
            (Mutator.CIGAR_DIFF, 3),
            (Mutator.CIGAR_INS, 1),
            (Mutator.CIGAR_MATCH, 20)
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
            (Mutator.CIGAR_SOFT_CLIP, 10),
            (Mutator.CIGAR_MATCH, 10),
            (Mutator.CIGAR_DEL, 1),
            (Mutator.CIGAR_MATCH, 9),
            (Mutator.CIGAR_INS, 1),
            (Mutator.CIGAR_MATCH, 20),
            (Mutator.CIGAR_SOFT_CLIP, 10)
        )
        return alignment
    
    def test_pileup(self):
        alignment = self.build_alignment()
        self.assertListEqual(['A'], get_base_pileup([SnvAlignment(alignment, 0)]))
        self.assertListEqual(['C'], get_base_pileup([SnvAlignment(alignment, 3)]))
        self.assertListEqual([], get_base_pileup([SnvAlignment(alignment, None)]))
        self.assertListEqual(['C'], get_base_pileup([SnvAlignment(alignment, 39)]))
        
        alignment = self.build_flanked_alignment()
        self.assertListEqual(['A'], get_base_pileup([SnvAlignment(alignment, 10)]))
        self.assertListEqual(['C'], get_base_pileup([SnvAlignment(alignment, 13)]))
    
    def test_base_setting(self):
        # regular alignment
        alignment = self.build_alignment()
        set_base(alignment=alignment, pos=0, base='C')
        set_base(alignment=alignment, pos=3, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        set_base(alignment=alignment, pos=10, base='C')
        set_base(alignment=alignment, pos=13, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
    
    def test_base_getting(self):
        # regular alignment
        alignment = self.build_alignment()
        self.assertEqual('A', get_base(alignment=alignment, pos=0))
        self.assertEqual('C', get_base(alignment=alignment, pos=3))
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        self.assertEqual('A', get_base(alignment=alignment, pos=10))
        self.assertEqual('C', get_base(alignment=alignment, pos=13))
    
    def test_pos_conversion(self):
        alignment = self.build_alignment()
        self.assertIsNone(ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(9, ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
        
        alignment = self.build_flanked_alignment()
        self.assertIsNone(ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(19, ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
    
    def test_cigar(self):
        alignment = self.build_alignment()
        Mutator.check_cigar_str(alignment=alignment, snv_pos=9)
        self.assertRaises(ValueError, lambda: Mutator.check_cigar_str(alignment=alignment, snv_pos=10))
        self.assertRaises(ValueError, lambda: Mutator.check_cigar_str(alignment=alignment, snv_pos=11))
    
    def test_mut_map(self):
        mut_map = Mutator.create_mut_map(alt_ac=[1, 1, 1, 1], ref_ac=[1, 1, 1, 1], rnd=Random(0))
        self.assertDictEqual({'A': 'G', 'T': 'T', 'G': 'C', 'C': 'A', 'N': 'N'}, mut_map)


if __name__ == '__main__':
    # python3 /usr/local/bin/nosetests -s /data/projects/varlock/scripts/mutator_test.py
    # python3 -m cProfile -s tottime /data/projects/varlock/scripts/mutator.py
    unittest.main()
