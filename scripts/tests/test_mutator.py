import filecmp
import unittest
from random import Random

import pysam

from varlock import *
from .random_mockup import RandomMockup


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
        self.assertListEqual(['A'], Mutator.get_base_pileup([SnvAlignment(alignment, 0)]))
        self.assertListEqual(['C'], Mutator.get_base_pileup([SnvAlignment(alignment, 3)]))
        self.assertListEqual([], Mutator.get_base_pileup([SnvAlignment(alignment, None)]))
        self.assertListEqual(['C'], Mutator.get_base_pileup([SnvAlignment(alignment, 39)]))
        
        alignment = self.build_flanked_alignment()
        self.assertListEqual(['A'], Mutator.get_base_pileup([SnvAlignment(alignment, 10)]))
        self.assertListEqual(['C'], Mutator.get_base_pileup([SnvAlignment(alignment, 13)]))
    
    def test_base_setting(self):
        # regular alignment
        alignment = self.build_alignment()
        Mutator.set_base(alignment=alignment, pos=0, base='C')
        Mutator.set_base(alignment=alignment, pos=3, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        Mutator.set_base(alignment=alignment, pos=10, base='C')
        Mutator.set_base(alignment=alignment, pos=13, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
    
    def test_base_getting(self):
        # regular alignment
        alignment = self.build_alignment()
        self.assertEqual('A', Mutator.get_base(alignment=alignment, pos=0))
        self.assertEqual('C', Mutator.get_base(alignment=alignment, pos=3))
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        self.assertEqual('A', Mutator.get_base(alignment=alignment, pos=10))
        self.assertEqual('C', Mutator.get_base(alignment=alignment, pos=13))
    
    def test_pos_conversion(self):
        alignment = self.build_alignment()
        self.assertIsNone(Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(9, Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
        
        alignment = self.build_flanked_alignment()
        self.assertIsNone(Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(19, Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
    
    def test_cigar(self):
        alignment = self.build_alignment()
        Mutator.check_cigar_str(alignment=alignment, snv_pos=9)
        self.assertRaises(ValueError, lambda: Mutator.check_cigar_str(alignment=alignment, snv_pos=10))
        self.assertRaises(ValueError, lambda: Mutator.check_cigar_str(alignment=alignment, snv_pos=11))
    
    def test_mut_map(self):
        mut_map = Mutator.create_mut_map(alt_ac=[1, 1, 1, 1], ref_ac=[1, 1, 1, 1], rnd=Random(0))
        self.assertDictEqual({'A': 'G', 'T': 'T', 'G': 'C', 'C': 'A', 'N': 'N'}, mut_map)
    
    def test_mutate(self):
        dir_path = 'tests/resources/mutator/mutate/'
        sam2bam(dir_path + 'input.sam', dir_path + 'input.bam')
        pysam.index(dir_path + 'input.bam')
        
        # EOF BAM case
        Vac.text2vac(dir_path + 'input_01.vac.txt', dir_path + 'input_01.vac')
        mut = MutatorCaller(rnd=RandomMockup())
        diff_file = mut.mutate(
            bam_filename=dir_path + 'input.bam',
            vac_filename=dir_path + 'input_01.vac',
            out_bam_filename=dir_path + 'output_01.bam'
        )
        with open(dir_path + 'output_01.diff', 'wb') as out_diff_file, diff_file:
            out_diff_file.write(diff_file.read())
        
        self.assertEqual(15, mut.stat(MutatorCaller.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(MutatorCaller.STAT_UNMAPPED_COUNT))
        self.assertEqual(12, mut.stat(MutatorCaller.STAT_OVERLAPPING_COUNT))
        self.assertEqual(5, mut.stat(MutatorCaller.STAT_MAX_COVERAGE))
        self.assertEqual(5, mut.stat(MutatorCaller.STAT_SNV_COUNT))
        self.assertEqual(12, mut.stat(MutatorCaller.STAT_MUT_COUNT))
        self.assertEqual(4, mut.stat(MutatorCaller.STAT_DIFF_COUNT))
        
        bam2sam(dir_path + 'output_01.bam', dir_path + 'output_01.sam')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired_01.sam', dir_path + 'output_01.sam'))
        Diff.diff2text(dir_path + 'output_01.diff', dir_path + 'output_01.diff.txt')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired_01.diff.txt', dir_path + 'output_01.diff.txt'))
        
        # EOF VAC case
        Vac.text2vac(dir_path + 'input_02.vac.txt', dir_path + 'input_02.vac')
        mut = MutatorCaller(rnd=RandomMockup())
        diff_file = mut.mutate(
            bam_filename=dir_path + 'input.bam',
            vac_filename=dir_path + 'input_02.vac',
            out_bam_filename=dir_path + 'output_02.bam'
        )
        with open(dir_path + 'output_02.diff', 'wb') as out_diff_file, diff_file:
            out_diff_file.write(diff_file.read())
        
        self.assertEqual(15, mut.stat(MutatorCaller.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(MutatorCaller.STAT_UNMAPPED_COUNT))
        self.assertEqual(7, mut.stat(MutatorCaller.STAT_OVERLAPPING_COUNT))
        self.assertEqual(3, mut.stat(MutatorCaller.STAT_MAX_COVERAGE))
        self.assertEqual(4, mut.stat(MutatorCaller.STAT_SNV_COUNT))
        self.assertEqual(7, mut.stat(MutatorCaller.STAT_MUT_COUNT))
        self.assertEqual(3, mut.stat(MutatorCaller.STAT_DIFF_COUNT))
        
        bam2sam(dir_path + 'output_02.bam', dir_path + 'output_02.sam')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired_02.sam', dir_path + 'output_02.sam'))
        Diff.diff2text(dir_path + 'output_02.diff', dir_path + 'output_02.diff.txt')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired_02.diff.txt', dir_path + 'output_02.diff.txt'))
        # TODO test case with more references (chromosomes)
    
    def test_unmutate(self):
        dir_path = 'tests/resources/mutator/unmutate/'
        
        # case 1
        pysam.index(dir_path + 'input_01.bam')
        Diff.text2diff(dir_path + 'input_01.diff.txt', dir_path + 'input_01.diff')
        mut = MutatorCaller(rnd=RandomMockup())
        mut.unmutate(
            bam_filename=dir_path + 'input_01.bam',
            diff_file=dir_path + 'input_01.diff',
            out_bam_filename=dir_path + 'output_01.bam'
        )
        self.assertEqual(15, mut.stat(MutatorCaller.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(MutatorCaller.STAT_UNMAPPED_COUNT))
        self.assertEqual(11, mut.stat(MutatorCaller.STAT_OVERLAPPING_COUNT))
        self.assertEqual(5, mut.stat(MutatorCaller.STAT_MAX_COVERAGE))
        self.assertEqual(12, mut.stat(MutatorCaller.STAT_MUT_COUNT))
        self.assertEqual(4, mut.stat(MutatorCaller.STAT_DIFF_COUNT))
        
        bam2sam(dir_path + 'output_01.bam', dir_path + 'output_01.sam')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired.sam', dir_path + 'output_01.sam'))
        
        # case 2
        pysam.index(dir_path + 'input_02.bam')
        Diff.text2diff(dir_path + 'input_02.diff.txt', dir_path + 'input_02.diff')
        mut = MutatorCaller(rnd=RandomMockup())
        mut.unmutate(
            bam_filename=dir_path + 'input_02.bam',
            diff_file=dir_path + 'input_02.diff',
            out_bam_filename=dir_path + 'output_02.bam'
        )
        self.assertEqual(15, mut.stat(MutatorCaller.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(MutatorCaller.STAT_UNMAPPED_COUNT))
        self.assertEqual(6, mut.stat(MutatorCaller.STAT_OVERLAPPING_COUNT))
        self.assertEqual(3, mut.stat(MutatorCaller.STAT_MAX_COVERAGE))
        self.assertEqual(7, mut.stat(MutatorCaller.STAT_MUT_COUNT))
        self.assertEqual(3, mut.stat(MutatorCaller.STAT_DIFF_COUNT))

        bam2sam(dir_path + 'output_02.bam', dir_path + 'output_02.sam')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired.sam', dir_path + 'output_02.sam'))
        
        # TODO more tests (ranges etc.)
        # TODO test case with more references (chromosomes)


if __name__ == '__main__':
    # python3 /usr/local/bin/nosetests -s /data/projects/varlock/scripts/mutator_test.py
    # python3 -m cProfile -s tottime /data/projects/varlock/scripts/mutator.py
    unittest.main()
