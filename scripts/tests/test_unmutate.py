import filecmp
import unittest

import pysam

from varlock.bam_mutator import BamMutator
from varlock.common import bam2sam
from varlock.diff import Diff
from varlock.vac import Vac
from .random_mockup import RandomMockup


class TestUnmutate(unittest.TestCase):
    RESOURCE_PATH = 'tests/resources/unmutate/'
    
    @classmethod
    def setUpClass(cls):
        cls.mut = BamMutator(rnd=RandomMockup())
    
    def test_unmutate_01(self):
        # case 1
        # sam2bam(self.UNMUT_DIR_PATH + 'input_01.sam', self.UNMUT_DIR_PATH + 'input_01.bam')
        pysam.index(self.RESOURCE_PATH + 'input_01.bam')
        Diff.text2diff(self.RESOURCE_PATH + 'input_01.diff.txt', self.RESOURCE_PATH + 'input_01.diff')
        with open(self.RESOURCE_PATH + 'input_01.diff', 'rb') as diff_file:
            self.mut.unmutate(
                bam_filename=self.RESOURCE_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_01.bam'
            )
        
        self.assertEqual(15, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        # not 12 as in "mutate" because one snv is same as the reference (no diff record)
        self.assertEqual(11, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(5, self.mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(12, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(4, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_01.bam', self.RESOURCE_PATH + 'output_01.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_01_02.sam', self.RESOURCE_PATH + 'output_01.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_02(self):
        # case 2
        # sam2bam(self.UNMUT_DIR_PATH + 'input_02.sam', self.UNMUT_DIR_PATH + 'input_02.bam')
        pysam.index(self.RESOURCE_PATH + 'input_02.bam')
        Diff.text2diff(self.RESOURCE_PATH + 'input_02.diff.txt', self.RESOURCE_PATH + 'input_02.diff')
        with open(self.RESOURCE_PATH + 'input_02.diff', 'rb') as diff_file:
            self.mut.unmutate(
                bam_filename=self.RESOURCE_PATH + 'input_02.bam',
                diff_file=diff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_02.bam'
            )
        
        self.assertEqual(15, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(6, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(3, self.mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(7, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(3, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_02.bam', self.RESOURCE_PATH + 'output_02.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_01_02.sam', self.RESOURCE_PATH + 'output_02.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_03(self):
        # range supplied
        pysam.index(self.RESOURCE_PATH + 'input_01.bam')
        Diff.text2diff(self.RESOURCE_PATH + 'input_01.diff.txt', self.RESOURCE_PATH + 'input_01.diff')
        with open(self.RESOURCE_PATH + 'input_01.diff', 'rb') as diff_file:
            self.mut.unmutate(
                bam_filename=self.RESOURCE_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_03.bam',
                start_ref_name='chr22',
                start_ref_pos=1000021,
                end_ref_name='chr22',
                end_ref_pos=1000061,
            )
        
        self.assertEqual(7, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(5, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(3, self.mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(5, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(2, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_03.bam', self.RESOURCE_PATH + 'output_03.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_03.sam', self.RESOURCE_PATH + 'output_03.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_04(self):
        # range covers unmapped read
        pysam.index(self.RESOURCE_PATH + 'input_01.bam')
        Diff.text2diff(self.RESOURCE_PATH + 'input_01.diff.txt', self.RESOURCE_PATH + 'input_01.diff')
        with open(self.RESOURCE_PATH + 'input_01.diff', 'rb') as diff_file:
            self.mut.unmutate(
                bam_filename=self.RESOURCE_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_04.bam',
                start_ref_name='chr22',
                start_ref_pos=1000054,
                end_ref_name='chr22',
                end_ref_pos=1000078,
            )
        
        self.assertEqual(4, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(3, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(2, self.mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(4, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(2, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_04.bam', self.RESOURCE_PATH + 'output_04.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_04.sam', self.RESOURCE_PATH + 'output_04.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_05(self):
        # include unmapped
        pysam.index(self.RESOURCE_PATH + 'input_01.bam')
        Diff.text2diff(self.RESOURCE_PATH + 'input_01.diff.txt', self.RESOURCE_PATH + 'input_01.diff')
        with open(self.RESOURCE_PATH + 'input_01.diff', 'rb') as diff_file:
            self.mut.unmutate(
                bam_filename=self.RESOURCE_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_05.bam',
                include_unmapped=True
            )
        
        self.assertEqual(19, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(11, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(5, self.mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(12, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(4, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_05.bam', self.RESOURCE_PATH + 'output_05.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_05.sam', self.RESOURCE_PATH + 'output_05.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_06(self):
        # include unmapped with range
        pysam.index(self.RESOURCE_PATH + 'input_01.bam')
        Diff.text2diff(self.RESOURCE_PATH + 'input_01.diff.txt', self.RESOURCE_PATH + 'input_01.diff')
        with open(self.RESOURCE_PATH + 'input_01.diff', 'rb') as diff_file:
            self.mut.unmutate(
                bam_filename=self.RESOURCE_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_06.bam',
                start_ref_name='chr22',
                start_ref_pos=1000021,
                end_ref_name='chr22',
                end_ref_pos=1000061,
                include_unmapped=True
            )
        
        self.assertEqual(11, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(5, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(3, self.mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(5, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(2, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_06.bam', self.RESOURCE_PATH + 'output_06.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_06.sam', self.RESOURCE_PATH + 'output_06.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_07(self):
        # include unmapped
        pysam.index(self.RESOURCE_PATH + 'input_01.bam')
        Diff.text2diff(self.RESOURCE_PATH + 'input_01.diff.txt', self.RESOURCE_PATH + 'input_01.diff')
        with open(self.RESOURCE_PATH + 'input_01.diff', 'rb') as diff_file:
            self.mut.unmutate(
                bam_filename=self.RESOURCE_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_07.bam',
                unmapped_only=True
            )
        
        self.assertEqual(4, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(0, self.mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(0, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        # only first diff is read
        self.assertEqual(1, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_07.bam', self.RESOURCE_PATH + 'output_07.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_07.sam', self.RESOURCE_PATH + 'output_07.sam')
        self.assertEqual(True, is_equal)


if __name__ == '__main__':
    # TODO more unmutate slice tests
    # TODO test cases with alignments of different length
    # TODO test cases with more references (chromosomes)
    unittest.main()
