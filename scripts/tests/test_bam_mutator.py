import filecmp
import unittest

import pysam

from varlock.bam_mutator import BamMutator
from varlock.common import sam2bam, bam2sam
from varlock.diff import Diff
from varlock.vac import Vac
from .random_mockup import RandomMockup


class TestBamMutator(unittest.TestCase):
    MUT_DIR_PATH = 'tests/resources/mutator/mutate/'
    UNMUT_DIR_PATH = 'tests/resources/mutator/unmutate/'
    
    @classmethod
    def setUpClass(cls):
        sam2bam(cls.MUT_DIR_PATH + 'input.sam', cls.MUT_DIR_PATH + 'input.bam')
        pysam.index(cls.MUT_DIR_PATH + 'input.bam')
    
    def test_mutate_01(self):
        # EOF BAM case
        Vac.text2vac(self.MUT_DIR_PATH + 'input_01.vac.txt', self.MUT_DIR_PATH + 'input_01.vac')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.MUT_DIR_PATH + 'output_01.diff', 'wb') as out_diff_file:
            mut.mutate(
                bam_filename=self.MUT_DIR_PATH + 'input.bam',
                vac_filename=self.MUT_DIR_PATH + 'input_01.vac',
                mut_bam_filename=self.MUT_DIR_PATH + 'output_01.bam',
                diff_file=out_diff_file
            )
            self.assertEqual(19, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(12, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(5, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(5, mut.stat(BamMutator.STAT_SNV_COUNT))
            self.assertEqual(12, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(4, mut.stat(BamMutator.STAT_DIFF_COUNT))
            
            bam2sam(self.MUT_DIR_PATH + 'output_01.bam', self.MUT_DIR_PATH + 'output_01.sam')
            is_equal = filecmp.cmp(self.MUT_DIR_PATH + 'desired_01.sam', self.MUT_DIR_PATH + 'output_01.sam')
            self.assertEqual(True, is_equal)
            Diff.diff2text(self.MUT_DIR_PATH + 'output_01.diff', self.MUT_DIR_PATH + 'output_01.diff.txt')
            is_equal = filecmp.cmp(self.MUT_DIR_PATH + 'desired_01.diff.txt', self.MUT_DIR_PATH + 'output_01.diff.txt')
            self.assertEqual(True, is_equal)
    
    def test_mutate_02(self):
        # EOF VAC case
        Vac.text2vac(self.MUT_DIR_PATH + 'input_02.vac.txt', self.MUT_DIR_PATH + 'input_02.vac')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.MUT_DIR_PATH + 'output_02.diff', 'wb') as out_diff_file:
            mut.mutate(
                bam_filename=self.MUT_DIR_PATH + 'input.bam',
                vac_filename=self.MUT_DIR_PATH + 'input_02.vac',
                mut_bam_filename=self.MUT_DIR_PATH + 'output_02.bam',
                diff_file=out_diff_file
            )
            
            self.assertEqual(19, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(7, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(3, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(4, mut.stat(BamMutator.STAT_SNV_COUNT))
            self.assertEqual(7, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(3, mut.stat(BamMutator.STAT_DIFF_COUNT))
            
            bam2sam(self.MUT_DIR_PATH + 'output_02.bam', self.MUT_DIR_PATH + 'output_02.sam')
            is_equal = filecmp.cmp(self.MUT_DIR_PATH + 'desired_02.sam', self.MUT_DIR_PATH + 'output_02.sam')
            self.assertEqual(True, is_equal)
            Diff.diff2text(self.MUT_DIR_PATH + 'output_02.diff', self.MUT_DIR_PATH + 'output_02.diff.txt')
            is_equal = filecmp.cmp(self.MUT_DIR_PATH + 'desired_02.diff.txt', self.MUT_DIR_PATH + 'output_02.diff.txt')
            self.assertEqual(True, is_equal)
    
    def test_unmutate_01(self):
        # case 1
        # sam2bam(self.UNMUT_DIR_PATH + 'input_01.sam', self.UNMUT_DIR_PATH + 'input_01.bam')
        pysam.index(self.UNMUT_DIR_PATH + 'input_01.bam')
        Diff.text2diff(self.UNMUT_DIR_PATH + 'input_01.diff.txt', self.UNMUT_DIR_PATH + 'input_01.diff')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.UNMUT_DIR_PATH + 'input_01.diff', 'rb') as diff_file:
            mut.unmutate(
                bam_filename=self.UNMUT_DIR_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.UNMUT_DIR_PATH + 'output_01.bam'
            )
            self.assertEqual(16, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            # not 12 as in "mutate" becouse one snv is same as the reference (no diff record)
            self.assertEqual(11, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(5, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(12, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(4, mut.stat(BamMutator.STAT_DIFF_COUNT))
            
            bam2sam(self.UNMUT_DIR_PATH + 'output_01.bam', self.UNMUT_DIR_PATH + 'output_01.sam')
            is_equal = filecmp.cmp(self.UNMUT_DIR_PATH + 'desired_01_02.sam', self.UNMUT_DIR_PATH + 'output_01.sam')
            self.assertEqual(True, is_equal)
    
    def test_unmutate_02(self):
        # case 2
        # sam2bam(self.UNMUT_DIR_PATH + 'input_02.sam', self.UNMUT_DIR_PATH + 'input_02.bam')
        pysam.index(self.UNMUT_DIR_PATH + 'input_02.bam')
        Diff.text2diff(self.UNMUT_DIR_PATH + 'input_02.diff.txt', self.UNMUT_DIR_PATH + 'input_02.diff')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.UNMUT_DIR_PATH + 'input_02.diff', 'rb') as diff_file:
            mut.unmutate(
                bam_filename=self.UNMUT_DIR_PATH + 'input_02.bam',
                diff_file=diff_file,
                out_bam_filename=self.UNMUT_DIR_PATH + 'output_02.bam'
            )
            self.assertEqual(16, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(6, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(3, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(7, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(3, mut.stat(BamMutator.STAT_DIFF_COUNT))
            
            bam2sam(self.UNMUT_DIR_PATH + 'output_02.bam', self.UNMUT_DIR_PATH + 'output_02.sam')
            is_equal = filecmp.cmp(self.UNMUT_DIR_PATH + 'desired_01_02.sam', self.UNMUT_DIR_PATH + 'output_02.sam')
            self.assertEqual(True, is_equal)
    
    def test_unmutate_03(self):
        # range supplied
        pysam.index(self.UNMUT_DIR_PATH + 'input_01.bam')
        Diff.text2diff(self.UNMUT_DIR_PATH + 'input_01.diff.txt', self.UNMUT_DIR_PATH + 'input_01.diff')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.UNMUT_DIR_PATH + 'input_01.diff', 'rb') as diff_file:
            mut.unmutate(
                bam_filename=self.UNMUT_DIR_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.UNMUT_DIR_PATH + 'output_03.bam',
                start_ref_name='chr22',
                start_ref_pos=1000021,
                end_ref_name='chr22',
                end_ref_pos=1000061,
            )
            self.assertEqual(7, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(5, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(3, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(5, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(2, mut.stat(BamMutator.STAT_DIFF_COUNT))
            
            bam2sam(self.UNMUT_DIR_PATH + 'output_03.bam', self.UNMUT_DIR_PATH + 'output_03.sam')
            is_equal = filecmp.cmp(self.UNMUT_DIR_PATH + 'desired_03.sam', self.UNMUT_DIR_PATH + 'output_03.sam')
            self.assertEqual(True, is_equal)
    
    def test_unmutate_04(self):
        # range covers unmapped read
        pysam.index(self.UNMUT_DIR_PATH + 'input_01.bam')
        Diff.text2diff(self.UNMUT_DIR_PATH + 'input_01.diff.txt', self.UNMUT_DIR_PATH + 'input_01.diff')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.UNMUT_DIR_PATH + 'input_01.diff', 'rb') as diff_file:
            mut.unmutate(
                bam_filename=self.UNMUT_DIR_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.UNMUT_DIR_PATH + 'output_04.bam',
                start_ref_name='chr22',
                start_ref_pos=1000054,
                end_ref_name='chr22',
                end_ref_pos=1000078,
            )
            self.assertEqual(5, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(3, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(2, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(4, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(2, mut.stat(BamMutator.STAT_DIFF_COUNT))
            
            bam2sam(self.UNMUT_DIR_PATH + 'output_04.bam', self.UNMUT_DIR_PATH + 'output_04.sam')
            is_equal = filecmp.cmp(self.UNMUT_DIR_PATH + 'desired_04.sam', self.UNMUT_DIR_PATH + 'output_04.sam')
            self.assertEqual(True, is_equal)
    
    def test_unmutate_05(self):
        # include unmapped
        pysam.index(self.UNMUT_DIR_PATH + 'input_01.bam')
        Diff.text2diff(self.UNMUT_DIR_PATH + 'input_01.diff.txt', self.UNMUT_DIR_PATH + 'input_01.diff')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.UNMUT_DIR_PATH + 'input_01.diff', 'rb') as diff_file:
            mut.unmutate(
                bam_filename=self.UNMUT_DIR_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.UNMUT_DIR_PATH + 'output_05.bam',
                include_unmapped=True
            )
            self.assertEqual(19, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(11, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(5, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(12, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(4, mut.stat(BamMutator.STAT_DIFF_COUNT))
            
            bam2sam(self.UNMUT_DIR_PATH + 'output_05.bam', self.UNMUT_DIR_PATH + 'output_05.sam')
            is_equal = filecmp.cmp(self.UNMUT_DIR_PATH + 'desired_05.sam', self.UNMUT_DIR_PATH + 'output_05.sam')
            self.assertEqual(True, is_equal)

    def test_unmutate_06(self):
        # include unmapped with range
        pysam.index(self.UNMUT_DIR_PATH + 'input_01.bam')
        Diff.text2diff(self.UNMUT_DIR_PATH + 'input_01.diff.txt', self.UNMUT_DIR_PATH + 'input_01.diff')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.UNMUT_DIR_PATH + 'input_01.diff', 'rb') as diff_file:
            mut.unmutate(
                bam_filename=self.UNMUT_DIR_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.UNMUT_DIR_PATH + 'output_06.bam',
                start_ref_name='chr22',
                start_ref_pos=1000021,
                end_ref_name='chr22',
                end_ref_pos=1000061,
                include_unmapped=True
            )
            self.assertEqual(11, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(5, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(3, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(5, mut.stat(BamMutator.STAT_MUT_COUNT))
            self.assertEqual(2, mut.stat(BamMutator.STAT_DIFF_COUNT))
        
            bam2sam(self.UNMUT_DIR_PATH + 'output_06.bam', self.UNMUT_DIR_PATH + 'output_06.sam')
            is_equal = filecmp.cmp(self.UNMUT_DIR_PATH + 'desired_06.sam', self.UNMUT_DIR_PATH + 'output_06.sam')
            self.assertEqual(True, is_equal)

    def test_unmutate_07(self):
        # include unmapped
        pysam.index(self.UNMUT_DIR_PATH + 'input_01.bam')
        Diff.text2diff(self.UNMUT_DIR_PATH + 'input_01.diff.txt', self.UNMUT_DIR_PATH + 'input_01.diff')
        mut = BamMutator(rnd=RandomMockup())
        with open(self.UNMUT_DIR_PATH + 'input_01.diff', 'rb') as diff_file:
            mut.unmutate(
                bam_filename=self.UNMUT_DIR_PATH + 'input_01.bam',
                diff_file=diff_file,
                out_bam_filename=self.UNMUT_DIR_PATH + 'output_07.bam',
                unmapped_only=True
            )
            self.assertEqual(4, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
            self.assertEqual(0, mut.stat(BamMutator.STAT_MUT_ALIGNMENT_COUNT))
            self.assertEqual(0, mut.stat(BamMutator.STAT_MAX_COVERAGE))
            self.assertEqual(0, mut.stat(BamMutator.STAT_MUT_COUNT))
            # only first diff is read
            self.assertEqual(1, mut.stat(BamMutator.STAT_DIFF_COUNT))
        
            bam2sam(self.UNMUT_DIR_PATH + 'output_07.bam', self.UNMUT_DIR_PATH + 'output_07.sam')
            is_equal = filecmp.cmp(self.UNMUT_DIR_PATH + 'desired_07.sam', self.UNMUT_DIR_PATH + 'output_07.sam')
            self.assertEqual(True, is_equal)


if __name__ == '__main__':
    # TODO more unmutate slice tests
    # TODO test cases with alignments of different length
    # TODO test cases with more references (chromosomes)
    unittest.main()
