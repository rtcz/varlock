import filecmp
import unittest

import pysam

from varlock.bdiff import BdiffIO
from varlock.bam_mutator import BamMutator
from varlock.common import bam2sam
from .random_mockup import RandomMockup


class TestUnmutate(unittest.TestCase):
    RESOURCE_PATH = 'tests/resources/unmutate/'
    
    @classmethod
    def setUpClass(cls):
        # TODO problem (bam->sam->bam) != bam
        # sam2bam(cls.RESOURCE_PATH + 'input.sam', cls.RESOURCE_PATH + 'input.bam')
        cls.mut = BamMutator(cls.RESOURCE_PATH + 'input.bam', RandomMockup())
        pysam.index(cls.RESOURCE_PATH + 'input.bam')
    
    def test_unmutate_01(self):
        # all mapped
        with BdiffIO.from_text_file(self.RESOURCE_PATH + 'input.diff.txt') as bdiff_file:
            self.mut.unmutate(
                bdiff_file=bdiff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_01.bam'
            )
        
        self.assertEqual(17, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        # one synonymous mutation is present
        self.assertEqual(15, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(16, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(6, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_01.bam', self.RESOURCE_PATH + 'output_01.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_01.sam', self.RESOURCE_PATH + 'output_01.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_03(self):
        # mapped within range
        with BdiffIO.from_text_file(self.RESOURCE_PATH + 'input.diff.txt') as bdiff_file:
            self.mut.unmutate(
                bdiff_file=bdiff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_03.bam',
                start_ref_name='chr2',
                start_ref_pos=1015,
                end_ref_name='chr2',
                end_ref_pos=1027,
            )
        
        self.assertEqual(9, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(9, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(8, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(3, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_03.bam', self.RESOURCE_PATH + 'output_03.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_03.sam', self.RESOURCE_PATH + 'output_03.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_05(self):
        # include unmapped
        with BdiffIO.from_text_file(self.RESOURCE_PATH + 'input.diff.txt') as bdiff_file:
            self.mut.unmutate(
                bdiff_file=bdiff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_05.bam',
                include_unmapped=True
            )
        
        self.assertEqual(21, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        # one synonymous mutation is present
        self.assertEqual(15, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(16, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(6, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_05.bam', self.RESOURCE_PATH + 'output_05.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_05.sam', self.RESOURCE_PATH + 'output_05.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_06(self):
        # include unmapped with range
        with BdiffIO.from_text_file(self.RESOURCE_PATH + 'input.diff.txt') as bdiff_file:
            self.mut.unmutate(
                bdiff_file=bdiff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_06.bam',
                start_ref_name='chr2',
                start_ref_pos=1015,
                end_ref_name='chr2',
                end_ref_pos=1027,
                include_unmapped=True
            )
        
        self.assertEqual(13, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(9, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(8, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(3, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_06.bam', self.RESOURCE_PATH + 'output_06.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_06.sam', self.RESOURCE_PATH + 'output_06.sam')
        self.assertEqual(True, is_equal)
    
    def test_unmutate_07(self):
        # include unmapped
        with BdiffIO.from_text_file(self.RESOURCE_PATH + 'input.diff.txt') as bdiff_file:
            self.mut.unmutate(
                bdiff_file=bdiff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output_07.bam',
                unmapped_only=True
            )
        
        self.assertEqual(4, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(0, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        # only the first diff record is read
        self.assertEqual(1, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(self.RESOURCE_PATH + 'output_07.bam', self.RESOURCE_PATH + 'output_07.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_07.sam', self.RESOURCE_PATH + 'output_07.sam')
        self.assertEqual(True, is_equal)


if __name__ == '__main__':
    # TODO more unmutate slice tests
    # TODO test cases with alignments of different length
    # TODO test cases with more references (chromosomes)
    unittest.main()
