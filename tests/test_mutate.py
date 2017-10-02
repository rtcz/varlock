import filecmp
import unittest

import pysam

from varlock.bdiff import BdiffIO
import varlock.common as cmn
from varlock.bam_mutator import BamMutator
from varlock.vac import Vac
from tests.random_mockup import RandomMockup


class TestMutate(unittest.TestCase):
    """
    Encrypted alignments are not tested here.
    In desired.sam they are replaced with those in actual output.sam to pass the test.
    """
    RESOURCE_PATH = 'tests/resources/mutate/'
    
    SECRET = bytes([255]) * 16
    
    @classmethod
    def setUpClass(cls):
        cmn.sam2bam(cls.RESOURCE_PATH + 'input.sam', cls.RESOURCE_PATH + 'input.bam')
        pysam.index(cls.RESOURCE_PATH + 'input.bam')
        cls.mut = BamMutator(cls.RESOURCE_PATH + 'input.bam', RandomMockup())
    
    def test_mutate_01(self):
        # EOF BAM case
        Vac.text2vac(self.RESOURCE_PATH + 'input_01.vac.txt', self.RESOURCE_PATH + 'input_01.vac')
        bdiff_file = self.mut.mutate(
            vac_filename=self.RESOURCE_PATH + 'input_01.vac',
            mut_bam_filename=self.RESOURCE_PATH + 'output_01.bam',
            secret=self.SECRET
        )
        
        self.assertEqual(21, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(16, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(8, self.mut.stat(BamMutator.STAT_VAC_COUNT))
        self.assertEqual(16, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(6, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        cmn.bam2sam(self.RESOURCE_PATH + 'output_01.bam', self.RESOURCE_PATH + 'output_01.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_01.sam', self.RESOURCE_PATH + 'output_01.sam')
        self.assertTrue(is_equal)
        
        BdiffIO.to_text_file(bdiff_file, self.RESOURCE_PATH + 'output_01.diff.txt')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_01.diff.txt', self.RESOURCE_PATH + 'output_01.diff.txt')
        self.assertTrue(is_equal)
    
    def test_mutate_02(self):
        # EOF VAC case
        Vac.text2vac(self.RESOURCE_PATH + 'input_02.vac.txt', self.RESOURCE_PATH + 'input_02.vac')
        bdiff_file = self.mut.mutate(
            vac_filename=self.RESOURCE_PATH + 'input_02.vac',
            mut_bam_filename=self.RESOURCE_PATH + 'output_02.bam',
            secret=self.SECRET
        )
        
        self.assertEqual(21, self.mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(13, self.mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(7, self.mut.stat(BamMutator.STAT_VAC_COUNT))
        self.assertEqual(14, self.mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(5, self.mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        cmn.bam2sam(self.RESOURCE_PATH + 'output_02.bam', self.RESOURCE_PATH + 'output_02.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_02.sam', self.RESOURCE_PATH + 'output_02.sam')
        self.assertTrue(is_equal)
        
        BdiffIO.to_text_file(bdiff_file, self.RESOURCE_PATH + 'output_02.diff.txt')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_02.diff.txt', self.RESOURCE_PATH + 'output_02.diff.txt')
        self.assertTrue(is_equal)


if __name__ == '__main__':
    # TODO test cases with alignments of different length
    # TODO test cases with more references (chromosomes)
    unittest.main()
