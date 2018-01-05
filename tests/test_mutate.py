import filecmp
import unittest

import pysam

import varlock_src.common as cmn
from tests.random import RandomMockup
from varlock_src.bam_mutator import BamMutator
from varlock_src.bdiff import BdiffIO
from varlock_src.random import VeryRandom
from varlock_src.vac import Vac


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
        cls._mut = BamMutator(cls.RESOURCE_PATH + 'input.bam')
        cls._rnd = VeryRandom(RandomMockup(0.5))
    
    def test_mutate_01(self):
        # EOF BAM case
        Vac.text2vac(self.RESOURCE_PATH + 'input_01.vac.txt', self.RESOURCE_PATH + 'input_01.vac')
        bdiff_file = self._mut.mutate(
            vac_filename=self.RESOURCE_PATH + 'input_01.vac',
            mut_bam_filename=self.RESOURCE_PATH + 'output_01.bam',
            secret=self.SECRET,
            rnd=self._rnd
        )
        
        self.assertEqual(21, self._mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(16, self._mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(8, self._mut.stat(BamMutator.STAT_VAC_COUNT))
        self.assertEqual(16, self._mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(6, self._mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        cmn.bam2sam(self.RESOURCE_PATH + 'output_01.bam', self.RESOURCE_PATH + 'output_01.sam')
        self.assertTrue(filecmp.cmp(
            self.RESOURCE_PATH + 'desired_01.sam',
            self.RESOURCE_PATH + 'output_01.sam')
        )
        
        BdiffIO.to_text_file(bdiff_file, self.RESOURCE_PATH + 'output_01.diff.txt')
        self.assertTrue(filecmp.cmp(
            self.RESOURCE_PATH + 'desired_01.diff.txt',
            self.RESOURCE_PATH + 'output_01.diff.txt')
        )
    
    def test_mutate_02(self):
        # EOF VAC case
        Vac.text2vac(self.RESOURCE_PATH + 'input_02.vac.txt', self.RESOURCE_PATH + 'input_02.vac')
        bdiff_file = self._mut.mutate(
            vac_filename=self.RESOURCE_PATH + 'input_02.vac',
            mut_bam_filename=self.RESOURCE_PATH + 'output_02.bam',
            secret=self.SECRET,
            rnd=self._rnd
        )
        
        self.assertEqual(21, self._mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(13, self._mut.stat(BamMutator.STAT_COVERING_COUNT))
        self.assertEqual(7, self._mut.stat(BamMutator.STAT_VAC_COUNT))
        self.assertEqual(14, self._mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(5, self._mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        cmn.bam2sam(self.RESOURCE_PATH + 'output_02.bam', self.RESOURCE_PATH + 'output_02.sam')
        self.assertTrue(filecmp.cmp(
            self.RESOURCE_PATH + 'desired_02.sam',
            self.RESOURCE_PATH + 'output_02.sam'
        ))
        
        BdiffIO.to_text_file(bdiff_file, self.RESOURCE_PATH + 'output_02.diff.txt')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired_02.diff.txt', self.RESOURCE_PATH + 'output_02.diff.txt')
        self.assertTrue(is_equal)


if __name__ == '__main__':
    # TODO test case with more references (chromosomes)
    unittest.main()
