import filecmp
import os
import unittest

import pysam

import src.common as cmn
from src.bam_mutator import BamMutator
from src.bdiff import BdiffIO
from src.vac import Vac
from tests.random_mockup import VeryRandomMockup


class TestMaskSnv(unittest.TestCase):
    """
    Encrypted alignments are not tested here.
    In desired.sam they are replaced with those in actual output.sam to pass the test.
    """
    RESOURCE_PATH = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'resources/mask_snv/')

    SECRET = bytes([255]) * 16

    # TODO test cases with low coverage where personal alleles can not be determined

    @classmethod
    def setUpClass(cls):
        cmn.sam2bam(cls.RESOURCE_PATH + 'input.sam', cls.RESOURCE_PATH + 'input.bam')
        pysam.index(cls.RESOURCE_PATH + 'input.bam')
        cls._mut = BamMutator(cls.RESOURCE_PATH + 'input.bam')
        cls._rng = VeryRandomMockup()

    def test_mask(self):
        # EOF BAM case
        Vac.text2vac(self.RESOURCE_PATH + 'input.vac.txt', self.RESOURCE_PATH + 'input.vac')
        bdiff_file = self._mut.mutate(
            vac_filename=self.RESOURCE_PATH + 'input.vac',
            mut_bam_filename=self.RESOURCE_PATH + 'output.bam',
            secret=self.SECRET,
            mut_p=0,
            rng=self._rng
        )
        cmn.bam2sam(self.RESOURCE_PATH + 'output.bam', self.RESOURCE_PATH + 'output.sam')
        BdiffIO.to_text_file(bdiff_file, self.RESOURCE_PATH + 'output.diff.txt')

        self.assertTrue(filecmp.cmp(
            self.RESOURCE_PATH + 'desired.sam',
            self.RESOURCE_PATH + 'output.sam'
        ))

        self.assertTrue(filecmp.cmp(
            self.RESOURCE_PATH + 'desired.diff.txt',
            self.RESOURCE_PATH + 'output.diff.txt'
        ))


if __name__ == '__main__':
    # TODO test case with more references (chromosomes)
    unittest.main()
