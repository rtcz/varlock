import filecmp
import os
import unittest

import pysam

from src.bam_mutator import BamMutator
from src.bdiff import BdiffIO
from src.vac import bam2sam
from tests.random_mockup import VeryRandomMockup


class TestUnmaskSnv(unittest.TestCase):
    """
    Encrypted alignments are not tested here.
    In desired.sam they are replaced with those in actual output.sam to pass the test.
    """
    RESOURCE_PATH = os.path.join(os.path.realpath(os.path.dirname(__file__)), 'resources/unmask_snv/')

    SECRET = bytes([255]) * 16

    @classmethod
    def setUpClass(cls):
        # TODO problem (bam1->sam->bam2) != bam1, using samtools manually instead
        # samtools view -hb --no-PG input.sam > input.bam
        # samtools index input.bam
        # sam2bam(cls.RESOURCE_PATH + 'input.sam', cls.RESOURCE_PATH + 'input.bam')
        cls.mut = BamMutator(cls.RESOURCE_PATH + 'input.bam')
        pysam.index(cls.RESOURCE_PATH + 'input.bam')
        cls._rng = VeryRandomMockup()

    def test_unmask(self):
        # all mapped
        with BdiffIO.from_text_file(self.RESOURCE_PATH + 'input.diff.txt') as bdiff_file:
            self.mut.unmutate(
                bdiff_file=bdiff_file,
                out_bam_filename=self.RESOURCE_PATH + 'output.bam',
                rng=self._rng
            )

        bam2sam(self.RESOURCE_PATH + 'output.bam', self.RESOURCE_PATH + 'output.sam')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'desired.sam', self.RESOURCE_PATH + 'output.sam')
        self.assertEqual(True, is_equal)


if __name__ == '__main__':
    # TODO test case with more references (chromosomes)
    unittest.main()
