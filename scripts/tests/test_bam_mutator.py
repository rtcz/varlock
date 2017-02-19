import filecmp
import unittest

import pysam

from varlock.bam_mutator import BamMutator
from varlock.common import sam2bam, bam2sam
from varlock.diff import Diff
from varlock.vac import Vac
from .random_mockup import RandomMockup


class TestBamMutator(unittest.TestCase):
    def test_mutate(self):
        dir_path = 'tests/resources/mutator/mutate/'
        sam2bam(dir_path + 'input.sam', dir_path + 'input.bam')
        pysam.index(dir_path + 'input.bam')
        
        # EOF BAM case
        Vac.text2vac(dir_path + 'input_01.vac.txt', dir_path + 'input_01.vac')
        mut = BamMutator(rnd=RandomMockup())
        diff_file = mut.mutate(
            bam_filename=dir_path + 'input.bam',
            vac_filename=dir_path + 'input_01.vac',
            out_bam_filename=dir_path + 'output_01.bam'
        )
        with open(dir_path + 'output_01.diff', 'wb') as out_diff_file, diff_file:
            out_diff_file.write(diff_file.read())
        
        self.assertEqual(15, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(BamMutator.STAT_UNMAPPED_COUNT))
        self.assertEqual(12, mut.stat(BamMutator.STAT_OVERLAPPING_COUNT))
        self.assertEqual(5, mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(5, mut.stat(BamMutator.STAT_SNV_COUNT))
        self.assertEqual(12, mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(4, mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(dir_path + 'output_01.bam', dir_path + 'output_01.sam')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired_01.sam', dir_path + 'output_01.sam'))
        Diff.diff2text(dir_path + 'output_01.diff', dir_path + 'output_01.diff.txt')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired_01.diff.txt', dir_path + 'output_01.diff.txt'))
        
        # EOF VAC case
        Vac.text2vac(dir_path + 'input_02.vac.txt', dir_path + 'input_02.vac')
        mut = BamMutator(rnd=RandomMockup())
        diff_file = mut.mutate(
            bam_filename=dir_path + 'input.bam',
            vac_filename=dir_path + 'input_02.vac',
            out_bam_filename=dir_path + 'output_02.bam'
        )
        with open(dir_path + 'output_02.diff', 'wb') as out_diff_file, diff_file:
            out_diff_file.write(diff_file.read())
        
        self.assertEqual(15, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(BamMutator.STAT_UNMAPPED_COUNT))
        self.assertEqual(7, mut.stat(BamMutator.STAT_OVERLAPPING_COUNT))
        self.assertEqual(3, mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(4, mut.stat(BamMutator.STAT_SNV_COUNT))
        self.assertEqual(7, mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(3, mut.stat(BamMutator.STAT_DIFF_COUNT))
        
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
        mut = BamMutator(rnd=RandomMockup())
        mut.unmutate(
            bam_filename=dir_path + 'input_01.bam',
            diff_file=dir_path + 'input_01.diff',
            out_bam_filename=dir_path + 'output_01.bam'
        )
        self.assertEqual(15, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(BamMutator.STAT_UNMAPPED_COUNT))
        self.assertEqual(11, mut.stat(BamMutator.STAT_OVERLAPPING_COUNT))
        self.assertEqual(5, mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(12, mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(4, mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(dir_path + 'output_01.bam', dir_path + 'output_01.sam')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired.sam', dir_path + 'output_01.sam'))
        
        # case 2
        pysam.index(dir_path + 'input_02.bam')
        Diff.text2diff(dir_path + 'input_02.diff.txt', dir_path + 'input_02.diff')
        mut = BamMutator(rnd=RandomMockup())
        mut.unmutate(
            bam_filename=dir_path + 'input_02.bam',
            diff_file=dir_path + 'input_02.diff',
            out_bam_filename=dir_path + 'output_02.bam'
        )
        self.assertEqual(15, mut.stat(BamMutator.STAT_ALIGNMENT_COUNT))
        self.assertEqual(0, mut.stat(BamMutator.STAT_UNMAPPED_COUNT))
        self.assertEqual(6, mut.stat(BamMutator.STAT_OVERLAPPING_COUNT))
        self.assertEqual(3, mut.stat(BamMutator.STAT_MAX_COVERAGE))
        self.assertEqual(7, mut.stat(BamMutator.STAT_MUT_COUNT))
        self.assertEqual(3, mut.stat(BamMutator.STAT_DIFF_COUNT))
        
        bam2sam(dir_path + 'output_02.bam', dir_path + 'output_02.sam')
        self.assertEqual(True, filecmp.cmp(dir_path + 'desired.sam', dir_path + 'output_02.sam'))
        
        # case 3
        pysam.index(dir_path + 'input_01.bam')
        Diff.text2diff(dir_path + 'input_01.diff.txt', dir_path + 'input_01.diff')
        mut = BamMutator(rnd=RandomMockup())
        mut.unmutate(
            bam_filename=dir_path + 'input_01.bam',
            diff_file=dir_path + 'input_01.diff',
            out_bam_filename=dir_path + 'output_03.bam'
        )
        # self.assertEqual(15, mut.stat(MutatorCaller.STAT_ALIGNMENT_COUNT))
        # self.assertEqual(0, mut.stat(MutatorCaller.STAT_UNMAPPED_COUNT))
        # self.assertEqual(6, mut.stat(MutatorCaller.STAT_OVERLAPPING_COUNT))
        # self.assertEqual(3, mut.stat(MutatorCaller.STAT_MAX_COVERAGE))
        # self.assertEqual(7, mut.stat(MutatorCaller.STAT_MUT_COUNT))
        # self.assertEqual(3, mut.stat(MutatorCaller.STAT_DIFF_COUNT))
        
        bam2sam(dir_path + 'output_03.bam', dir_path + 'output_03.sam')
        # self.assertEqual(True, filecmp.cmp(dir_path + 'desired.sam', dir_path + 'output_02.sam'))
        
        
        # TODO test case with more references (chromosomes)


if __name__ == '__main__':
    # python3 /usr/local/bin/nosetests -s /data/projects/varlock/scripts/mutator_test.py
    # python3 -m cProfile -s tottime /data/projects/varlock/scripts/mutator.py
    unittest.main()
