import unittest

from varlock import *
from .context import *
from .random_mockup import RandomMockup


class TestCommon(unittest.TestCase):
    def test_multi_random(self):
        self.assertEqual(2, multi_random([1, 1, 1, 1], RandomMockup()))
        self.assertEqual(0, multi_random([1, 0, 0, 0], RandomMockup()))
        self.assertEqual(3, multi_random([0, 0, 0, 1], RandomMockup()))
        self.assertEqual(2, multi_random([0, 1, 1, 0], RandomMockup()))
        self.assertEqual(1, multi_random([1, 1, 0, 0], RandomMockup()))
        self.assertEqual(3, multi_random([0, 0, 1, 1], RandomMockup()))
        self.assertEqual(0, multi_random([0, 0, 0, 0], RandomMockup()))
        self.assertEqual(1, multi_random([1, 1], RandomMockup()))
        self.assertEqual(1, multi_random([0, 1], RandomMockup()))
        self.assertEqual(0, multi_random([1, 0], RandomMockup()))
        self.assertEqual(0, multi_random([0, 0], RandomMockup()))
        self.assertEqual(0, multi_random([0], RandomMockup()))
        self.assertEqual(0, multi_random([0.8, 0.05, 0.15, 0], RandomMockup()))
        self.assertEqual(3, multi_random([1, 0, 0, 1], RandomMockup()))
    
    def test_fai_parsing(self):
        with pysam.AlignmentFile(SAM_FILENAME, "r") as sam_file:
            fai_list = get_fai_list(sam_file)
            reference = fai_list[21]
            
            self.assertEqual(reference.name, 'chr22')
            self.assertEqual(reference.start, 2829728720)
            self.assertEqual(reference.length, 51304566)
    
    def test_test(self):
        self.assertListEqual([3, 2, 1, 1], count_bases(['A', 'A', 'G', 'T', 'C', 'T', 'A']))
        self.assertListEqual([0, 0, 2, 1], count_bases(['G', 'C', 'G']))
    
    def test_index(self):
        with pysam.AlignmentFile(SAM_FILENAME, "r") as sam_file:
            fai_list = get_fai_list(sam_file)
            fai_dict = fai_list2dict(fai_list)
            
            self.assertEqual(0, pos2index('chr1', 0, fai_dict))
            self.assertEqual(249250620, pos2index('chr1', 249250620, fai_dict))
            self.assertEqual(249250621, pos2index('chr2', 0, fai_dict))
            self.assertEqual(492449993, pos2index('chr2', 243199372, fai_dict))
            self.assertEqual(492449994, pos2index('chr3', 0, fai_dict))
            
            self.assertTupleEqual(('chr1', 0), index2pos(0, fai_list))
            self.assertTupleEqual(('chr1', 249250620), index2pos(249250620, fai_list))
            self.assertTupleEqual(('chr2', 0), index2pos(249250621, fai_list))
            self.assertTupleEqual(('chr2', 243199372), index2pos(492449993, fai_list))
            self.assertTupleEqual(('chr3', 0), index2pos(492449994, fai_list))


if __name__ == '__main__':
    unittest.main()
