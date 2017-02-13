import unittest

from varlock import *
from .random_mockup import RandomMockup


class TestCommon(unittest.TestCase):
    @staticmethod
    def create_fai():
        with pysam.AlignmentFile('resources/common/input.sam', "r") as sam_file:
            return FastaIndex(sam_file)
    
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
        fai = self.create_fai()
        ref = fai.list[21]
        
        self.assertEqual(ref.name, 'chr22')
        self.assertEqual(ref.start, 2829728720)
        self.assertEqual(ref.length, 51304566)
    
    def test_test(self):
        self.assertListEqual([3, 2, 1, 1], count_bases(['A', 'A', 'G', 'T', 'C', 'T', 'A']))
        self.assertListEqual([0, 0, 2, 1], count_bases(['G', 'C', 'G']))
    
    def test_index(self):
        fai = self.create_fai()
        
        self.assertEqual(0, fai.pos2index('chr1', 0))
        self.assertEqual(249250620, fai.pos2index('chr1', 249250620))
        self.assertEqual(249250621, fai.pos2index('chr2', 0))
        self.assertEqual(492449993, fai.pos2index('chr2', 243199372))
        self.assertEqual(492449994, fai.pos2index('chr3', 0))
        
        self.assertTupleEqual(('chr1', 0), fai.index2pos(0))
        self.assertTupleEqual(('chr1', 249250620), fai.index2pos(249250620))
        self.assertTupleEqual(('chr2', 0), fai.index2pos(249250621))
        self.assertTupleEqual(('chr2', 243199372), fai.index2pos(492449993))
        self.assertTupleEqual(('chr3', 0), fai.index2pos(492449994))


if __name__ == '__main__':
    unittest.main()
