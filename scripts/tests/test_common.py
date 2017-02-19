import unittest

from varlock import Mutator
from varlock.common import multi_random
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
    
    def test_base_count(self):
        self.assertListEqual([3, 2, 1, 1], Mutator.count_bases(['A', 'A', 'G', 'T', 'C', 'T', 'A']))
        self.assertListEqual([0, 0, 2, 1], Mutator.count_bases(['G', 'C', 'G']))


if __name__ == '__main__':
    unittest.main()
