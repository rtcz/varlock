import unittest

from varlock_src.cigar import Cigar


class TestCigar(unittest.TestCase):
    def test_subrange(self):
        self.assertListEqual([(Cigar.OP_MATCH, 2)], Cigar.del_subrange([(Cigar.OP_MATCH, 3)], 0, 1))
        self.assertListEqual([(Cigar.OP_MATCH, 1)], Cigar.del_subrange([(Cigar.OP_MATCH, 3)], 0, 2))
        self.assertListEqual([], Cigar.del_subrange([(Cigar.OP_MATCH, 3)], 0, 3))
        self.assertListEqual([(Cigar.OP_MATCH, 2)], Cigar.del_subrange([(Cigar.OP_MATCH, 3)], 1, 2))
        self.assertListEqual([(Cigar.OP_MATCH, 2)], Cigar.del_subrange([(Cigar.OP_MATCH, 3)], 2, 3))
        self.assertListEqual([(Cigar.OP_MATCH, 3)], Cigar.del_subrange([(Cigar.OP_MATCH, 3)], 3, 4))
        
        self.assertListEqual(
            [(Cigar.OP_MATCH, 2), (Cigar.OP_INS, 3)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 3)], 0, 1)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 2)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 3)], 5, 6)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 2), (Cigar.OP_INS, 2)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 3)], 2, 4)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 1), (Cigar.OP_INS, 2)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 3)], 1, 4)
        )
        
        self.assertListEqual(
            [(Cigar.OP_MATCH, 2), (Cigar.OP_INS, 2), (Cigar.OP_MATCH, 2)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 2), (Cigar.OP_INS, 3), (Cigar.OP_MATCH, 2)], 3, 4)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 4)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 2), (Cigar.OP_MATCH, 3)], 2, 6)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 6)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 2), (Cigar.OP_MATCH, 3)], 3, 5)
        )
        
        self.assertListEqual(
            [(Cigar.OP_MATCH, 6)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 2), (Cigar.OP_MATCH, 3)], 3, 5)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 4)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_INS, 3), (Cigar.OP_MATCH, 3)], 1, 6)
        )
        
        self.assertListEqual(
            [(Cigar.OP_MATCH, 3)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 2), (Cigar.OP_DEL, 2), (Cigar.OP_MATCH, 2)], 1, 2)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 2)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 2), (Cigar.OP_DEL, 2), (Cigar.OP_MATCH, 2)], 0, 2)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 3), (Cigar.OP_DEL, 2), (Cigar.OP_MATCH, 2)],
            Cigar.del_subrange([(Cigar.OP_MATCH, 3), (Cigar.OP_DEL, 2), (Cigar.OP_MATCH, 3)], 3, 4)
        )
    
    def test_place_op(self):
        self.assertRaises(AssertionError, lambda: Cigar.place_op([], 0, Cigar.OP_HARD_CLIP, 1))
        self.assertRaises(AssertionError, lambda: Cigar.place_op([], 0, Cigar.OP_SOFT_CLIP, 1))
        
        self.assertListEqual(
            [(Cigar.OP_MATCH, 3)],
            Cigar.place_op([(Cigar.OP_MATCH, 2)], 0, Cigar.OP_MATCH, 1)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 3)],
            Cigar.place_op([(Cigar.OP_MATCH, 2)], 1, Cigar.OP_MATCH, 1)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 3)],
            Cigar.place_op([(Cigar.OP_MATCH, 2)], 2, Cigar.OP_MATCH, 1)
        )
        
        self.assertListEqual(
            [(Cigar.OP_INS, 1), (Cigar.OP_MATCH, 2)],
            Cigar.place_op([(Cigar.OP_MATCH, 2)], 0, Cigar.OP_INS, 1)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 1), (Cigar.OP_INS, 1), (Cigar.OP_MATCH, 1)],
            Cigar.place_op([(Cigar.OP_MATCH, 2)], 1, Cigar.OP_INS, 1)
        )
        self.assertListEqual(
            [(Cigar.OP_MATCH, 2), (Cigar.OP_INS, 1)],
            Cigar.place_op([(Cigar.OP_MATCH, 2)], 2, Cigar.OP_INS, 1)
        )
    
    def test_variant(self):
        self.assertRaises(AssertionError, lambda: Cigar.variant('', ''))
        
        self.assertListEqual([(Cigar.OP_MATCH, 1)], Cigar.variant('G', 'G'))
        self.assertListEqual([(Cigar.OP_DEL, 1)], Cigar.variant('G', ''))
        self.assertListEqual([(Cigar.OP_MATCH, 1), (Cigar.OP_INS, 1)], Cigar.variant('G', 'GG'))
        self.assertListEqual([(Cigar.OP_MATCH, 1), (Cigar.OP_INS, 2)], Cigar.variant('G', 'GGG'))
        
        self.assertListEqual([(Cigar.OP_MATCH, 1)], Cigar.variant('G', 'T'))
        self.assertListEqual([(Cigar.OP_DEL, 1), (Cigar.OP_INS, 2)], Cigar.variant('G', 'TT'))
        self.assertListEqual([(Cigar.OP_DEL, 2), (Cigar.OP_INS, 1)], Cigar.variant('GG', 'T'))
        self.assertListEqual([(Cigar.OP_DEL, 2), (Cigar.OP_INS, 2)], Cigar.variant('GG', 'TT'))
        
        self.assertListEqual([(Cigar.OP_MATCH, 3)], Cigar.variant('GGG', 'GTG'))
        self.assertListEqual([(Cigar.OP_MATCH, 3)], Cigar.variant('GGG', 'GTG'))
        self.assertListEqual([(Cigar.OP_MATCH, 6)], Cigar.variant('GGGGGG', 'GTGTGT'))
