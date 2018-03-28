import re
import unittest

from varlock_src.cigar import Cigar


class TestCigar(unittest.TestCase):
    @staticmethod
    def cigar2tuples(cigar: str):
        return [(Cigar.OP_MAP[op], int(count)) for (count, op) in re.findall(re.compile('([0-9]+)([A-Z]+)'), cigar)]
    
    def test_cigarstr2tuples(self):
        self.assertListEqual(
            [(Cigar.OP_MATCH, 66), (Cigar.OP_DEL, 3), (Cigar.OP_INS, 1), (Cigar.OP_MATCH, 41)],
            self.cigar2tuples('66M3D1I41M')
        )
    
    def test_subrange(self):
        self.assertListEqual(self.cigar2tuples('2M'), Cigar.del_subrange(self.cigar2tuples('3M'), 0, 1))
        self.assertListEqual(self.cigar2tuples('1M'), Cigar.del_subrange(self.cigar2tuples('3M'), 0, 2))
        self.assertListEqual([], Cigar.del_subrange(self.cigar2tuples('3M'), 0, 3))
        self.assertListEqual(self.cigar2tuples('2M'), Cigar.del_subrange(self.cigar2tuples('3M'), 1, 2))
        self.assertListEqual(self.cigar2tuples('2M'), Cigar.del_subrange(self.cigar2tuples('3M'), 2, 3))
        self.assertListEqual(self.cigar2tuples('3M'), Cigar.del_subrange(self.cigar2tuples('3M'), 3, 4))
        
        self.assertListEqual(self.cigar2tuples('2M3I'), Cigar.del_subrange(self.cigar2tuples('3M3I'), 0, 1))
        self.assertListEqual(self.cigar2tuples('3M2I'), Cigar.del_subrange(self.cigar2tuples('3M3I'), 5, 6))
        self.assertListEqual(self.cigar2tuples('2M1I'), Cigar.del_subrange(self.cigar2tuples('3M3I'), 2, 5))
        self.assertListEqual(self.cigar2tuples('2M2I'), Cigar.del_subrange(self.cigar2tuples('3M3I'), 2, 4))
        self.assertListEqual(self.cigar2tuples('1M2I'), Cigar.del_subrange(self.cigar2tuples('3M3I'), 1, 4))
        
        self.assertListEqual(self.cigar2tuples('2M2I2M'), Cigar.del_subrange(self.cigar2tuples('2M3I2M'), 3, 4))
        self.assertListEqual(self.cigar2tuples('4M'), Cigar.del_subrange(self.cigar2tuples('3M2I3M'), 2, 6))
        self.assertListEqual(self.cigar2tuples('6M'), Cigar.del_subrange(self.cigar2tuples('3M2I3M'), 3, 5))
        self.assertListEqual(self.cigar2tuples('4M'), Cigar.del_subrange(self.cigar2tuples('3M3I3M'), 1, 6))
        
        self.assertListEqual(self.cigar2tuples('3M'), Cigar.del_subrange(self.cigar2tuples('2M2D2M'), 1, 2))
        self.assertListEqual(self.cigar2tuples('2M2D1M'), Cigar.del_subrange(self.cigar2tuples('2M2D2M'), 2, 3))
        self.assertListEqual(self.cigar2tuples('2M'), Cigar.del_subrange(self.cigar2tuples('2M2D2M'), 0, 2))
        self.assertListEqual(self.cigar2tuples('3M2D2M'), Cigar.del_subrange(self.cigar2tuples('3M2D3M'), 3, 4))
        
        self.assertListEqual(
            self.cigar2tuples('66M3D41M'), Cigar.del_subrange(self.cigar2tuples('66M3D1M2D41M'), 66, 67)
        )
        self.assertListEqual(
            self.cigar2tuples('1I101M'), Cigar.del_subrange(self.cigar2tuples('2M5I101M'), 0, 6)
        )
    
    def test_place_op(self):
        self.assertRaises(AssertionError, lambda: Cigar.place_op([], 0, Cigar.OP_HARD_CLIP, 1))
        self.assertRaises(AssertionError, lambda: Cigar.place_op([], 0, Cigar.OP_SOFT_CLIP, 1))
        
        self.assertListEqual(
            self.cigar2tuples('3M'), Cigar.place_op(self.cigar2tuples('2M'), 0, Cigar.OP_MATCH, 1)
        )
        self.assertListEqual(
            self.cigar2tuples('3M'), Cigar.place_op(self.cigar2tuples('2M'), 2, Cigar.OP_MATCH, 1)
        )
        self.assertListEqual(
            self.cigar2tuples('1I2M'), Cigar.place_op(self.cigar2tuples('2M'), 0, Cigar.OP_INS, 1)
        )
        self.assertListEqual(
            self.cigar2tuples('1M1I1M'), Cigar.place_op(self.cigar2tuples('2M'), 1, Cigar.OP_INS, 1)
        )
        self.assertListEqual(
            self.cigar2tuples('2M1I'), Cigar.place_op(self.cigar2tuples('2M'), 2, Cigar.OP_INS, 1)
        )
        self.assertListEqual(
            self.cigar2tuples('66M3D44M'), Cigar.place_op(self.cigar2tuples('66M3D41M'), 66, Cigar.OP_MATCH, 3)
        )
        self.assertListEqual(
            self.cigar2tuples('1M1I2M'), Cigar.place_op(self.cigar2tuples('1I2M'), 0, Cigar.OP_MATCH, 1)
        )
    
    def test_variant(self):
        self.assertRaises(AssertionError, lambda: Cigar.compute('', ''))
        
        self.assertListEqual([(Cigar.OP_MATCH, 1)], Cigar.compute('G', 'G'))
        self.assertListEqual([(Cigar.OP_DEL, 1)], Cigar.compute('G', ''))
        self.assertListEqual([(Cigar.OP_MATCH, 1), (Cigar.OP_INS, 1)], Cigar.compute('G', 'GG'))
        self.assertListEqual([(Cigar.OP_MATCH, 1), (Cigar.OP_INS, 2)], Cigar.compute('G', 'GGG'))
        
        self.assertListEqual([(Cigar.OP_MATCH, 1)], Cigar.compute('G', 'T'))
        self.assertListEqual([(Cigar.OP_DEL, 1), (Cigar.OP_INS, 2)], Cigar.compute('G', 'TT'))
        self.assertListEqual([(Cigar.OP_DEL, 2), (Cigar.OP_INS, 1)], Cigar.compute('GG', 'T'))
        self.assertListEqual([(Cigar.OP_DEL, 2), (Cigar.OP_INS, 2)], Cigar.compute('GG', 'TT'))
        
        self.assertListEqual([(Cigar.OP_MATCH, 3)], Cigar.compute('GGG', 'GTG'))
        self.assertListEqual([(Cigar.OP_MATCH, 3)], Cigar.compute('GGG', 'GTG'))
        self.assertListEqual([(Cigar.OP_MATCH, 6)], Cigar.compute('GGGGGG', 'GTGTGT'))
        
        self.assertListEqual(self.cigar2tuples('1M2D'), Cigar.compute('CTT', 'C'))
        self.assertListEqual(self.cigar2tuples('3M'), Cigar.compute('CTT', 'CTT'))
