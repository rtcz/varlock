import unittest

from varlock_src.cigar import Cigar


class TestCigar(unittest.TestCase):
    def test_str2tuples(self):
        self.assertListEqual(
            [
                (Cigar.OP_MATCH_ID, 66),
                (Cigar.OP_DEL_ID, 3),
                (Cigar.OP_INS_ID, 1),
                (Cigar.OP_MATCH_ID, 41)
            ],
            Cigar.str2tuples('66M3D1I41M')
        )
    
    def test_tuples2str(self):
        self.assertEqual(
            '66M3D1I41M',
            Cigar.tuples2str([
                (Cigar.OP_MATCH_ID, 66),
                (Cigar.OP_DEL_ID, 3),
                (Cigar.OP_INS_ID, 1),
                (Cigar.OP_MATCH_ID, 41)
            ])
        )
        
    def test_shrink(self):
        self.assertEqual('2M3I1M1D2M', Cigar.shrink('MMIIIMDMM'))
    
    def test_expand(self):
        self.assertEqual('MMIIIMDMM', Cigar.expand('2M3I1M1D2M'))
        
    def test_replace_seq(self):
        self.assertEqual('3M1D2M', Cigar.mask('4M1I2M', 2, 'C', 'CG'))
        self.assertEqual('2S3M1D2M', Cigar.mask('2S4M1I2M', 4, 'C', 'CG'))
        self.assertEqual('2I3M1D2M', Cigar.mask('2I4M1I2M', 4, 'C', 'CG'))
        self.assertEqual('1M2D3M1D2M', Cigar.mask('1M2D4M1I2M', 3, 'C', 'CG'))
        self.assertEqual('61M2D34M', Cigar.mask('62M1D34M', 60, 'C', 'CG'))
        
        # TODO tests
    
    def test_subrange(self):
        self.assertListEqual(Cigar.str2tuples('2M'), Cigar.del_subrange(Cigar.str2tuples('3M'), 0, 1))
        self.assertListEqual(Cigar.str2tuples('1M'), Cigar.del_subrange(Cigar.str2tuples('3M'), 0, 2))
        self.assertListEqual([], Cigar.del_subrange(Cigar.str2tuples('3M'), 0, 3))
        self.assertListEqual(Cigar.str2tuples('2M'), Cigar.del_subrange(Cigar.str2tuples('3M'), 1, 2))
        self.assertListEqual(Cigar.str2tuples('2M'), Cigar.del_subrange(Cigar.str2tuples('3M'), 2, 3))
        self.assertListEqual(Cigar.str2tuples('3M'), Cigar.del_subrange(Cigar.str2tuples('3M'), 3, 4))
        
        self.assertListEqual(Cigar.str2tuples('2M3I'), Cigar.del_subrange(Cigar.str2tuples('3M3I'), 0, 1))
        self.assertListEqual(Cigar.str2tuples('3M2I'), Cigar.del_subrange(Cigar.str2tuples('3M3I'), 5, 6))
        self.assertListEqual(Cigar.str2tuples('2M1I'), Cigar.del_subrange(Cigar.str2tuples('3M3I'), 2, 5))
        self.assertListEqual(Cigar.str2tuples('2M2I'), Cigar.del_subrange(Cigar.str2tuples('3M3I'), 2, 4))
        self.assertListEqual(Cigar.str2tuples('1M2I'), Cigar.del_subrange(Cigar.str2tuples('3M3I'), 1, 4))
        
        self.assertListEqual(Cigar.str2tuples('2M2I2M'), Cigar.del_subrange(Cigar.str2tuples('2M3I2M'), 3, 4))
        self.assertListEqual(Cigar.str2tuples('4M'), Cigar.del_subrange(Cigar.str2tuples('3M2I3M'), 2, 6))
        self.assertListEqual(Cigar.str2tuples('6M'), Cigar.del_subrange(Cigar.str2tuples('3M2I3M'), 3, 5))
        self.assertListEqual(Cigar.str2tuples('4M'), Cigar.del_subrange(Cigar.str2tuples('3M3I3M'), 1, 6))
        
        self.assertListEqual(Cigar.str2tuples('3M'), Cigar.del_subrange(Cigar.str2tuples('2M2D2M'), 1, 2))
        self.assertListEqual(Cigar.str2tuples('2M2D1M'), Cigar.del_subrange(Cigar.str2tuples('2M2D2M'), 2, 3))
        self.assertListEqual(Cigar.str2tuples('2M'), Cigar.del_subrange(Cigar.str2tuples('2M2D2M'), 0, 2))
        self.assertListEqual(Cigar.str2tuples('3M2D2M'), Cigar.del_subrange(Cigar.str2tuples('3M2D3M'), 3, 4))
        
        self.assertListEqual(
            Cigar.str2tuples('66M3D41M'), Cigar.del_subrange(Cigar.str2tuples('66M3D1M2D41M'), 66, 67)
        )
        self.assertListEqual(
            Cigar.str2tuples('1I101M'), Cigar.del_subrange(Cigar.str2tuples('2M5I101M'), 0, 6)
        )
        
        self.assertListEqual(
            Cigar.str2tuples('12S94M'),
            Cigar.del_subrange(Cigar.str2tuples('12S62M1D34M'), 72, 74)
        )
        self.assertListEqual(
            Cigar.str2tuples('12S94M'),
            Cigar.del_subrange(Cigar.str2tuples('12S61M1D34M'), 72, 73)
        )
    
    # TODO
    # def wrapper(self, cigar_f, cigar_str, *args):
    #     return self.tuples2cigar(cigar_f(Cigar.str2tuples(cigar_str) + args))
    
    def test_place_op(self):
        self.assertRaises(AssertionError, lambda: Cigar.place_op([], 0, Cigar.OP_HARD_CLIP_ID, 1))
        self.assertRaises(AssertionError, lambda: Cigar.place_op([], 0, Cigar.OP_SOFT_CLIP_ID, 1))
        
        self.assertListEqual(
            Cigar.str2tuples('3M'), Cigar.place_op(Cigar.str2tuples('2M'), 0, Cigar.OP_MATCH_ID, 1)
        )
        self.assertListEqual(
            Cigar.str2tuples('3M'), Cigar.place_op(Cigar.str2tuples('2M'), 2, Cigar.OP_MATCH_ID, 1)
        )
        self.assertListEqual(
            Cigar.str2tuples('1I2M'), Cigar.place_op(Cigar.str2tuples('2M'), 0, Cigar.OP_INS_ID, 1)
        )
        self.assertListEqual(
            Cigar.str2tuples('1M1I1M'), Cigar.place_op(Cigar.str2tuples('2M'), 1, Cigar.OP_INS_ID, 1)
        )
        self.assertListEqual(
            Cigar.str2tuples('2M1I'), Cigar.place_op(Cigar.str2tuples('2M'), 2, Cigar.OP_INS_ID, 1)
        )
        self.assertListEqual(
            Cigar.str2tuples('66M3D44M'), Cigar.place_op(Cigar.str2tuples('66M3D41M'), 66, Cigar.OP_MATCH_ID, 3)
        )
        self.assertListEqual(
            Cigar.str2tuples('1M1I2M'), Cigar.place_op(Cigar.str2tuples('1I2M'), 0, Cigar.OP_MATCH_ID, 1)
        )
        
        self.assertListEqual(
            Cigar.str2tuples('12S95M'),
            Cigar.place_op(Cigar.str2tuples('12S94M'), 72, Cigar.OP_MATCH_ID, 1)
        )
        self.assertListEqual(
            Cigar.str2tuples('12S61M1D34M'),
            Cigar.place_op(Cigar.str2tuples('12S95M'), 73, Cigar.OP_DEL_ID, 1)
        )
        self.assertListEqual(
            Cigar.str2tuples('12S96M'),
            Cigar.place_op(Cigar.str2tuples('12S94M'), 72, Cigar.OP_MATCH_ID, 2)
        )
    
    def test_variant(self):
        self.assertRaises(AssertionError, lambda: Cigar.compute('', ''))
        
        self.assertListEqual([(Cigar.OP_MATCH_ID, 1)], Cigar.compute('G', 'G'))
        self.assertListEqual([(Cigar.OP_DEL_ID, 1)], Cigar.compute('G', ''))
        self.assertListEqual([(Cigar.OP_MATCH_ID, 1), (Cigar.OP_INS_ID, 1)], Cigar.compute('G', 'GG'))
        self.assertListEqual([(Cigar.OP_MATCH_ID, 1), (Cigar.OP_INS_ID, 2)], Cigar.compute('G', 'GGG'))
        
        self.assertListEqual([(Cigar.OP_MATCH_ID, 1)], Cigar.compute('G', 'T'))
        self.assertListEqual([(Cigar.OP_DEL_ID, 1), (Cigar.OP_INS_ID, 2)], Cigar.compute('G', 'TT'))
        self.assertListEqual([(Cigar.OP_DEL_ID, 2), (Cigar.OP_INS_ID, 1)], Cigar.compute('GG', 'T'))
        self.assertListEqual([(Cigar.OP_DEL_ID, 2), (Cigar.OP_INS_ID, 2)], Cigar.compute('GG', 'TT'))
        
        self.assertListEqual([(Cigar.OP_MATCH_ID, 3)], Cigar.compute('GGG', 'GTG'))
        self.assertListEqual([(Cigar.OP_MATCH_ID, 3)], Cigar.compute('GGG', 'GTG'))
        self.assertListEqual([(Cigar.OP_MATCH_ID, 6)], Cigar.compute('GGGGGG', 'GTGTGT'))
        
        self.assertListEqual(Cigar.str2tuples('1M2D'), Cigar.compute('CTT', 'C'))
        # self.assertListEqual(Cigar.str2tuples('1M1D'), Cigar.compute('CG', 'C'))
        # self.assertListEqual(Cigar.str2tuples('2M'), Cigar.compute('CG', 'CG'))
