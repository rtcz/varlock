import unittest

from varlock_src.cigar import Cigar


class TestCigar(unittest.TestCase):
    # def test_str2tuples(self):
    #     self.assertListEqual(
    #         [
    #             (Cigar.OP_MATCH_ID, 4),
    #             (Cigar.OP_DEL_ID, 3),
    #             (Cigar.OP_MATCH_ID, 1),
    #             (Cigar.OP_INS_ID, 1),
    #             (Cigar.OP_MATCH_ID, 2)
    #         ],
    #         Cigar.str2tuples('4M3D1M1I2M')
    #     )
    #
    # def test_tuples2str(self):
    #     self.assertEqual(
    #         '4M3D1M1I2M',
    #         Cigar.tuples2str([
    #             (Cigar.OP_MATCH_ID, 4),
    #             (Cigar.OP_DEL_ID, 3),
    #             (Cigar.OP_MATCH_ID, 1),
    #             (Cigar.OP_INS_ID, 1),
    #             (Cigar.OP_MATCH_ID, 2)
    #         ])
    #     )
    
    def test_shrink(self):
        self.assertEqual('2M3I1M1D2M', Cigar.shrink('MMIIIMDMM'))
    
    def test_expand(self):
        self.assertEqual('MMIIIMDMM', Cigar.expand('2M3I1M1D2M'))
    
    def test_tuples2exp_str(self):
        self.assertEqual(
            'MMMMDDDMIMM',
            Cigar.tuples2exp_str([
                (Cigar.OP_MATCH_ID, 4),
                (Cigar.OP_DEL_ID, 3),
                (Cigar.OP_MATCH_ID, 1),
                (Cigar.OP_INS_ID, 1),
                (Cigar.OP_MATCH_ID, 2)
            ])
        )
    
    def test_exp_str2tuples(self):
        self.assertEqual(
            [
                (Cigar.OP_MATCH_ID, 4),
                (Cigar.OP_DEL_ID, 3),
                (Cigar.OP_MATCH_ID, 1),
                (Cigar.OP_INS_ID, 1),
                (Cigar.OP_MATCH_ID, 2)
            ],
            Cigar.exp_str2tuples('MMMMDDDMIMM')
        )
    
    def test_mask(self):
        # TODO test INDEL place
        
        self.assertEqual('MMMIMIMM', Cigar.mask('MMMMIMM', 2, 'CG', 'C'))
        
        self.assertEqual('MMMDMM', Cigar.mask('MMMMIMM', 2, 'C', 'CG'))
        self.assertEqual('SSMMMDMM', Cigar.mask('SSMMMMIMM', 4, 'C', 'CG'))
        self.assertEqual('IIMMMDMM', Cigar.mask('IIMMMMIMM', 4, 'C', 'CG'))
        self.assertEqual('MDDMMMDMM', Cigar.mask('MDDMMMMIMM', 3, 'C', 'CG'))
        self.assertEqual('MMMDDMMM', Cigar.mask('MMMMDMMM', 2, 'C', 'CG'))
        
        self.assertEqual('MD', Cigar.mask('', 0, 'C', 'CG'))
        
        self.assertEqual('MD', Cigar.mask('MM', 0, 'C', 'CG'))
        self.assertEqual('MIM', Cigar.mask('MM', 0, 'CG', 'C'))
        self.assertEqual('MM', Cigar.mask('MM', 0, 'CG', 'CG'))
        
        self.assertEqual('MMD', Cigar.mask('MM', 1, 'C', 'CG'))
        self.assertEqual('MMI', Cigar.mask('MM', 1, 'CG', 'C'))
        self.assertEqual('MMM', Cigar.mask('MM', 1, 'CG', 'CG'))
        
        self.assertEqual('MMMD', Cigar.mask('MM', 2, 'C', 'CG'))
        self.assertEqual('MIMD', Cigar.mask('MI', 2, 'C', 'CG'))
        self.assertEqual('MMMI', Cigar.mask('MM', 2, 'CG', 'C'))
        self.assertEqual('MMMM', Cigar.mask('MM', 2, 'CG', 'CG'))
        
        self.assertRaises(ValueError, lambda: Cigar.mask('MD', 2, 'C', 'CG'))
        self.assertRaises(ValueError, lambda: Cigar.mask('MM', 3, 'C', 'CG'))
    
    def test_variant(self):
        self.assertEqual('', Cigar.variant('', ''))
        
        self.assertEqual('M', Cigar.variant('G', 'G'))
        self.assertEqual('D', Cigar.variant('', 'G'))
        self.assertEqual('MI', Cigar.variant('GG', 'G'))
        self.assertEqual('MII', Cigar.variant('GGG', 'G'))
        
        self.assertEqual('MI', Cigar.variant('TT', 'G'))
        self.assertEqual('MD', Cigar.variant('T', 'GG'))
        self.assertEqual('MM', Cigar.variant('TT', 'GG'))
        
        self.assertEqual('MMM', Cigar.variant('GTG', 'GGG'))
        self.assertEqual('MMM', Cigar.variant('GTG', 'GGG'))
        self.assertEqual('MMMMMM', Cigar.variant('GTGTGT', 'GGGGGG'))
        
        self.assertEqual('MDD', Cigar.variant('C', 'CTT'))
        self.assertEqual('MD', Cigar.variant('C', 'CG'))
        self.assertEqual('MM', Cigar.variant('CG', 'CG'))
