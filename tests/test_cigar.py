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
    
    def test_seq_pos2cigar_pos(self):
        self.assertEqual(1, Cigar.seq_pos2cigar_pos('HM', 0))
        self.assertEqual(0, Cigar.seq_pos2cigar_pos('M', 0))
        self.assertEqual(None, Cigar.seq_pos2cigar_pos('HHM', 1))
        self.assertEqual(1, Cigar.seq_pos2cigar_pos('MM', 1))
        self.assertEqual(3, Cigar.seq_pos2cigar_pos('SSMMM', 3))
        self.assertEqual(None, Cigar.seq_pos2cigar_pos('M', 1))
        self.assertEqual(5, Cigar.seq_pos2cigar_pos('MIIMDMM', 4))
    
    # TODO deprecated
    def test_mask(self):
        self.assertEqual('MMMIMIMM', Cigar.mask('MMMMIMM', 2, 'CG', 'C'))
        
        self.assertEqual('MMMDMM', Cigar.mask('MMMMIMM', 2, 'C', 'CG'))
        self.assertEqual('SSMMMDMM', Cigar.mask('SSMMMMIMM', 4, 'C', 'CG'))
        self.assertEqual('HMIM', Cigar.mask('HMM', 0, 'CG', 'C'))
        self.assertEqual('IIMMMDMM', Cigar.mask('IIMMMMIMM', 4, 'C', 'CG'))
        self.assertEqual('MDDMMMDMM', Cigar.mask('MDDMMMMIMM', 3, 'C', 'CG'))
        self.assertEqual('MMMDDMMM', Cigar.mask('MMMMDMMM', 2, 'C', 'CG'))
        
        self.assertRaises(ValueError, lambda: Cigar.mask('', 0, 'C', 'CG'))
        
        self.assertEqual('MD', Cigar.mask('MM', 0, 'C', 'CG'))
        self.assertEqual('MIM', Cigar.mask('MM', 0, 'CG', 'C'))
        self.assertEqual('MM', Cigar.mask('MM', 0, 'CG', 'CG'))
        
        self.assertEqual('MMD', Cigar.mask('MM', 1, 'C', 'CG'))
        self.assertEqual('MMI', Cigar.mask('MM', 1, 'CG', 'C'))
        self.assertEqual('MMM', Cigar.mask('MM', 1, 'CG', 'CG'))
        
        self.assertEqual('MDMD', Cigar.mask('MDM', 1, 'C', 'CG'))
        self.assertEqual('MDMI', Cigar.mask('MDM', 1, 'CG', 'C'))
        self.assertEqual('MDMM', Cigar.mask('MDM', 1, 'CG', 'CG'))
        
        self.assertEqual('MIISS', Cigar.mask('MSS', 0, 'CAT', 'C'))
        self.assertEqual('MIIII', Cigar.mask('MII', 0, 'CAT', 'C'))
        
        self.assertRaises(ValueError, lambda: Cigar.mask('MM', 2, 'CG', 'CG'))
        self.assertRaises(ValueError, lambda: Cigar.mask('MD', 1, 'CG', 'CG'))
        self.assertRaises(ValueError, lambda: Cigar.mask('MI', 2, 'CG', 'CG'))
    
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
    
    def test_matching_allele(self):
        self.assertTupleEqual(('A', 'M'), Cigar._matching_allele('TTAAAC', 'MMMMMM', ['A', 'AAAC'], 0, 2))
        self.assertTupleEqual(('AAAC', 'MIII'), Cigar._matching_allele('TTAAAC', 'MMMIII', ['A', 'AAAC'], 0, 2))
        
        # TODO resolve partial match
        # self.assertTupleEqual(('AAAC', 'MIII'), Cigar._matching_allele('TTAAA', 'MMMII', ['A', 'AAAC'], 0, 2))
        # temporary solution
        self.assertRaises(ValueError, lambda: Cigar._matching_allele('TTAAA', 'MMMII', ['A', 'AAAC'], 0, 2))
        
        self.assertTupleEqual(('AA', 'MI'), Cigar._matching_allele('AA', 'MI', ['A', 'AA', 'AC'], 0, 0))
        self.assertTupleEqual(('AC', 'MI'), Cigar._matching_allele('AC', 'MI', ['A', 'AA', 'AC'], 0, 0))
    
    def test_replace_allele(self):
        self.assertTupleEqual(('TTAAAC', 'MMMMMM'), Cigar.replace_allele('TTAAAC', 'MMMMMM', ['A', 'AAAC'], 0, 0, 2))
        self.assertTupleEqual(('TTA', 'MMM'), Cigar.replace_allele('TTAAAC', 'MMMIII', ['A', 'AAAC'], 0, 0, 2))
        
        self.assertTupleEqual(('TTA', 'MMM'), Cigar.replace_allele('TTAAAC', 'MMMIII', ['A', 'AAAC'], 0, 0, 2))
    
    # TODO deprecated
    def test_matching_alleles(self):
        self.assertListEqual(['A'], Cigar.matching_alleles(Cigar.expand('6M'), 2, ['A', 'AAAC'], 0))
        self.assertListEqual(['AAAC'], Cigar.matching_alleles(Cigar.expand('3M3I'), 2, ['A', 'AAAC'], 0))
        
        # temporary solution
        self.assertListEqual(['AAAC'], Cigar.matching_alleles(Cigar.expand('1M3I'), 0, ['A', 'AAAC'], 0))
        self.assertListEqual([], Cigar.matching_alleles(Cigar.expand('1M2I'), 0, ['A', 'AAAC'], 0))
        
        # result has descending order
        self.assertListEqual(['AC', 'AA'], Cigar.matching_alleles(Cigar.expand('1M1I'), 0, ['A', 'AA', 'AC'], 0))
