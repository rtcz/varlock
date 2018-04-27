import unittest

from varlock_src.cigar import Cigar, NotFoundError


class TestCigar(unittest.TestCase):
    def test_shrink(self):
        self.assertEqual('2M3I1M1D2M', Cigar.shrink('MMIIIMDMM'))
    
    def test_expand(self):
        self.assertEqual('MMIIIMDMM', Cigar.expand('2M3I1M1D2M'))
    
    def test_tuples2exp_str(self):
        self.assertEqual(
            'MMMMDDDMIMM',
            Cigar.tuples2exp_str([
                (Cigar.OP2ID[Cigar.OP_MATCH], 4),
                (Cigar.OP2ID[Cigar.OP_DEL], 3),
                (Cigar.OP2ID[Cigar.OP_MATCH], 1),
                (Cigar.OP2ID[Cigar.OP_INS], 1),
                (Cigar.OP2ID[Cigar.OP_MATCH], 2)
            ])
        )
    
    def test_exp_str2tuples(self):
        self.assertEqual(
            [
                (Cigar.OP2ID[Cigar.OP_MATCH], 4),
                (Cigar.OP2ID[Cigar.OP_DEL], 3),
                (Cigar.OP2ID[Cigar.OP_MATCH], 1),
                (Cigar.OP2ID[Cigar.OP_INS], 1),
                (Cigar.OP2ID[Cigar.OP_MATCH], 2)
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
    
    def test_variant(self):
        self.assertEqual('', Cigar.allele('', ''))
        
        self.assertEqual('M', Cigar.allele('G', 'G'))
        self.assertEqual('D', Cigar.allele('', 'G'))
        self.assertEqual('MI', Cigar.allele('GG', 'G'))
        self.assertEqual('MII', Cigar.allele('GGG', 'G'))
        
        self.assertEqual('MI', Cigar.allele('TT', 'G'))
        self.assertEqual('MD', Cigar.allele('T', 'GG'))
        self.assertEqual('MM', Cigar.allele('TT', 'GG'))
        
        self.assertEqual('MMM', Cigar.allele('GTG', 'GGG'))
        self.assertEqual('MMM', Cigar.allele('GTG', 'GGG'))
        self.assertEqual('MMMMMM', Cigar.allele('GTGTGT', 'GGGGGG'))
        
        self.assertEqual('MDD', Cigar.allele('C', 'CTT'))
        self.assertEqual('MD', Cigar.allele('C', 'CG'))
        self.assertEqual('MM', Cigar.allele('CG', 'CG'))
    
    def test_matching_allele(self):
        self.assertTupleEqual(('A', 'M'), Cigar.matching_allele('TTAAAC', 'MMMMMM', ['A', 'AAAC'], 'A', 2))
        self.assertTupleEqual(('AAAC', 'MIII'), Cigar.matching_allele('TTAAAC', 'MMMIII', ['A', 'AAAC'], 'A', 2))
        
        self.assertTupleEqual(('AA', 'MI'), Cigar.matching_allele('AA', 'MI', ['A', 'AA', 'AC'], 'A', 0))
        self.assertTupleEqual(('AC', 'MI'), Cigar.matching_allele('AC', 'MI', ['A', 'AA', 'AC'], 'A', 0))
        
        # both alleles are covered by CIGAR while one is not covered by sequence
        self.assertTupleEqual(('A', 'MD'), Cigar.matching_allele(
            seq='A',
            exp_cigar='MD',
            alleles=['A', 'AT'],
            ref_allele='AT',
            seq_pos=0,
            cigar_pos=0
        ))
        
        # longest allele can be ruled out due to partial CIGAR mismatch
        self.assertTupleEqual(('C', 'M'), Cigar.matching_allele(
            seq='CA',
            exp_cigar='MM',
            alleles=['C', 'CAT'],
            ref_allele='C',
            seq_pos=0,
            cigar_pos=0
        ))
        
        # longest allele can not be ruled out due to partial CIGAR match
        self.assertRaises(NotFoundError, lambda: Cigar.matching_allele(
            seq='C',
            exp_cigar='M',
            alleles=['C', 'CAT'],
            ref_allele='C',
            seq_pos=0,
            cigar_pos=0
        ))
        
        # ERR015528.25489857 1272715
        
        # TGCTTGTTTAAACCATTGTATTACCAACAGTTGTGTGTTGATTGTGTTTGTTTTTATGAAGTTTCCCAGGGGATGGTGTTCCCTGAGTTTTAGTCTAAGGCATAGATA
        # 101M4I3M
        
        # TGCTTGTTTAAACCATTGTATTACCAACAGTTGTGTGTTGATTGTGTTTGTTTTTATGAAGTTTCCCAGGGGATGGTGTTCCCTGAGTTTTAGTCTAAGGCATA
        # 104M
