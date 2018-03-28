import unittest
from array import array

import pysam

from varlock_src.common import max_match_cigar
from varlock_src.variant import AlignedVariant
from varlock_src.variant import qual_array2str
from varlock_src.variant import qual_str2array


class TestVariant(unittest.TestCase):
    @staticmethod
    def build_alignment(seq: str, qual_seq: str = None, cigar: str = None) -> pysam.AlignedSegment:
        """
        SAM/BAM files may include extra flanking bases that are not part of the alignment.
        These bases may be the result of the Smith-Waterman or other algorithms,
        which may not require alignments that begin at the first residue or end at the last.
        In addition, extra sequencing adapters, multiplex identifiers,
        and low-quality bases that were not considered for alignment may have been retained.
        :param seq: query sequence
        :param qual_seq: query quality sequence
        :param cigar: cigar sequence
        """
        alignment = pysam.AlignedSegment()
        alignment.query_name = 'test'
        alignment.query_sequence = seq
        alignment.reference_start = 1000
        
        if qual_seq is not None:
            alignment.query_qualities = qual_str2array(qual_seq)
        
        if cigar is None:
            cigar = '%dM' % len(seq)
        
        alignment.cigarstring = cigar
        
        return alignment
    
    def test_qual_convertors(self):
        qual_ar = array('B',
                        [34, 34, 34, 34, 34, 37, 37, 38, 38, 38, 37, 38, 37, 37, 36, 37, 38, 34, 38, 38, 38, 38, 38, 38,
                         38, 37, 37, 37, 38, 38, 37, 38, 38, 38, 38])
        qual_str = 'CCCCCFFGGGFGFFEFGCGGGGGGGFFFGGFGGGG'
        
        self.assertEqual(qual_str, qual_array2str(qual_ar))
        self.assertEqual(qual_ar, qual_str2array(qual_str))
    
    def test_seq(self):
        # insertion:
        align = self.build_alignment("CCAGCTGTGGCAGGGGACAGACATTGTACGGTGGG", "CCCCCFFGGGFGFFEFGCGGGGGGGFFFGGFGGGG")
        av = AlignedVariant(align, 11, 12, 'A')
        av.seq = 'AGTTTT'
        
        self.assertEqual(av.alignment.query_sequence, "CCAGCTGTGGCAGTTTTGGGGACAGACATTGTACGGTGGG")
        self.assertEqual(qual_array2str(av.alignment.query_qualities), "CCCCCFFGGGFGFFFFFFFEFGCGGGGGGGFFFGGFGGGG")
        self.assertEqual(av.alignment.cigarstring, "12M5I23M")
        
        # deletion
        align = self.build_alignment("CCAGCTGTGGCAGTTTTGGGGACAGACATTGTACGGTGGG",
                                     "CCCCCFFGGGFGFFFFFFFEFGCGGGGGGGFFFGGFGGGG")
        av = AlignedVariant(align, 11, 17, 'AGTTTT')
        av.seq = 'A'
        
        self.assertEqual(av.alignment.query_sequence, "CCAGCTGTGGCAGGGGACAGACATTGTACGGTGGG")
        self.assertEqual(qual_array2str(av.alignment.query_qualities), "CCCCCFFGGGFGFFEFGCGGGGGGGFFFGGFGGGG")
        self.assertEqual(av.alignment.cigarstring, "12M5D23M")
        
        # another insertion
        align = self.build_alignment("TCATCTGAATCCACTGGGGAATGGGACGATTTTG", "GGGFGGGGFCGGGGGGGGDGGGGGGGGGGCCCCC")
        av = AlignedVariant(align, 4, 5, 'C')
        av.seq = 'CTTT'
        
        self.assertEqual(av.alignment.query_sequence, "TCATCTTTTGAATCCACTGGGGAATGGGACGATTTTG")
        self.assertEqual(qual_array2str(av.alignment.query_qualities), "GGGFGGGGGGGFCGGGGGGGGDGGGGGGGGGGCCCCC")
        self.assertEqual(av.alignment.cigarstring, "5M3I29M")
        
        # another one
        align = self.build_alignment("AGTTTCCTTTTTGTCAGAGACAAGGTCTCCCTACG", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCC")
        av = AlignedVariant(align, 26, 28, 'TC')
        av.seq = 'T'
        
        self.assertEqual(av.alignment.query_sequence, "AGTTTCCTTTTTGTCAGAGACAAGGTTCCCTACG")
        self.assertEqual(av.alignment.cigarstring, "27M1D7M")
        
        # double indels
        align = self.build_alignment("CGGACTCTGGGCACTTGGGCCTCAGCTTTCAGAGC", "CC<<CGGGE<<FGGFFG@<FD>FFGGDCEGFCFAF")
        av = AlignedVariant(align, 15, 17, 'TG')
        av.seq = 'T'
        
        self.assertEqual(av.alignment.query_sequence, 'CGGACTCTGGGCACTTGGCCTCAGCTTTCAGAGC')
        self.assertEqual(av.alignment.cigarstring, "16M1D18M")
        
        av = AlignedVariant(av.alignment, 18, 19, 'C')
        av.seq = 'CAG'
        
        self.assertEqual(av.alignment.query_sequence, 'CGGACTCTGGGCACTTGGCAGCTCAGCTTTCAGAGC')
        self.assertEqual(av.alignment.cigarstring, "16M1D3M2I15M")
        
        # and back:
        av = AlignedVariant(av.alignment, 15, 16, 'TG')
        av.seq = 'TG'
        
        self.assertEqual(av.alignment.query_sequence, 'CGGACTCTGGGCACTTGGGCAGCTCAGCTTTCAGAGC')
        self.assertEqual(av.alignment.cigarstring, "20M2I15M")
        
        av = AlignedVariant(av.alignment, 19, 22, 'C')
        av.seq = 'C'
        
        self.assertEqual(av.alignment.query_sequence, 'CGGACTCTGGGCACTTGGGCCTCAGCTTTCAGAGC')
        self.assertEqual(av.alignment.cigarstring, "35M")
    
    def test_max_match_cigar(self):
        alignment = self.build_alignment("AGTTTCCTTTTTGTCAGAGACAAGGTCTCCCTACG")
        alignment.reference_start = 1583399
        
        # TODO
        # print(align.cigarstring)
        self.assertEqual(max_match_cigar(alignment, 25, 1583423, ['T', 'TC'], 'T'), 1)
        
        alignment = self.build_alignment(
            'AGCAGAAGAAATAAAGATTAACTTTGACCTTGATATAACACAAAAGTCTAAAATTCACCATGGATA'
            'CTT'
            'CATGTTATTAATGTTACTCTTTTTAGCTGATTTCTTAGG',
            None, '66M3D1M2D41M')
        
        self.assertEqual(max_match_cigar(alignment, 66, 67, ['CTT', 'C'], 'CTT'), 1)
