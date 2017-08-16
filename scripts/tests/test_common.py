import unittest
import pysam
from bitarray import bitarray
from bitstring import BitArray

from varlock.cigar import Cigar
import varlock.common as cmn
from varlock.po import SnvAlignment
from tests.random_mockup import RandomMockup
from random import Random


class TestCommon(unittest.TestCase):
    def test_multi_random(self):
        self.assertEqual(2, cmn.multi_random([1, 1, 1, 1], RandomMockup()))
        self.assertEqual(0, cmn.multi_random([1, 0, 0, 0], RandomMockup()))
        self.assertEqual(3, cmn.multi_random([0, 0, 0, 1], RandomMockup()))
        self.assertEqual(2, cmn.multi_random([0, 1, 1, 0], RandomMockup()))
        self.assertEqual(1, cmn.multi_random([1, 1, 0, 0], RandomMockup()))
        self.assertEqual(3, cmn.multi_random([0, 0, 1, 1], RandomMockup()))
        self.assertEqual(0, cmn.multi_random([0, 0, 0, 0], RandomMockup()))
        self.assertEqual(1, cmn.multi_random([1, 1], RandomMockup()))
        self.assertEqual(1, cmn.multi_random([0, 1], RandomMockup()))
        self.assertEqual(0, cmn.multi_random([1, 0], RandomMockup()))
        self.assertEqual(0, cmn.multi_random([0, 0], RandomMockup()))
        self.assertEqual(0, cmn.multi_random([0], RandomMockup()))
        self.assertEqual(0, cmn.multi_random([0.8, 0.05, 0.15, 0], RandomMockup()))
        self.assertEqual(3, cmn.multi_random([1, 0, 0, 1], RandomMockup()))
    
    def test_base_count(self):
        self.assertListEqual([3, 2, 1, 1], cmn.count_bases(['A', 'A', 'G', 'T', 'C', 'T', 'A']))
        self.assertListEqual([0, 0, 2, 1], cmn.count_bases(['G', 'C', 'G']))
    
    @staticmethod
    def build_alignment():
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = "ATGC" * 10
        alignment.reference_start = 1000
        alignment.cigartuples = (
            (Cigar.CIGAR_MATCH, 10),
            (Cigar.CIGAR_DEL, 1),
            (Cigar.CIGAR_EQUAL, 3),
            (Cigar.CIGAR_DIFF, 3),
            (Cigar.CIGAR_INS, 1),
            (Cigar.CIGAR_MATCH, 20)
        )
        return alignment
    
    @staticmethod
    def build_flanked_alignment():
        """
        SAM/BAM files may include extra flanking bases that are not part of the alignment.
        These bases may be the result of the Smith-Waterman or other algorithms,
        which may not require alignments that begin at the first residue or end at the last.
        In addition, extra sequencing adapters, multiplex identifiers,
        and low-quality bases that were not considered for alignment may have been retained.
        """
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = "X" * 10 + "ATGC" * 10 + "X" * 10
        alignment.reference_start = 1000
        alignment.cigartuples = (
            (Cigar.CIGAR_SOFT_CLIP, 10),
            (Cigar.CIGAR_MATCH, 10),
            (Cigar.CIGAR_DEL, 1),
            (Cigar.CIGAR_MATCH, 9),
            (Cigar.CIGAR_INS, 1),
            (Cigar.CIGAR_MATCH, 20),
            (Cigar.CIGAR_SOFT_CLIP, 10)
        )
        return alignment
    
    def test_pileup(self):
        alignment = self.build_alignment()
        self.assertListEqual(['A'], cmn.get_base_pileup([SnvAlignment(alignment, 0)]))
        self.assertListEqual(['C'], cmn.get_base_pileup([SnvAlignment(alignment, 3)]))
        self.assertListEqual([], cmn.get_base_pileup([SnvAlignment(alignment, None)]))
        self.assertListEqual(['C'], cmn.get_base_pileup([SnvAlignment(alignment, 39)]))
        
        alignment = self.build_flanked_alignment()
        self.assertListEqual(['A'], cmn.get_base_pileup([SnvAlignment(alignment, 10)]))
        self.assertListEqual(['C'], cmn.get_base_pileup([SnvAlignment(alignment, 13)]))
    
    def test_base_setting(self):
        # regular alignment
        alignment = self.build_alignment()
        cmn.set_base(alignment=alignment, pos=0, base='C')
        cmn.set_base(alignment=alignment, pos=3, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        cmn.set_base(alignment=alignment, pos=10, base='C')
        cmn.set_base(alignment=alignment, pos=13, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
    
    def test_base_getting(self):
        # regular alignment
        alignment = self.build_alignment()
        self.assertEqual('A', cmn.get_base(alignment=alignment, pos=0))
        self.assertEqual('C', cmn.get_base(alignment=alignment, pos=3))
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        self.assertEqual('A', cmn.get_base(alignment=alignment, pos=10))
        self.assertEqual('C', cmn.get_base(alignment=alignment, pos=13))
    
    def test_pos_conversion(self):
        alignment = self.build_alignment()
        self.assertIsNone(cmn.ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(9, cmn.ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
        
        alignment = self.build_flanked_alignment()
        self.assertIsNone(cmn.ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(19, cmn.ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
    
    def test_cigar(self):
        alignment = self.build_alignment()
        Cigar.validate(alignment=alignment, snv_pos=9)
        self.assertRaises(ValueError, lambda: Cigar.validate(alignment=alignment, snv_pos=10))
        self.assertRaises(ValueError, lambda: Cigar.validate(alignment=alignment, snv_pos=11))
    
    def test_mut_map(self):
        mut_map = cmn.create_mut_map(alt_ac=[1, 1, 1, 1], ref_ac=[1, 1, 1, 1], rnd=Random(0))
        self.assertDictEqual({'A': 'G', 'T': 'T', 'G': 'C', 'C': 'A', 'N': 'N'}, mut_map)
    
    @staticmethod
    def bits2bytes(bit_str: str):
        return bitarray(bit_str).tobytes()
    
    @staticmethod
    def bytes2bits(byte_list: bytes):
        return BitArray(byte_list).bin
    
    def test_stream_cipher(self):
        # ATGC is 00 01 10 11
        # homogenic key
        self.assertEqual('ATGC', cmn.stream_cipher('ATGC', self.bits2bytes('00000000')))
        self.assertEqual('CGTA', cmn.stream_cipher('ATGC', self.bits2bytes('11111111')))
        self.assertEqual('ATGC', cmn.stream_cipher('CGTA', self.bits2bytes('11111111')))
        
        # heterogenic key
        self.assertEqual('CAAC', cmn.stream_cipher('ATGC', self.bits2bytes('11011000')))
        self.assertEqual('ATGC', cmn.stream_cipher('CAAC', self.bits2bytes('11011000')))
        
        # unknown base
        self.assertEqual('CANC', cmn.stream_cipher('ATNC', self.bits2bytes('11011000')))
        self.assertEqual('ATNC', cmn.stream_cipher('CANC', self.bits2bytes('11011000')))
        
        # longer input
        self.assertEqual('CAACCANC', cmn.stream_cipher('ATGCATNC', self.bits2bytes('1101100011011000')))
        self.assertEqual('ATGCATNC', cmn.stream_cipher('CAACCANC', self.bits2bytes('1101100011011000')))
        self.assertEqual('CAACCA', cmn.stream_cipher('ATGCAT', self.bits2bytes('1101100011011000')))
        self.assertEqual('ATGCAT', cmn.stream_cipher('CAACCA', self.bits2bytes('1101100011011000')))
        
        # repeated key
        self.assertEqual('CAACCANC', cmn.stream_cipher('ATGCATNC', self.bits2bytes('11011000')))
        self.assertEqual('ATGCATNC', cmn.stream_cipher('CAACCANC', self.bits2bytes('11011000')))
        self.assertEqual('CAACCA', cmn.stream_cipher('ATGCAT', self.bits2bytes('11011000')))
        self.assertEqual('ATGCAT', cmn.stream_cipher('CAACCA', self.bits2bytes('11011000')))
    
    def test_seq2bytes(self):
        # ATGC is 00 01 10 11
        self.assertEqual('00011011', self.bytes2bits(cmn.seq2bytes('ATGC')))
        self.assertEqual('11100100', self.bytes2bits(cmn.seq2bytes('CGTA')))
        self.assertEqual('00010000', self.bytes2bits(cmn.seq2bytes('AT')))
        self.assertEqual('10110000', self.bytes2bits(cmn.seq2bytes('GC')))
        self.assertEqual('0001101110110000', self.bytes2bits(cmn.seq2bytes('ATGCGC')))
        # self.assertEqual(bitarray('00110111').tobytes(), seq2bytes('ATGCGC'))
    
    def test_bytes2seq(self):
        self.assertEqual('ATGC', cmn.bytes2seq(self.bits2bytes('00011011'), 4))
        self.assertEqual('CGTA', cmn.bytes2seq(self.bits2bytes('11100100'), 4))
        self.assertEqual('AT', cmn.bytes2seq(self.bits2bytes('00010000'), 2))
        self.assertEqual('GC', cmn.bytes2seq(self.bits2bytes('10110000'), 2))
        self.assertEqual('ATGCGC', cmn.bytes2seq(self.bits2bytes('0001101110110000'), 6))


if __name__ == '__main__':
    unittest.main()
