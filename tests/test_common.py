import unittest
from random import Random

from bitarray import bitarray
from bitstring import BitArray

import varlock_src.common as cmn
from tests.random import RandomMockup
from varlock_src.random import VeryRandom


class TestCommon(unittest.TestCase):
    def test_multi_random(self):
        rnd = VeryRandom(RandomMockup(0.5))
        self.assertEqual(2, rnd.multirand_index([1, 1, 1, 1]))
        self.assertEqual(0, rnd.multirand_index([1, 0, 0, 0]))
        self.assertEqual(3, rnd.multirand_index([0, 0, 0, 1]))
        self.assertEqual(2, rnd.multirand_index([0, 1, 1, 0]))
        self.assertEqual(1, rnd.multirand_index([1, 1, 0, 0]))
        self.assertEqual(3, rnd.multirand_index([0, 0, 1, 1]))
        self.assertEqual(0, rnd.multirand_index([0, 0, 0, 0]))
        self.assertEqual(1, rnd.multirand_index([1, 1]))
        self.assertEqual(1, rnd.multirand_index([0, 1]))
        self.assertEqual(0, rnd.multirand_index([1, 0]))
        self.assertEqual(0, rnd.multirand_index([0, 0]))
        self.assertEqual(0, rnd.multirand_index([0]))
        self.assertEqual(0, rnd.multirand_index([0.8, 0.05, 0.15, 0]))
        self.assertEqual(3, rnd.multirand_index([1, 0, 0, 1]))
    
    def test_base_count(self):
        self.assertListEqual([3, 2, 1, 1], cmn.base_counts(['A', 'A', 'G', 'T', 'C', 'T', 'A']))
        self.assertListEqual([0, 0, 2, 1], cmn.base_counts(['G', 'C', 'G']))
    
    def test_snv_mut_map(self):
        mut_map = cmn.snv_mut_map(private_freqs=[1, 1, 1, 1], public_freqs=[1, 1, 1, 1], rnd=VeryRandom(Random(0)))
        self.assertDictEqual({'A': 'G', 'T': 'T', 'G': 'C', 'C': 'A'}, mut_map)
        
        mut_map = cmn.snv_mut_map(private_freqs=[3, 0, 1, 0], public_freqs=[2, 1, 4, 0],
                                  rnd=VeryRandom(RandomMockup(0.5)))
        self.assertDictEqual({'A': 'G', 'T': 'T', 'G': 'A', 'C': 'C'}, mut_map)
    
    def test_indel_mut_map(self):
        mut_map = cmn.indel_mut_map(
            private_freq_map={'G': 2},
            public_freq_map={'GCG': 1, 'G': 0},
            rnd=VeryRandom(RandomMockup(0.5))
        )
        self.assertDictEqual({'G': 'GCG', 'GCG': 'G'}, mut_map)
        
        mut_map = cmn.indel_mut_map(
            private_freq_map={'AGT': 2, 'A': 1},
            public_freq_map={'AG': 2, 'A': 1, 'AGT': 0},
            rnd=VeryRandom(RandomMockup(0.5))
        )
        self.assertDictEqual({'A': 'A', 'AG': 'AGT', 'AGT': 'AG'}, mut_map)
    
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
    
    def test_byte2base_perm(self):
        self.assertListEqual(['A', 'A', 'A', 'A'], cmn.byte2base_perm(0b00000000))
        self.assertListEqual(['A', 'T', 'G', 'C'], cmn.byte2base_perm(0b11100100))
        self.assertListEqual(['C', 'G', 'T', 'A'], cmn.byte2base_perm(0b00011011))
    
    def test_base2byte_perm(self):
        self.assertEqual(0b00000000, cmn.base_perm2byte(['A', 'A', 'A', 'A']))
        self.assertEqual(0b11100100, cmn.base_perm2byte(['A', 'T', 'G', 'C']))
        self.assertEqual(0b00011011, cmn.base_perm2byte(['C', 'G', 'T', 'A']))


if __name__ == '__main__':
    unittest.main()
