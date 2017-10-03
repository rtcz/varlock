import unittest
import pysam

from varlock.fasta_index import FastaIndex


class TestFastaIndex(unittest.TestCase):
    # TODO test more methods
    
    @staticmethod
    def create_fai():
        with pysam.AlignmentFile('tests/resources/fasta_index/input.sam', "r") as sam_file:
            return FastaIndex(sam_file.header)
    
    def test_fai_parsing(self):
        fai = self.create_fai()
        ref = fai[1]
        
        self.assertEqual(ref.name, 'chr2')
        self.assertEqual(ref.start, 10000)
        self.assertEqual(ref.length, 10000)
    
    def test_index(self):
        fai = self.create_fai()
        
        self.assertEqual(0, fai.pos2index('chr1', 0))
        self.assertEqual(1000, fai.pos2index('chr1', 1000))
        self.assertEqual(10000, fai.pos2index('chr2', 0))
        self.assertEqual(11000, fai.pos2index('chr2', 1000))
        
        self.assertTupleEqual(('chr1', 0), fai.index2pos(0))
        self.assertTupleEqual(('chr1', 1000), fai.index2pos(1000))
        self.assertTupleEqual(('chr2', 0), fai.index2pos(10000))
        self.assertTupleEqual(('chr2', 1000), fai.index2pos(11000))


if __name__ == '__main__':
    unittest.main()
