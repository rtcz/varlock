import unittest
import pysam

from varlock.fasta_index import FastaIndex


class TestFastaIndex(unittest.TestCase):
    # TODO test more methods
    
    @staticmethod
    def create_fai():
        with pysam.AlignmentFile('tests/resources/common/input.sam', "r") as sam_file:
            return FastaIndex(sam_file.header)
    
    def test_fai_parsing(self):
        fai = self.create_fai()
        ref = fai[21]
        
        self.assertEqual(ref.name, 'chr22')
        self.assertEqual(ref.start, 2829728720)
        self.assertEqual(ref.length, 51304566)
    
    def test_index(self):
        fai = self.create_fai()
        
        self.assertEqual(0, fai.pos2index('chr1', 0))
        self.assertEqual(249250620, fai.pos2index('chr1', 249250620))
        self.assertEqual(249250621, fai.pos2index('chr2', 0))
        self.assertEqual(492449993, fai.pos2index('chr2', 243199372))
        self.assertEqual(492449994, fai.pos2index('chr3', 0))
        
        self.assertTupleEqual(('chr1', 0), fai.index2pos(0))
        self.assertTupleEqual(('chr1', 249250620), fai.index2pos(249250620))
        self.assertTupleEqual(('chr2', 0), fai.index2pos(249250621))
        self.assertTupleEqual(('chr2', 243199372), fai.index2pos(492449993))
        self.assertTupleEqual(('chr3', 0), fai.index2pos(492449994))


if __name__ == '__main__':
    unittest.main()
