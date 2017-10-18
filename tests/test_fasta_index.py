import unittest
import pysam

from varlock.fasta_index import FastaIndex


class TestFastaIndex(unittest.TestCase):
    def setUp(self):
        with pysam.AlignmentFile('tests/resources/fasta_index/input.sam', "r") as sam_file:
            self._fai = FastaIndex(sam_file.header)
    
    def test_fai_parsing(self):
        ref = self._fai[1]
        
        self.assertEqual(ref.name, 'chr2')
        self.assertEqual(ref.start, 10000)
        self.assertEqual(ref.length, 10000)
    
    def test_pos2index(self):
        self.assertEqual(0, self._fai.pos2index('chr1', 0))
        self.assertEqual(1000, self._fai.pos2index('chr1', 1000))
        self.assertEqual(10000, self._fai.pos2index('chr2', 0))
        self.assertEqual(11000, self._fai.pos2index('chr2', 1000))
    
    def test_index2pos(self):
        self.assertTupleEqual(('chr1', 0), self._fai.index2pos(0))
        self.assertTupleEqual(('chr1', 1000), self._fai.index2pos(1000))
        self.assertTupleEqual(('chr2', 0), self._fai.index2pos(10000))
        self.assertTupleEqual(('chr2', 1000), self._fai.index2pos(11000))
    
    def test_resolve_start_index(self):
        self.assertRaises(ValueError, lambda: self._fai.resolve_start_index(None, 0))
        self.assertEqual(0, self._fai.resolve_start_index(None, None))
        self.assertEqual(0, self._fai.resolve_start_index('chr1', None))
        self.assertEqual(10, self._fai.resolve_start_index('chr1', 10))
        self.assertEqual(10000, self._fai.resolve_start_index('chr2', None))
        self.assertEqual(10010, self._fai.resolve_start_index('chr2', 10))
    
    def test_resolve_end_index(self):
        self.assertRaises(ValueError, lambda: self._fai.resolve_end_index(None, 0))
        self.assertEqual(19999, self._fai.resolve_end_index(None, None))
        self.assertEqual(9999, self._fai.resolve_end_index('chr1', None))
        self.assertEqual(10, self._fai.resolve_end_index('chr1', 10))
        self.assertEqual(19999, self._fai.resolve_end_index('chr2', None))
        self.assertEqual(10010, self._fai.resolve_end_index('chr2', 10))

    def test_index(self):
        self.assertEqual(0, self._fai.first_index())
        self.assertEqual(19999, self._fai.last_index())
        
    def test_ref(self):
        self.assertEqual('chr1', self._fai.first_ref())
        self.assertEqual('chr2', self._fai.last_ref())


if __name__ == '__main__':
    unittest.main()
