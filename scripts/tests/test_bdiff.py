import io
import unittest

from varlock.bdiff import BdiffFile
from varlock.diff import Diff


class TestBbiff(unittest.TestCase):
    # SECRET = bytes([255] * Diff.SECRET_SIZE)
    # CHECKSUM = b'0123456789ABCDEF'
    RESOURCE_PATH = 'tests/resources/bdiff/'
    
    @classmethod
    def setUpClass(cls):
        with BdiffFile(cls.RESOURCE_PATH + 'input.bdiff', 'w') as bdiff_file:
            bdiff_file.write_snv(10011, ('C', 'A', 'T', 'G'))
            bdiff_file.write_snv(10021, ('G', 'A', 'T', 'C'))
            bdiff_file.write_indel(10032, ['AT', 'ATT', 'A'])
            bdiff_file.write_snv(10041, ('T', 'G', 'A', 'C'))
            bdiff_file.write_indel(10052, ['T', 'TT'])
            bdiff_file.write_indel(10062, ['GCG', 'G', 'GCGCG'])
            bdiff_file.write_snv(10071, ('G', 'C', 'A', 'T'))
    
    def test_io(self):
        with BdiffFile(self.RESOURCE_PATH + 'input.bdiff', 'r') as bdiff_file:
            self.assertEqual(bdiff_file.filename, self.RESOURCE_PATH + 'input.bdiff')
            self.assertEqual(bdiff_file.mode, 'r')
            self.assertDictEqual({}, bdiff_file.header)
            self.assertEqual(4, bdiff_file.snv_count)
            self.assertEqual(3, bdiff_file.indel_count)
            self.assertTupleEqual((10011, ('C', 'A', 'T', 'G')), bdiff_file.read_snv())
            self.assertTupleEqual((10021, ('G', 'A', 'T', 'C')), bdiff_file.read_snv())
            self.assertTupleEqual((10041, ('T', 'G', 'A', 'C')), bdiff_file.read_snv())
            self.assertTupleEqual((10071, ('G', 'C', 'A', 'T')), bdiff_file.read_snv())
            self.assertRaises(EOFError, lambda: bdiff_file.read_snv())
            
            self.assertTupleEqual((10032, ['AT', 'ATT', 'A']), bdiff_file.read_indel())
            self.assertTupleEqual((10052, ['T', 'TT']), bdiff_file.read_indel())
            self.assertTupleEqual((10062, ['GCG', 'G', 'GCGCG']), bdiff_file.read_indel())
            self.assertRaises(EOFError, lambda: bdiff_file.read_indel())
    
    def test_seek_index(self):
        # TODO
        self.assertEqual(self.rec_offset(0), Diff.seek_closest_index(diff_file, 2830728740))
        self.assertEqual(self.rec_offset(0), Diff.seek_closest_index(diff_file, 2830728741))
        self.assertEqual(self.rec_offset(1), Diff.seek_closest_index(diff_file, 2830728780))
        self.assertEqual(self.rec_offset(1), Diff.seek_closest_index(diff_file, 2830728781))
        self.assertEqual(self.rec_offset(2), Diff.seek_closest_index(diff_file, 2830728782))
        self.assertEqual(self.rec_offset(2), Diff.seek_closest_index(diff_file, 2830728786))
        self.assertEqual(self.rec_offset(3), Diff.seek_closest_index(diff_file, 2830728790))
        self.assertEqual(self.rec_offset(3), Diff.seek_closest_index(diff_file, 2830728806))
        self.assertRaises(IndexError, lambda: Diff.seek_closest_index(diff_file, 2830728807))

        self.assertRaises(IndexError, lambda: Diff.seek_closest_index(diff_file, 2830728740, False))
        self.assertEqual(self.rec_offset(0), Diff.seek_closest_index(diff_file, 2830728741, False))
        self.assertEqual(self.rec_offset(0), Diff.seek_closest_index(diff_file, 2830728780, False))
        self.assertEqual(self.rec_offset(1), Diff.seek_closest_index(diff_file, 2830728781, False))
        self.assertEqual(self.rec_offset(1), Diff.seek_closest_index(diff_file, 2830728782, False))
        self.assertEqual(self.rec_offset(2), Diff.seek_closest_index(diff_file, 2830728786, False))
        self.assertEqual(self.rec_offset(2), Diff.seek_closest_index(diff_file, 2830728790, False))
        self.assertEqual(self.rec_offset(3), Diff.seek_closest_index(diff_file, 2830728806, False))
        self.assertEqual(self.rec_offset(3), Diff.seek_closest_index(diff_file, 2830728807, False))
    
    @staticmethod
    def rec_offset(record_id):
        return record_id * Diff.SNV_RECORD_SIZE + Diff.HEADER_SIZE
    
    def test_seek_range(self):
        # build diff
        diff_file = self.__build_diff_file()
        # invalid range
        self.assertRaises(AssertionError, lambda: Diff.seek_subrange(diff_file, 2000000000, 1000000000))
        
        # outer range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 1000000000, 3000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 1000000000, 2000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 3000000000, 4000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 1000000000, 2830728740))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 2830728807, 3000000000))
        
        # inner range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(2)),
                              Diff.seek_subrange(diff_file, 2830728741, 2830728781))
        self.assertTupleEqual((self.rec_offset(1), self.rec_offset(3)),
                              Diff.seek_subrange(diff_file, 2830728781, 2830728786))
        self.assertTupleEqual((self.rec_offset(2), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 2830728786, 2830728806))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(1)),
                              Diff.seek_subrange(diff_file, 2830728741, 2830728780))
        self.assertTupleEqual((self.rec_offset(1), self.rec_offset(3)),
                              Diff.seek_subrange(diff_file, 2830728780, 2830728790))
        # right intersect range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(1)),
                              Diff.seek_subrange(diff_file, 2830728741, 2830728741))
        self.assertTupleEqual((self.rec_offset(1), self.rec_offset(2)),
                              Diff.seek_subrange(diff_file, 2830728781, 2830728781))
        self.assertTupleEqual((self.rec_offset(3), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 2830728806, 2830728806))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 2830728741, 3000000000))
        self.assertTupleEqual((self.rec_offset(1), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 2830728781, 3000000000))
        self.assertTupleEqual((self.rec_offset(2), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 2830728786, 3000000000))
        self.assertTupleEqual((self.rec_offset(3), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 2830728806, 3000000000))
        # left intersect range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(1)),
                              Diff.seek_subrange(diff_file, 1000000000, 2830728741))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(2)),
                              Diff.seek_subrange(diff_file, 1000000000, 2830728781))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(3)),
                              Diff.seek_subrange(diff_file, 1000000000, 2830728786))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(4)),
                              Diff.seek_subrange(diff_file, 1000000000, 2830728806))
    
    def test_header_range(self):
        diff_file = io.BytesIO()
        self.assertRaises(AssertionError, lambda: Diff.write_header(diff_file, 20, 10, self.CHECKSUM))
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, 10, 20, self.CHECKSUM)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        Diff.validate_header_range(diff_file)
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, 10, 10, self.CHECKSUM)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.validate_header_range(diff_file)
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, 11, 20, self.CHECKSUM)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        self.assertRaises(ValueError, lambda: Diff.validate_header_range(diff_file))
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, 10, 19, self.CHECKSUM)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        self.assertRaises(ValueError, lambda: Diff.validate_header_range(diff_file))
    
    def test_truncate(self):
        diff_file = self.__build_diff_file()
        truncated_file = Diff.truncate(diff_file)
        self.assertTupleEqual(
            (self.CHECKSUM, 2000000000, 3000000000, self.SECRET),
            Diff.read_header(truncated_file)
        )
        self.assertRaises(EOFError, lambda: Diff.read_record(truncated_file))
    
    def test_slice(self):
        diff_file = self.__build_diff_file()
        
        # invalid range
        self.assertRaises(AssertionError, lambda: Diff.slice(diff_file, 2000000000, 1000000000))
        
        # out of range
        self.assertRaises(IndexError, lambda: Diff.slice(diff_file, 2000000000, 2000000000))
        self.assertRaises(IndexError, lambda: Diff.slice(diff_file, 3000000000, 3000000000))
        
        # erase secret key
        sliced_file = Diff.slice(diff_file, 2000000000, 3000000000, False)
        self.assertTupleEqual(
            (self.CHECKSUM, 2000000000, 3000000000, Diff.SECRET_PLACEHOLDER),
            Diff.read_header(sliced_file)
        )
        
        # exact range
        sliced_file = Diff.slice(diff_file, 2830728741, 2830728806)
        self.assertTupleEqual(
            (self.CHECKSUM, 2830728741, 2830728806, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        
        # outer range
        sliced_file = Diff.slice(diff_file, 2000000000, 3000000000)
        self.assertTupleEqual(
            (self.CHECKSUM, 2000000000, 3000000000, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        
        # inner range
        sliced_file = Diff.slice(diff_file, 2830728781, 2830728786)
        self.assertTupleEqual(
            (self.CHECKSUM, 2830728781, 2830728786, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        
        sliced_file = Diff.slice(diff_file, 2830728780, 2830728790)
        self.assertTupleEqual(
            (self.CHECKSUM, 2830728780, 2830728790, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        
        # right intersect range
        sliced_file = Diff.slice(diff_file, 2830728806, 3000000000)
        self.assertTupleEqual(
            (self.CHECKSUM, 2830728806, 3000000000, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        
        sliced_file = Diff.slice(diff_file, 2830728790, 3000000000)
        self.assertTupleEqual(
            (self.CHECKSUM, 2830728790, 3000000000, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        
        # left intersect range
        sliced_file = Diff.slice(diff_file, 2000000000, 2830728741)
        self.assertTupleEqual(
            (self.CHECKSUM, 2000000000, 2830728741, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        
        sliced_file = Diff.slice(diff_file, 2000000000, 2830728750)
        self.assertTupleEqual(
            (self.CHECKSUM, 2000000000, 2830728750, self.SECRET),
            Diff.read_header(sliced_file)
        )
        self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))


if __name__ == '__main__':
    unittest.main()
