import io
import unittest

from varlock import *


class TestDiff(unittest.TestCase):
    @staticmethod
    def __build_diff_file():
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, b'0123456789ABCDEF', 2000000000, 3000000000)
        Diff.write_record(diff_file, 2830728741, ('C', 'A', 'T', 'G'))
        Diff.write_record(diff_file, 2830728781, ('G', 'A', 'T', 'C'))
        Diff.write_record(diff_file, 2830728786, ('T', 'G', 'A', 'C'))
        Diff.write_record(diff_file, 2830728806, ('G', 'C', 'A', 'T'))
        return diff_file
    
    def test_io(self):
        diff_file = self.__build_diff_file()
        self.assertTupleEqual((b'0123456789ABCDEF', 2000000000, 3000000000), Diff.read_header(diff_file))
        self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(diff_file))
        self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(diff_file))
        self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(diff_file))
        self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(diff_file))
        self.assertRaises(EOFError, lambda: Diff.read_record(diff_file))
    
    def test_seek_index(self):
        diff_file = self.__build_diff_file()
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
        return record_id * Diff.RECORD_LENGTH + Diff.HEADER_LENGTH
    
    def test_seek_range(self):
        # build diff
        diff_file = self.__build_diff_file()
        # diff inside range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(4)), Diff.seek_subrange(diff_file, 1000000000, 3000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 1000000000, 2000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 3000000000, 4000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 1000000000, 2830728740))
        self.assertRaises(IndexError, lambda: Diff.seek_subrange(diff_file, 2830728807, 3000000000))
        # diff covers range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(2)), Diff.seek_subrange(diff_file, 2830728741, 2830728781))
        self.assertTupleEqual((self.rec_offset(1), self.rec_offset(3)), Diff.seek_subrange(diff_file, 2830728781, 2830728786))
        self.assertTupleEqual((self.rec_offset(2), self.rec_offset(4)), Diff.seek_subrange(diff_file, 2830728786, 2830728806))
        # diff starts inside range and ends after range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(1)), Diff.seek_subrange(diff_file, 2830728741, 2830728741))
        self.assertTupleEqual((self.rec_offset(1), self.rec_offset(2)), Diff.seek_subrange(diff_file, 2830728781, 2830728781))
        self.assertTupleEqual((self.rec_offset(3), self.rec_offset(4)), Diff.seek_subrange(diff_file, 2830728806, 2830728806))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(4)), Diff.seek_subrange(diff_file, 2830728741, 3000000000))
        self.assertTupleEqual((self.rec_offset(1), self.rec_offset(4)), Diff.seek_subrange(diff_file, 2830728781, 3000000000))
        self.assertTupleEqual((self.rec_offset(2), self.rec_offset(4)), Diff.seek_subrange(diff_file, 2830728786, 3000000000))
        self.assertTupleEqual((self.rec_offset(3), self.rec_offset(4)), Diff.seek_subrange(diff_file, 2830728806, 3000000000))
        # diff starts before range and ends inside range
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(1)), Diff.seek_subrange(diff_file, 1000000000, 2830728741))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(2)), Diff.seek_subrange(diff_file, 1000000000, 2830728781))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(3)), Diff.seek_subrange(diff_file, 1000000000, 2830728786))
        self.assertTupleEqual((self.rec_offset(0), self.rec_offset(4)), Diff.seek_subrange(diff_file, 1000000000, 2830728806))
    
    def test_header_range(self):
        diff_file = io.BytesIO()
        self.assertRaises(AssertionError, lambda: Diff.write_header(diff_file, b'0123456789ABCDEF', 20, 10))
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, b'0123456789ABCDEF', 10, 20)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        Diff.validate_header_range(diff_file)
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, b'0123456789ABCDEF', 10, 10)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.validate_header_range(diff_file)
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, b'0123456789ABCDEF', 11, 20)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        self.assertRaises(ValueError, lambda: Diff.validate_header_range(diff_file))
        
        diff_file = io.BytesIO()
        Diff.write_header(diff_file, b'0123456789ABCDEF', 10, 19)
        Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        self.assertRaises(ValueError, lambda: Diff.validate_header_range(diff_file))
    
    def test_slice(self):
        # TODO
        diff_file = self.__build_diff_file()
        
        self.assertRaises(EOFError, lambda: Diff.slice(diff_file, 2830728781, 2830728786))
        
        sliced_diff = Diff.slice(diff_file, 2830728781, 2830728786)
        self.assertTupleEqual((b'0123456789ABCDEF', 2830728781, 2830728786), Diff.read_header(sliced_diff))
        self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_diff))
        self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_diff))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_diff))


if __name__ == '__main__':
    unittest.main()
