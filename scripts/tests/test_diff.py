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
    
    def test_seek_pos(self):
        diff_file = self.__build_diff_file()
        # position before diff - go to EOF
        self.assertEqual(self.__position(4), Diff.seek_pos(diff_file, 2000000000))
        self.assertEqual(self.__position(0), Diff.seek_pos(diff_file, 2830728741))
        self.assertEqual(self.__position(1), Diff.seek_pos(diff_file, 2830728781))
        self.assertEqual(self.__position(2), Diff.seek_pos(diff_file, 2830728786))
        self.assertEqual(self.__position(3), Diff.seek_pos(diff_file, 2830728806))
        # position after diff - go to EOF
        self.assertEqual(self.__position(4), Diff.seek_pos(diff_file, 3000000000))
    
    @staticmethod
    def __position(record_id):
        return record_id * Diff.RECORD_LENGTH + Diff.HEADER_LENGTH
    
    def test_seek_range_start(self):
        # build diff
        diff_file = self.__build_diff_file()
        # diff inside range
        self.assertEqual(self.__position(0), Diff.seek_range(diff_file, 1000000000, 3000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_range(diff_file, 1000000000, 2000000000))
        self.assertRaises(IndexError, lambda: Diff.seek_range(diff_file, 3000000000, 4000000000))
        # diff covers range
        self.assertEqual(self.__position(0), Diff.seek_range(diff_file, 2830728741, 2830728781))
        self.assertEqual(self.__position(1), Diff.seek_range(diff_file, 2830728781, 2830728786))
        self.assertEqual(self.__position(2), Diff.seek_range(diff_file, 2830728786, 2830728806))
        # diff starts inside range and ends after range
        self.assertEqual(self.__position(0), Diff.seek_range(diff_file, 2830728741, 3000000000))
        self.assertEqual(self.__position(1), Diff.seek_range(diff_file, 2830728781, 3000000000))
        self.assertEqual(self.__position(2), Diff.seek_range(diff_file, 2830728786, 3000000000))
        self.assertEqual(self.__position(3), Diff.seek_range(diff_file, 2830728806, 3000000000))
        # diff starts before range and ends inside range
        self.assertEqual(self.__position(0), Diff.seek_range(diff_file, 0, 2830728741))
        self.assertEqual(self.__position(0), Diff.seek_range(diff_file, 0, 2830728781))
        self.assertEqual(self.__position(0), Diff.seek_range(diff_file, 0, 2830728786))
        self.assertEqual(self.__position(0), Diff.seek_range(diff_file, 0, 2830728806))
    
    def test_header_range(self):
        diff_file = io.BytesIO()
        self.assertRaises(ValueError, lambda: Diff.write_header(diff_file, b'0123456789ABCDEF', 20, 10))
        
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
        diff_file = io.BytesIO()
        sliced_diff = Diff.slice(diff_file, 0, 0)
        
        self.assertRaises(EOFError, lambda: Diff.slice(diff_file, ))
        
        self.assertTupleEqual((b'0123456789ABCDEF', 2000000000, 3000000000), Diff.read_header(sliced_diff))
        self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_diff))
        self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_diff))
        self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_diff))
        self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_diff))
        self.assertRaises(EOFError, lambda: Diff.read_record(sliced_diff))


if __name__ == '__main__':
    unittest.main()
