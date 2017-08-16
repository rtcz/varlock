import unittest

from varlock.bdiff import BdiffIO


class TestBdiff(unittest.TestCase):
    # SECRET = bytes([255] * Diff.SECRET_SIZE)
    # CHECKSUM = b'0123456789ABCDEF'
    RESOURCE_PATH = 'tests/resources/bdiff/'
    
    @classmethod
    def setUpClass(cls):
        bdiff = BdiffIO()
        bdiff.write_snv(1010, ('C', 'A', 'T', 'G'))
        bdiff.write_snv(1020, ('G', 'A', 'T', 'C'))
        bdiff.write_indel(1030, [(10, 'AT'), (5, 'ATT'), (5, 'A')])
        bdiff.write_snv(1040, ('T', 'G', 'A', 'C'))
        bdiff.write_indel(1050, [(5, 'T'), (1, 'TT')])
        bdiff.write_indel(1060, [(10, 'GCG'), (5, 'G'), (1, 'GCGCG')])
        bdiff.write_snv(1070, ('G', 'C', 'A', 'T'))
        cls._bdiff_file = bdiff.file()
    
    def test_from_text_file(self):
        # TODO
        pass
    
    def test_to_text_file(self):
        # TODO
        pass
    
    def test_io(self):
        bdiff = BdiffIO(self._bdiff_file)
        self.assertTrue(bdiff.is_read_mode)
        self.assertDictEqual({}, bdiff.header)
        self.assertEqual(4, bdiff.snv_count)
        self.assertEqual(3, bdiff.indel_count)
        self.assertTupleEqual((1010, ('C', 'A', 'T', 'G')), bdiff.read_record())
        self.assertTupleEqual((1020, ('G', 'A', 'T', 'C')), bdiff.read_record())
        self.assertTupleEqual((1030, [(10, 'AT'), (5, 'ATT'), (5, 'A')]), bdiff.read_record())
        self.assertTupleEqual((1040, ('T', 'G', 'A', 'C')), bdiff.read_record())
        self.assertTupleEqual((1050, [(5, 'T'), (1, 'TT')]), bdiff.read_record())
        self.assertTupleEqual((1060, [(10, 'GCG'), (5, 'G'), (1, 'GCGCG')]), bdiff.read_record())
        self.assertTupleEqual((1070, ('G', 'C', 'A', 'T')), bdiff.read_record())
        self.assertRaises(EOFError, lambda: bdiff.read_record())
    
    def _record2pos(self, record_id: int):
        """
        :param record_id: 1-based
        :return:
        """
        cur_pos = self._bdiff_file.tell()
        bdiff = BdiffIO(self._bdiff_file)
        
        curr_id = 1
        while curr_id < record_id:
            bdiff.read_record()
            curr_id += 1
        
        record_pos = self._bdiff_file.tell()
        self._bdiff_file.seek(cur_pos)
        return record_pos
    
    def test_empty_file(self):
        # empty file
        bdiff = BdiffIO().readcopy()
        self.assertTrue(bdiff.is_empty())
        
        self.assertIsNone(bdiff.first_index)
        self.assertIsNone(bdiff.last_index)
        
        self.assertIsNone(bdiff.tell_index_gte(0))
        self.assertIsNone(bdiff.tell_index_gte(1000))
        
        self.assertIsNone(bdiff.tell_index_lte(0))
        self.assertIsNone(bdiff.tell_index_lte(1000))
    
    def test_single_file(self):
        wbdiff = BdiffIO()
        wbdiff.write_snv(1010, ('C', 'A', 'T', 'G'))
        rbdiff = wbdiff.readcopy()
        self.assertFalse(rbdiff.is_empty())
        
        self.assertEqual(1010, rbdiff.first_index)
        self.assertEqual(1010, rbdiff.last_index)
        
        self.assertEqual(rbdiff.header_size, rbdiff.tell_index_gte(1000))
        self.assertIsNone(rbdiff.tell_index_gte(2000))
        
        self.assertIsNone(rbdiff.tell_index_lte(1000))
        self.assertEqual(rbdiff.header_size, rbdiff.tell_index_lte(2000))
    
    def test_tell_index_gte(self):
        bdiff = BdiffIO(self._bdiff_file)
        
        # one before first
        self.assertEqual(self._record2pos(1), bdiff.tell_index_gte(1009))
        # first
        self.assertEqual(self._record2pos(1), bdiff.tell_index_gte(1010))
        # one after first
        self.assertEqual(self._record2pos(2), bdiff.tell_index_gte(1011))
        # between first and second
        self.assertEqual(self._record2pos(2), bdiff.tell_index_gte(1015))
        # one before second
        self.assertEqual(self._record2pos(2), bdiff.tell_index_gte(1019))
        # second
        self.assertEqual(self._record2pos(2), bdiff.tell_index_gte(1020))
        # one after second
        self.assertEqual(self._record2pos(3), bdiff.tell_index_gte(1021))
        # one before third
        self.assertEqual(self._record2pos(3), bdiff.tell_index_gte(1029))
        # third
        self.assertEqual(self._record2pos(3), bdiff.tell_index_gte(1030))
        # one after third
        self.assertEqual(self._record2pos(4), bdiff.tell_index_gte(1031))
        # one before last
        self.assertEqual(self._record2pos(7), bdiff.tell_index_gte(1069))
        # last
        self.assertEqual(self._record2pos(7), bdiff.tell_index_gte(1070))
        # one after last
        self.assertIsNone(bdiff.tell_index_gte(1071))
    
    def test_tell_index_lte(self):
        bdiff = BdiffIO(self._bdiff_file)
        
        # one before first
        self.assertIsNone(bdiff.tell_index_lte(1008))
        # first
        self.assertEqual(self._record2pos(1), bdiff.tell_index_lte(1010))
        # one after first
        self.assertEqual(self._record2pos(1), bdiff.tell_index_lte(1011))
        # between first and second
        self.assertEqual(self._record2pos(1), bdiff.tell_index_lte(1015))
        # on before second
        self.assertEqual(self._record2pos(1), bdiff.tell_index_lte(1019))
        # second
        self.assertEqual(self._record2pos(2), bdiff.tell_index_lte(1020))
        # one after second
        self.assertEqual(self._record2pos(2), bdiff.tell_index_lte(1021))
        # one before third
        self.assertEqual(self._record2pos(2), bdiff.tell_index_lte(1029))
        # third
        self.assertEqual(self._record2pos(3), bdiff.tell_index_lte(1030))
        # one after third
        self.assertEqual(self._record2pos(3), bdiff.tell_index_lte(1031))
        # one before last
        self.assertEqual(self._record2pos(6), bdiff.tell_index_lte(1069))
        # last
        self.assertEqual(self._record2pos(7), bdiff.tell_index_lte(1070))
        # one after last
        self.assertEqual(self._record2pos(7), bdiff.tell_index_lte(1071))
    
    def _range2pos(self, lower_index: int, upper_index: int):
        return self._record2pos(lower_index), self._record2pos(upper_index)
    
    def test_tell_range(self):
        bdiff = BdiffIO(self._bdiff_file)
        # invalid range
        self.assertRaises(AssertionError, lambda: bdiff.tell_range(2000, 1000))
        
        # outer range
        self.assertRaises(IndexError, lambda: bdiff.tell_range(0, 1000))
        self.assertRaises(IndexError, lambda: bdiff.tell_range(2000, 3000))
        self.assertRaises(IndexError, lambda: bdiff.tell_range(0, 1009))
        self.assertRaises(IndexError, lambda: bdiff.tell_range(1071, 2000))
        
        bdiff = BdiffIO(self._bdiff_file)
        
        # inner range
        self.assertTupleEqual(self._range2pos(1, 7), bdiff.tell_range(1000, 2000))
        self.assertTupleEqual(self._range2pos(1, 7), bdiff.tell_range(1010, 1070))
        self.assertTupleEqual(self._range2pos(1, 2), bdiff.tell_range(1010, 1020))
        self.assertTupleEqual(self._range2pos(2, 3), bdiff.tell_range(1020, 1030))
        self.assertTupleEqual(self._range2pos(3, 4), bdiff.tell_range(1030, 1040))
        self.assertTupleEqual(self._range2pos(1, 1), bdiff.tell_range(1010, 1019))
        self.assertTupleEqual(self._range2pos(1, 1), bdiff.tell_range(1010, 1010))
        self.assertTupleEqual(self._range2pos(2, 6), bdiff.tell_range(1011, 1069))
        
        # right outer range
        self.assertTupleEqual(self._range2pos(1, 7), bdiff.tell_range(1010, 3000))
        self.assertTupleEqual(self._range2pos(2, 7), bdiff.tell_range(1020, 3000))
        self.assertTupleEqual(self._range2pos(3, 7), bdiff.tell_range(1030, 3000))
        
        # left outer range
        self.assertTupleEqual(self._range2pos(1, 1), bdiff.tell_range(1000, 1010))
        self.assertTupleEqual(self._range2pos(1, 2), bdiff.tell_range(1000, 1020))
        self.assertTupleEqual(self._range2pos(1, 3), bdiff.tell_range(1000, 1030))
        
        # def test_header_range(self):
        #     diff_file = io.BytesIO()
        #     self.assertRaises(AssertionError, lambda: Diff.write_header(diff_file, 20, 10, self.CHECKSUM))
        #
        #     diff_file = io.BytesIO()
        #     Diff.write_header(diff_file, 10, 20, self.CHECKSUM)
        #     Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        #     Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        #     Diff.validate_header_range(diff_file)
        #
        #     diff_file = io.BytesIO()
        #     Diff.write_header(diff_file, 10, 10, self.CHECKSUM)
        #     Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        #     Diff.validate_header_range(diff_file)
        #
        #     diff_file = io.BytesIO()
        #     Diff.write_header(diff_file, 11, 20, self.CHECKSUM)
        #     Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        #     Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        #     self.assertRaises(ValueError, lambda: Diff.validate_header_range(diff_file))
        #
        #     diff_file = io.BytesIO()
        #     Diff.write_header(diff_file, 10, 19, self.CHECKSUM)
        #     Diff.write_record(diff_file, 10, ('C', 'A', 'T', 'G'))
        #     Diff.write_record(diff_file, 20, ('G', 'A', 'T', 'C'))
        #     self.assertRaises(ValueError, lambda: Diff.validate_header_range(diff_file))
        
        # def test_truncate(self):
        #     diff_file = self.__build_diff_file()
        #     truncated_file = Diff.truncate(diff_file)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2000000000, 3000000000, self.SECRET),
        #         Diff.read_header(truncated_file)
        #     )
        #     self.assertRaises(EOFError, lambda: Diff.read_record(truncated_file))
        
        # def test_slice(self):
        #     bdiff = BdiffWrapper(self._bdiff_file)
        #
        #     # invalid range
        #     self.assertRaises(AssertionError, lambda: Diff.slice(diff_file, 2000000000, 1000000000))
        #
        #     # out of range
        #     self.assertRaises(IndexError, lambda: Diff.slice(diff_file, 2000000000, 2000000000))
        #     self.assertRaises(IndexError, lambda: Diff.slice(diff_file, 3000000000, 3000000000))
        #
        #     # erase secret key
        #     sliced_file = Diff.slice(diff_file, 2000000000, 3000000000, False)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2000000000, 3000000000, Diff.SECRET_PLACEHOLDER),
        #         Diff.read_header(sliced_file)
        #     )
        #
        #     # exact range
        #     sliced_file = Diff.slice(diff_file, 2830728741, 2830728806)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2830728741, 2830728806, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        #
        #     # outer range
        #     sliced_file = Diff.slice(diff_file, 2000000000, 3000000000)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2000000000, 3000000000, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        #
        #     # inner range
        #     sliced_file = Diff.slice(diff_file, 2830728781, 2830728786)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2830728781, 2830728786, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        #
        #     sliced_file = Diff.slice(diff_file, 2830728780, 2830728790)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2830728780, 2830728790, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728781, ('G', 'A', 'T', 'C')), Diff.read_record(sliced_file))
        #     self.assertTupleEqual((2830728786, ('T', 'G', 'A', 'C')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        #
        #     # right intersect range
        #     sliced_file = Diff.slice(diff_file, 2830728806, 3000000000)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2830728806, 3000000000, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        #
        #     sliced_file = Diff.slice(diff_file, 2830728790, 3000000000)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2830728790, 3000000000, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728806, ('G', 'C', 'A', 'T')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        #
        #     # left intersect range
        #     sliced_file = Diff.slice(diff_file, 2000000000, 2830728741)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2000000000, 2830728741, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))
        #
        #     sliced_file = Diff.slice(diff_file, 2000000000, 2830728750)
        #     self.assertTupleEqual(
        #         (self.CHECKSUM, 2000000000, 2830728750, self.SECRET),
        #         Diff.read_header(sliced_file)
        #     )
        #     self.assertTupleEqual((2830728741, ('C', 'A', 'T', 'G')), Diff.read_record(sliced_file))
        #     self.assertRaises(EOFError, lambda: Diff.read_record(sliced_file))


if __name__ == '__main__':
    unittest.main()
