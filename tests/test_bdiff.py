import unittest

from varlock.bdiff import BdiffIO


class TestBdiff(unittest.TestCase):
    # SECRET = bytes([255] * Diff.SECRET_SIZE)
    # CHECKSUM = b'0123456789ABCDEF'
    
    _header = {
        BdiffIO.FROM_INDEX: 1000,
        BdiffIO.TO_INDEX: 2000
    }
    
    @classmethod
    def setUpClass(cls):
        bdiff = BdiffIO()
        bdiff._write_snv(1010, 0, ('C', 'A', 'T', 'G'))
        bdiff._write_snv(1020, 1, ('G', 'A', 'T', 'C'))
        bdiff._write_indel(1030, 0, ['AT', 'ATT', 'A'])
        bdiff._write_snv(1040, 2, ('T', 'G', 'A', 'C'))
        bdiff._write_indel(1050, 1, ['T', 'TT'])
        bdiff._write_indel(1060, 2, ['GCG', 'G', 'GCGCG'])
        bdiff._write_snv(1070, 3, ('G', 'C', 'A', 'T'))
        
        cls._bdiff_file = bdiff.file(cls._header)
    
    def test_from_text_file(self):
        # TODO
        pass
    
    def test_to_text_file(self):
        # TODO
        pass
    
    def test_seq_perm(self):
        # TODO
        pass
    
    def test_file_index(self):
        write_io = BdiffIO(index_resolution=3)
        write_io._write_snv(1000, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1010, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1020, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1030, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1040, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1050, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1060, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1070, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1080, 0, ('C', 'A', 'T', 'G'))
        write_io._write_snv(1090, 0, ('C', 'A', 'T', 'G'))
        
        read_io = BdiffIO(write_io.file())
        
        self.assertEqual(1000, read_io.first_index)
        self.assertEqual(1090, read_io.last_index)
        self.assertEqual(10, read_io.snv_count)
        self.assertEqual(0, read_io.indel_count)
        self.assertEqual(True, read_io.is_read_mode)
        self.assertEqual(3, read_io.index_resolution)
        self.assertListEqual([(1020, 14), (1050, 35), (1080, 56)], read_io.file_index)
        
        self.assertEqual(self._index2pos(read_io, 1000), read_io._indexed_pos(0))
        self.assertEqual(self._index2pos(read_io, 1000), read_io._indexed_pos(1000))
        self.assertEqual(self._index2pos(read_io, 1000), read_io._indexed_pos(1019))
        self.assertEqual(self._index2pos(read_io, 1020), read_io._indexed_pos(1020))
        self.assertEqual(self._index2pos(read_io, 1020), read_io._indexed_pos(1021))
        self.assertEqual(self._index2pos(read_io, 1020), read_io._indexed_pos(1049))
        self.assertEqual(self._index2pos(read_io, 1050), read_io._indexed_pos(1050))
        self.assertEqual(self._index2pos(read_io, 1000), read_io._indexed_pos(1010))
    
    def test_io(self):
        bdiff = BdiffIO(self._bdiff_file)
        self.assertTrue(bdiff.is_read_mode)
        self.assertDictEqual(self._header, bdiff.header)
        self.assertEqual(4, bdiff.snv_count)
        self.assertEqual(3, bdiff.indel_count)
        self.assertTupleEqual((1010, 0, ('C', 'A', 'T', 'G')), bdiff._read_record())
        self.assertTupleEqual((1020, 1, ('G', 'A', 'T', 'C')), bdiff._read_record())
        self.assertTupleEqual((1030, 0, ['AT', 'ATT', 'A']), bdiff._read_record())
        self.assertTupleEqual((1040, 2, ('T', 'G', 'A', 'C')), bdiff._read_record())
        self.assertTupleEqual((1050, 1, ['T', 'TT']), bdiff._read_record())
        self.assertTupleEqual((1060, 2, ['GCG', 'G', 'GCGCG']), bdiff._read_record())
        self.assertTupleEqual((1070, 3, ('G', 'C', 'A', 'T')), bdiff._read_record())
        self.assertRaises(EOFError, lambda: bdiff._read_record())
    
    @staticmethod
    def _index2pos(bdiff: BdiffIO, search_index: int):
        """
        Record with search_index must be present.
        :param bdiff:
        :param search_index:
        :return:
        """
        saved_pos = bdiff.tell()
        pos = saved_pos
        while True:
            index, ref_id, alts = bdiff._read_record()
            if index == search_index:
                break
            pos = bdiff.tell()
        
        bdiff.seek(saved_pos)
        return pos
    
    def _record2pos(self, record_id: int):
        """
        Record with record_id must be present.
        :param record_id: 1-based
        :return:
        """
        saved_pos = self._bdiff_file.tell()
        bdiff = BdiffIO(self._bdiff_file)
        
        curr_id = 1
        while curr_id < record_id:
            bdiff._read_record()
            curr_id += 1
        
        record_pos = self._bdiff_file.tell()
        self._bdiff_file.seek(saved_pos)
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
        wbdiff._write_snv(1010, 0, ('C', 'A', 'T', 'G'))
        rbdiff = wbdiff.readcopy()
        self.assertFalse(rbdiff.is_empty())
        
        self.assertEqual(1010, rbdiff.first_index)
        self.assertEqual(1010, rbdiff.last_index)
        
        self.assertEqual(rbdiff.data_offset, rbdiff.tell_index_gte(1000))
        self.assertIsNone(rbdiff.tell_index_gte(2000))
        
        self.assertIsNone(rbdiff.tell_index_lte(1000))
        self.assertEqual(rbdiff.data_offset, rbdiff.tell_index_lte(2000))
    
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
    
    def test_file(self):
        """
        Testing only branch with inner method BdiffIO._file_from_slice.
        """
        bdiff = BdiffIO(self._bdiff_file)
        
        # invalid range
        self.assertRaises(AssertionError, lambda: bdiff.file({BdiffIO.FROM_INDEX: 2000, BdiffIO.TO_INDEX: 1000}, False))
        
        # out of range
        self.assertRaises(IndexError, lambda: bdiff.file({BdiffIO.FROM_INDEX: 1071, BdiffIO.TO_INDEX: 2000}, False))
        self.assertRaises(IndexError, lambda: bdiff.file({BdiffIO.FROM_INDEX: 0, BdiffIO.TO_INDEX: 1009}, False))
        
        # full range
        bdiff_io = BdiffIO(bdiff.file({BdiffIO.FROM_INDEX: 0, BdiffIO.TO_INDEX: 2000}, False))
        self.assertDictEqual({BdiffIO.FROM_INDEX: 0, BdiffIO.TO_INDEX: 2000}, bdiff_io.header)
        self.assertEqual(1010, bdiff_io.first_index)
        self.assertEqual(1070, bdiff_io.last_index)
        self.assertEqual(4, bdiff_io.snv_count)
        self.assertEqual(3, bdiff_io.indel_count)
        
        # exact range
        bdiff_io = BdiffIO(bdiff.file({BdiffIO.FROM_INDEX: 1010, BdiffIO.TO_INDEX: 1070}, False))
        self.assertDictEqual({BdiffIO.FROM_INDEX: 1010, BdiffIO.TO_INDEX: 1070}, bdiff_io.header)
        self.assertEqual(1010, bdiff_io.first_index)
        self.assertEqual(1070, bdiff_io.last_index)
        self.assertEqual(4, bdiff_io.snv_count)
        self.assertEqual(3, bdiff_io.indel_count)
        
        # inner range
        bdiff_io = BdiffIO(bdiff.file({BdiffIO.FROM_INDEX: 1020, BdiffIO.TO_INDEX: 1060}, False))
        self.assertDictEqual({BdiffIO.FROM_INDEX: 1020, BdiffIO.TO_INDEX: 1060}, bdiff_io.header)
        self.assertEqual(1020, bdiff_io.first_index)
        self.assertEqual(1060, bdiff_io.last_index)
        self.assertEqual(2, bdiff_io.snv_count)
        self.assertEqual(3, bdiff_io.indel_count)
        
        # left intersect range
        bdiff_io = BdiffIO(bdiff.file({BdiffIO.FROM_INDEX: 0, BdiffIO.TO_INDEX: 1020}, False))
        self.assertDictEqual({BdiffIO.FROM_INDEX: 0, BdiffIO.TO_INDEX: 1020}, bdiff_io.header)
        self.assertEqual(1010, bdiff_io.first_index)
        self.assertEqual(1020, bdiff_io.last_index)
        self.assertEqual(2, bdiff_io.snv_count)
        self.assertEqual(0, bdiff_io.indel_count)
        
        bdiff_io = BdiffIO(bdiff.file({BdiffIO.FROM_INDEX: 0, BdiffIO.TO_INDEX: 1010}, False))
        self.assertDictEqual({BdiffIO.FROM_INDEX: 0, BdiffIO.TO_INDEX: 1010}, bdiff_io.header)
        self.assertEqual(1010, bdiff_io.first_index)
        self.assertEqual(1010, bdiff_io.last_index)
        self.assertEqual(1, bdiff_io.snv_count)
        self.assertEqual(0, bdiff_io.indel_count)
        
        # right intersect range
        bdiff_io = BdiffIO(bdiff.file({BdiffIO.FROM_INDEX: 1060, BdiffIO.TO_INDEX: 2000}, False))
        self.assertDictEqual({BdiffIO.FROM_INDEX: 1060, BdiffIO.TO_INDEX: 2000}, bdiff_io.header)
        self.assertEqual(1060, bdiff_io.first_index)
        self.assertEqual(1070, bdiff_io.last_index)
        self.assertEqual(1, bdiff_io.snv_count)
        self.assertEqual(1, bdiff_io.indel_count)
        
        bdiff_io = BdiffIO(bdiff.file({BdiffIO.FROM_INDEX: 1070, BdiffIO.TO_INDEX: 2000}, False))
        self.assertDictEqual({BdiffIO.FROM_INDEX: 1070, BdiffIO.TO_INDEX: 2000}, bdiff_io.header)
        self.assertEqual(1070, bdiff_io.first_index)
        self.assertEqual(1070, bdiff_io.last_index)
        self.assertEqual(1, bdiff_io.snv_count)
        self.assertEqual(0, bdiff_io.indel_count)


if __name__ == '__main__':
    unittest.main()
