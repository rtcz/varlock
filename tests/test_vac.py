import filecmp
import io
import os
import unittest

from varlock.vac import Vac
from varlock.fasta_index import FastaIndex
from varlock.bam import open_bam


class TestVac(unittest.TestCase):
    RESOURCE_PATH = 'tests/resources/vac/'
    
    # TODO test empty and single record file
    
    @classmethod
    def vac(cls):
        with open_bam(cls.RESOURCE_PATH + 'input.sam', 'rb') as sam_file:
            sam_header = sam_file.header
        
        return Vac(FastaIndex(sam_header))
    
    @classmethod
    def __build_empty_vac_file(cls):
        # TODO
        pass
    
    @classmethod
    def __build_vac_file(cls):
        vac_file = io.BytesIO()
        Vac.write_header(vac_file, 3, 3)
        Vac._write_snv_record(vac_file, 2000000000, 0, (3, 2, 1, 0))
        Vac._write_snv_record(vac_file, 2000000002, 1, (0, 1, 2, 3))
        Vac._write_snv_record(vac_file, 2000000004, 2, (1, 1, 1, 1))
        Vac._write_indel_record(vac_file, 2000000001, [(10, 'A'), (1, 'ATCG')])
        Vac._write_indel_record(vac_file, 2000000003, [(10, 'AT'), (1, 'ATCGT')])
        Vac._write_indel_record(vac_file, 2000000005, [(10, 'AAAA'), (10, 'ATCG'), (0, 'A')])
        vac_file.seek(0)
        return vac_file
    
    def test_io(self):
        vac_file = self.__build_vac_file()
        self.assertTupleEqual((3, 3), Vac.read_header(vac_file))
        self.assertTupleEqual((2000000000, 0, (3, 2, 1, 0)), Vac.read_snv_record(vac_file))
        self.assertTupleEqual((2000000002, 1, (0, 1, 2, 3)), Vac.read_snv_record(vac_file))
        self.assertTupleEqual((2000000004, 2, (1, 1, 1, 1)), Vac.read_snv_record(vac_file))
        self.assertTupleEqual((2000000001, [10, 1], ['A', 'ATCG']), Vac.read_indel_record(vac_file))
        self.assertTupleEqual((2000000003, [10, 1], ['AT', 'ATCGT']), Vac.read_indel_record(vac_file))
        self.assertTupleEqual((2000000005, [10, 10, 0], ['AAAA', 'ATCG', 'A']), Vac.read_indel_record(vac_file))
        self.assertRaises(EOFError, lambda: Vac.read_snv_record(vac_file))
    
    def test_compact_base_count(self):
        vac = self.vac()
        self.assertListEqual([10000, 1000, 100, 10], vac.compact_base_count([10000, 1000, 100, 10]))
        self.assertListEqual([65535, 6554, 655, 66], vac.compact_base_count([100000, 10000, 1000, 100]))
    
    def test_parse_ac(self):
        self.assertListEqual([5, 3], Vac.parse_ac(['AC=5,3', 'AN=5008', 'NS=2504']))
    
    def test_parse_an(self):
        self.assertEqual(5008, Vac.parse_an(['AC=5,3', 'AN=5008', 'NS=2504']))
    
    def test_is_snv(self):
        self.assertTrue(Vac.is_snv(['A', 'T']))
        self.assertTrue(Vac.is_snv(['A', 'T', 'G', 'C']))
        
        self.assertFalse(Vac.is_snv(['A', 'TT', 'G', 'C']))
        self.assertFalse(Vac.is_snv(['AA', 'T', 'G', 'C']))
        self.assertFalse(Vac.is_snv(['AA', 'TT', 'GG', 'CC']))
        self.assertFalse(Vac.is_snv(['AA', 'TT', 'GG', 'CC']))
        
        self.assertFalse(Vac.is_snv(['A', '.']))
        
        self.assertFalse(Vac.is_snv(['A', '*']))
        self.assertFalse(Vac.is_snv(['A', '*', 'T']))
        self.assertFalse(Vac.is_snv(['A', 'T', '*']))
        
        self.assertFalse(Vac.is_snv(['N', 'A']))
        self.assertFalse(Vac.is_snv(['A', 'N']))
        self.assertFalse(Vac.is_snv(['A', 'T', 'N']))
        self.assertFalse(Vac.is_snv(['A', 'N', 'T']))
        
        self.assertFalse(Vac.is_snv(['A', '<CN1>', '<CN2>']))
        self.assertFalse(Vac.is_snv(['A', '<INS:ME:ALU>']))
    
    def test_is_indel(self):
        self.assertFalse(Vac.is_indel(['A', 'T']))
        self.assertFalse(Vac.is_indel(['A', 'T', 'G', 'C']))
        
        self.assertTrue(Vac.is_indel(['A', 'TT', 'G', 'C']))
        self.assertTrue(Vac.is_indel(['A', 'T', 'GG', 'CC']))
        self.assertTrue(Vac.is_indel(['AA', 'T', 'G', 'C']))
        self.assertTrue(Vac.is_indel(['AA', 'TT', 'GG', 'CC']))
        
        self.assertFalse(Vac.is_indel(['AA', '.']))
        
        self.assertFalse(Vac.is_indel(['AA', '*']))
        self.assertFalse(Vac.is_indel(['AA', '*', 'T']))
        self.assertFalse(Vac.is_indel(['AA', 'T', '*']))
        self.assertFalse(Vac.is_indel(['A', '*', 'TT']))
        self.assertFalse(Vac.is_indel(['A', 'TT', '*']))
        
        self.assertFalse(Vac.is_indel(['N', 'AA']))
        self.assertFalse(Vac.is_indel(['AA', 'N']))
        self.assertFalse(Vac.is_indel(['AA', 'T', 'N']))
        self.assertFalse(Vac.is_indel(['AA', 'N', 'T']))
        self.assertFalse(Vac.is_indel(['A', 'TT', 'N']))
        self.assertFalse(Vac.is_indel(['A', 'N', 'TT']))
        self.assertFalse(Vac.is_indel(['A', 'TN', 'G']))
        self.assertFalse(Vac.is_indel(['A', 'G', 'TN']))
        
        self.assertFalse(Vac.is_indel(['AA', '<CN1>', '<CN2>']))
        self.assertFalse(Vac.is_indel(['AA', '<INS:ME:ALU>']))
    
    def test_parse_allele_counts(self):
        # ref: str, info_ac: list, info_an: int, alt_list: list
        vac = self.vac()
        self.assertListEqual(
            [889, 100, 10, 1],
            vac.parse_allele_counts('AN=1000;AC=100,10,1')
        )
        self.assertListEqual(
            [5000, 8, 0],
            vac.parse_allele_counts('AN=5008;AC=8,0')
        )
        self.assertListEqual(
            [4008, 1000],
            vac.parse_allele_counts('AN=5008;AC=1000')
        )
    
    def test_vcf2vac(self):
        vac = self.vac()
        vac_filename = self.RESOURCE_PATH + 'output.vac'
        with open(self.RESOURCE_PATH + 'input.vcf', 'rt') as vcf_file, \
                open(vac_filename, 'wb') as vac_file:
            vac.vcf2vac(vcf_file, vac_file)
        
        vac.vac2text(self.RESOURCE_PATH + 'output.vac', self.RESOURCE_PATH + 'output.vac.txt')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'output.vac.txt', self.RESOURCE_PATH + 'desired.vac.txt')
        self.assertEqual(True, is_equal)
        self.assertFalse(os.path.isfile(vac_filename + Vac.SNV_TEMP_EXT))
        self.assertFalse(os.path.isfile(vac_filename + Vac.INDEL_TEMP_EXT))


if __name__ == '__main__':
    unittest.main()
