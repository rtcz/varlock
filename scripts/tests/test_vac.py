import filecmp
import io
import unittest

from varlock import Vac, FastaIndex
from varlock.bam import open_bam


class TestVac(unittest.TestCase):
    RESOURCE_PATH = 'tests/resources/vac/'
    INFO_LIST = ['AC=5,3', 'AN=5008', 'NS=2504']
    
    @classmethod
    def vac(cls):
        with open_bam(cls.RESOURCE_PATH + 'input.sam', 'rb') as sam_file:
            return Vac(FastaIndex(sam_file, keep_chr=False))
    
    @classmethod
    def __build_vac_file(cls):
        vac_file = io.BytesIO()
        Vac.write_record(vac_file, 2000000000, (3, 2, 1, 0))
        Vac.write_record(vac_file, 2000000001, (0, 1, 2, 3))
        Vac.write_record(vac_file, 2000000002, (1, 1, 1, 1))
        vac_file.seek(0)
        return vac_file
    
    def test_io(self):
        vac_file = self.__build_vac_file()
        self.assertTupleEqual((2000000000, (3, 2, 1, 0)), Vac.read_record(vac_file))
        self.assertTupleEqual((2000000001, (0, 1, 2, 3)), Vac.read_record(vac_file))
        self.assertTupleEqual((2000000002, (1, 1, 1, 1)), Vac.read_record(vac_file))
        self.assertRaises(EOFError, lambda: Vac.read_record(vac_file))
    
    def test_compact_base_count(self):
        vac = self.vac()
        self.assertListEqual([10000, 1000, 100, 10], vac.compact_base_count([10000, 1000, 100, 10]))
        self.assertListEqual([65535, 6554, 655, 66], vac.compact_base_count([100000, 10000, 1000, 100]))
    
    # AC=5;AF=0.000998403;AN=5008;NS=2504;DP=20274;EAS_AF=0;AMR_AF=0.0014;AFR_AF=0;EUR_AF=0.004;SAS_AF=0;AA=.
    
    def test_parse_ac(self):
        self.assertListEqual([5, 3], Vac.parse_ac(self.INFO_LIST))
    
    def test_parse_an(self):
        self.assertEqual(5008, Vac.parse_an(self.INFO_LIST))
    
    def test_is_snp(self):
        self.assertTrue(Vac.is_snp('A', ['T']))
        self.assertTrue(Vac.is_snp('A', ['T', 'G', 'C']))
        self.assertFalse(Vac.is_snp('A', ['TT', 'G', 'C']))
        self.assertFalse(Vac.is_snp('AA', ['T', 'G', 'C']))
        self.assertFalse(Vac.is_snp('AA', ['TT', 'GG', 'CC']))
        self.assertFalse(Vac.is_snp('AA', ['TT', 'GG', 'CC']))
        self.assertFalse(Vac.is_snp('A', ['.']))
    
    def test_ac_map(self):
        # ref: str, info_ac: list, info_an: int, alt_list: list
        vac = self.vac()
        self.assertDictEqual(
            {'A': 889, 'T': 100, 'G': 10, 'C': 1},
            vac.ac_map('A', ['T', 'G', 'C'], [100, 10, 1], 1000)
        )
        
        # FIXME reference is zero - is this valid case ?
        self.assertDictEqual(
            {'A': 8, 'T': 0, 'G': 0, 'C': 5000},
            vac.ac_map('T', ['A', 'G', 'C'], [8, 0, 5000], 5008)
        )
    
    def test_vcf2vac(self):
        vac = self.vac()
        with open(self.RESOURCE_PATH + 'input.vcf', 'rt') as vcf_file, \
                open(self.RESOURCE_PATH + 'output.vac', 'wb') as out_vac_file:
            vac.vcf2vac(vcf_file, out_vac_file)
        
        vac.vac2text(self.RESOURCE_PATH + 'output.vac', self.RESOURCE_PATH + 'output.vac.txt')
        is_equal = filecmp.cmp(self.RESOURCE_PATH + 'output.vac.txt', self.RESOURCE_PATH + 'desired.vac.txt')
        self.assertEqual(True, is_equal)


if __name__ == '__main__':
    unittest.main()
