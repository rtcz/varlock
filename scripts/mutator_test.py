import unittest
from random import Random

import pysam

from mutator import Mutator, SnvAlignment


class MutatorTest(unittest.TestCase):
    BAM_FILENAME = "/data/projects/varlock/scripts/mapping/1020b_S23_chr22.bam"
    FAI_FILENAME = "/data/projects/varlock/scripts/test/hg19.fa.fai"
    FAI2_FILENAME = "/data/projects/varlock/scripts/test/hs37d5.fa.fai"
    VAC_FILENAME = "/data/projects/varlock/scripts/test/in/chr22.vac"
    OUT_FILENAME = "/data/projects/varlock/scripts/test/1020b_S23_chr22_mut1.bam"
    
    RND_OUT = [0.8444218515250481,
               0.7579544029403025,
               0.420571580830845,
               0.25891675029296335,
               0.5112747213686085,
               0.4049341374504143,
               0.7837985890347726,
               0.30331272607892745,
               0.4765969541523558,
               0.5833820394550312]
    
    @staticmethod
    def build_alignment():
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = "ATGC" * 10
        alignment.reference_start = 1000
        alignment.cigartuples = (
            (Mutator.CIGAR_MATCH, 10),
            (Mutator.CIGAR_DEL, 1),
            (Mutator.CIGAR_EQUAL, 3),
            (Mutator.CIGAR_DIFF, 3),
            (Mutator.CIGAR_INS, 1),
            (Mutator.CIGAR_MATCH, 20)
        )
        return alignment
    
    @staticmethod
    def build_flanked_alignment():
        """
        SAM/BAM files may include extra flanking bases that are not part of the alignment.
        These bases may be the result of the Smith-Waterman or other algorithms,
        which may not require alignments that begin at the first residue or end at the last.
        In addition, extra sequencing adapters, multiplex identifiers,
        and low-quality bases that were not considered for alignment may have been retained.
        """
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = "X" * 10 + "ATGC" * 10 + "X" * 10
        alignment.reference_start = 1000
        alignment.cigartuples = (
            (Mutator.CIGAR_SOFT_CLIP, 10),
            (Mutator.CIGAR_MATCH, 10),
            (Mutator.CIGAR_DEL, 1),
            (Mutator.CIGAR_MATCH, 9),
            (Mutator.CIGAR_INS, 1),
            (Mutator.CIGAR_MATCH, 20),
            (Mutator.CIGAR_SOFT_CLIP, 10)
        )
        return alignment
    
    # def setUp(self):
    #     print("TestUM:setup() before each test method")
    #
    # def tearDown(self):
    #     print("TestUM:teardown() after each test method")
    #
    # @classmethod
    # def setUpClass(cls):
    #     print("setup_class() before any methods in this class")
    #
    # @classmethod
    # def tearDownClass(cls):
    #     print("teardown_class() after any methods in this class")
    
    def test_multi_random(self):
        self.assertEqual(3, Mutator.multi_random([1, 1, 1, 1], Random(0)))
        self.assertEqual(0, Mutator.multi_random([1, 0, 0, 0], Random(0)))
        self.assertEqual(3, Mutator.multi_random([0, 0, 0, 1], Random(0)))
        self.assertEqual(2, Mutator.multi_random([0, 1, 1, 0], Random(0)))
        self.assertEqual(1, Mutator.multi_random([1, 1, 0, 0], Random(0)))
        self.assertEqual(3, Mutator.multi_random([0, 0, 1, 1], Random(0)))
        self.assertEqual(3, Mutator.multi_random([0, 0, 0, 0], Random(0)))
        self.assertEqual(1, Mutator.multi_random([1, 1], Random(0)))
        self.assertEqual(1, Mutator.multi_random([0, 1], Random(0)))
        self.assertEqual(0, Mutator.multi_random([1, 0], Random(0)))
        self.assertEqual(1, Mutator.multi_random([0, 0], Random(0)))
        self.assertEqual(0, Mutator.multi_random([0], Random(0)))
        self.assertEqual(2, Mutator.multi_random([0.8, 0.5, 1.5, 0], Random(0)))
        self.assertEqual(3, Mutator.multi_random([1, 0, 0, 1], Random(0)))
    
    def test_fai_parsing(self):
        fai_list = Mutator.parse_fai(self.FAI_FILENAME)
        fai_line = fai_list[21]
        desired_fail_line = {"name": '22', "length": 51304566, "offset": 2886323449, "linebases": 50}
        self.assertDictEqual(desired_fail_line, fai_line)
        
        fai_list = Mutator.parse_fai(self.FAI2_FILENAME)
        fai_line = fai_list[21]
        desired_fail_line = {"name": '22', "length": 51304566, "offset": 2876892038, "linebases": 60}
        self.assertDictEqual(desired_fail_line, fai_line)
    
    def test_test(self):
        self.assertListEqual([3, 2, 1, 1], Mutator.count_bases(['A', 'A', 'G', 'T', 'C', 'T', 'A']))
        self.assertListEqual([0, 0, 2, 1], Mutator.count_bases(['G', 'C', 'G']))
        
        self.assertEqual('22', Mutator.strip_chr('chr22'))
        self.assertEqual('22', Mutator.strip_chr('22'))
    
    def test_index(self):
        mut = Mutator(self.FAI_FILENAME)
        
        self.assertEqual(0, mut.pos2index('chr1', 0))
        self.assertEqual(0, mut.pos2index('1', 0))
        self.assertEqual(249250620, mut.pos2index('1', 249250620))
        self.assertEqual(249250621, mut.pos2index('2', 0))
        self.assertEqual(492449993, mut.pos2index('2', 243199372))
        self.assertEqual(492449994, mut.pos2index('3', 0))
        
        self.assertTupleEqual(('1', 0), mut.index2pos(0))
        self.assertTupleEqual(('1', 249250620), mut.index2pos(249250620))
        self.assertTupleEqual(('2', 0), mut.index2pos(249250621))
        self.assertTupleEqual(('2', 243199372), mut.index2pos(492449993))
        self.assertTupleEqual(('3', 0), mut.index2pos(492449994))
    
    def test_pileup(self):
        alignment = self.build_alignment()
        self.assertListEqual(['A'], Mutator.get_base_pileup([SnvAlignment(alignment, 0)]))
        self.assertListEqual(['C'], Mutator.get_base_pileup([SnvAlignment(alignment, 3)]))
        self.assertListEqual([], Mutator.get_base_pileup([SnvAlignment(alignment, None)]))
        self.assertListEqual(['C'], Mutator.get_base_pileup([SnvAlignment(alignment, 39)]))
        
        alignment = self.build_flanked_alignment()
        self.assertListEqual(['A'], Mutator.get_base_pileup([SnvAlignment(alignment, 10)]))
        self.assertListEqual(['C'], Mutator.get_base_pileup([SnvAlignment(alignment, 13)]))
    
    def test_base_setting(self):
        # regular alignment
        alignment = self.build_alignment()
        Mutator.set_base(alignment=alignment, pos=0, base='C')
        Mutator.set_base(alignment=alignment, pos=3, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        Mutator.set_base(alignment=alignment, pos=10, base='C')
        Mutator.set_base(alignment=alignment, pos=13, base='A')
        self.assertEqual('C', alignment.query_alignment_sequence[0])
        self.assertEqual('A', alignment.query_alignment_sequence[3])
    
    def test_base_getting(self):
        # regular alignment
        alignment = self.build_alignment()
        self.assertEqual('A', Mutator.get_base(alignment=alignment, pos=0))
        self.assertEqual('C', Mutator.get_base(alignment=alignment, pos=3))
        
        # flanked alignment
        alignment = self.build_flanked_alignment()
        self.assertEqual('A', Mutator.get_base(alignment=alignment, pos=10))
        self.assertEqual('C', Mutator.get_base(alignment=alignment, pos=13))
    
    def test_pos_conversion(self):
        alignment = self.build_alignment()
        self.assertIsNone(Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(9, Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
        
        alignment = self.build_flanked_alignment()
        self.assertIsNone(Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1010))
        self.assertEqual(19, Mutator.ref_pos2seq_pos(alignment=alignment, ref_pos=1009))
    
    def test_cigar(self):
        alignment = self.build_alignment()
        Mutator.check_cigar_str(alignment=alignment, snv_pos=9)
        self.assertRaises(ValueError, lambda: Mutator.check_cigar_str(alignment=alignment, snv_pos=10))
        self.assertRaises(ValueError, lambda: Mutator.check_cigar_str(alignment=alignment, snv_pos=11))
    
    def test_mut_map(self):
        mut_map = Mutator.create_mut_map(alt_ac=[1, 1, 1, 1], ref_ac=[1, 1, 1, 1], rnd=Random(0))
        self.assertDictEqual({'A': 'G', 'T': 'T', 'G': 'C', 'C': 'A', 'N': 'N'}, mut_map)
    
    def test_method(self):
        mut = Mutator(self.FAI_FILENAME, rnd=Random(0))
        # mut.mutate(vac_filename=self.VAC_FILENAME, bam_filename=self.BAM_FILENAME, out_filename=self.OUT_FILENAME)
        # TODO
        pass


if __name__ == '__main__':
    unittest.main()
