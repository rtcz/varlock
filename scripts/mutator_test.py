import filecmp
import struct
import unittest
from random import Random

import pysam

from mutator import Mutator, SnvAlignment


class MutatorTest(unittest.TestCase):
    FAI_FILENAME = "/data/projects/varlock/scripts/test/hg19.fa.fai"
    FAI2_FILENAME = "/data/projects/varlock/scripts/test/hs37d5.fa.fai"
    
    TEXT_VAC_FILENAMES = [
        "/data/projects/varlock/scripts/test/input_01.vac.txt",
        "/data/projects/varlock/scripts/test/input_02.vac.txt"
    ]
    VAC_FILENAMES = [
        "/data/projects/varlock/scripts/test/input_01.vac",
        "/data/projects/varlock/scripts/test/input_02.vac",
    ]
    
    IN_SAM_FILENAME = "/data/projects/varlock/scripts/test/input.sam"
    IN_BAM_FILENAME = "/data/projects/varlock/scripts/test/input.bam"
    
    DESIRED_SAM_FILENAMES = [
        "/data/projects/varlock/scripts/test/desired_01.sam",
        "/data/projects/varlock/scripts/test/desired_02.sam"
    ]
    OUT_BAM_FILENAMES = [
        "/data/projects/varlock/scripts/test/output_01.bam",
        "/data/projects/varlock/scripts/test/output_02.bam"
    ]
    OUT_SAM_FILENAMES = [
        "/data/projects/varlock/scripts/test/output_01.sam",
        "/data/projects/varlock/scripts/test/output_02.sam"
    ]
    
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
        reference = fai_list[21]
        self.assertEqual(reference.ref_id, 21)
        self.assertEqual(reference.name, '22')
        self.assertEqual(reference.start, 2829728720)
        self.assertEqual(reference.length, 51304566)
        
        fai_list = Mutator.parse_fai(self.FAI2_FILENAME)
        reference = fai_list[21]
        self.assertEqual(reference.ref_id, 21)
        self.assertEqual(reference.name, '22')
        self.assertEqual(reference.start, 2829728720)
        self.assertEqual(reference.length, 51304566)
    
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
        # EOF BAM case
        mut = Mutator(self.FAI_FILENAME, rnd=RandomMockup(), verbose=False)
        mut.mutate(vac_filename=self.VAC_FILENAMES[0], bam_filename=self.IN_BAM_FILENAME,
                   out_filename=self.OUT_BAM_FILENAMES[0])
        self.bam2sam(self.OUT_BAM_FILENAMES[0], self.OUT_SAM_FILENAMES[0])
        self.assertEqual(True, filecmp.cmp(self.OUT_SAM_FILENAMES[0], self.DESIRED_SAM_FILENAMES[0]))
        self.assertEqual(15, mut.alignment_counter)
        self.assertEqual(0, mut.unmapped_counter)
        self.assertEqual(12, mut.overlapping_counter)
        self.assertEqual(5, mut.max_snv_alignments)
        self.assertEqual(5, mut.snv_counter)
        self.assertEqual(13, mut.mut_counter)
        self.assertEqual(12, mut.ns_mut_counter)
        
        # EOF VAC case
        mut = Mutator(self.FAI_FILENAME, rnd=RandomMockup(), verbose=False)
        mut.mutate(vac_filename=self.VAC_FILENAMES[1], bam_filename=self.IN_BAM_FILENAME,
                   out_filename=self.OUT_BAM_FILENAMES[1])
        self.bam2sam(self.OUT_BAM_FILENAMES[1], self.OUT_SAM_FILENAMES[1])
        self.assertEqual(True, filecmp.cmp(self.OUT_SAM_FILENAMES[1], self.DESIRED_SAM_FILENAMES[1]))
        self.assertEqual(15, mut.alignment_counter)
        self.assertEqual(0, mut.unmapped_counter)
        self.assertEqual(7, mut.overlapping_counter)
        self.assertEqual(3, mut.max_snv_alignments)
        self.assertEqual(4, mut.snv_counter)
        self.assertEqual(8, mut.mut_counter)
        self.assertEqual(7, mut.ns_mut_counter)
    
    @staticmethod
    def text2vac(text_filename, vac_filename):
        with open(text_filename, "r") as text_file, \
                open(vac_filename, "wb") as vac_file:
            for line in text_file:
                data = list(map(int, line.strip().split("\t")))
                vac_file.write(struct.pack('<IHHHH', *data))
    
    @staticmethod
    def vac2list(vac_filename):
        data = []
        with open(vac_filename, "rb") as vac_file:
            while True:
                byte_string = vac_file.read(12)
                if len(byte_string) == 0:
                    break
                data.append(struct.unpack('<IHHHH', byte_string))
        
        return data
    
    @staticmethod
    def sam2bam(sam_filename, bam_filename):
        with pysam.AlignmentFile(sam_filename, "r") as sam_file, \
                pysam.AlignmentFile(bam_filename, "wb", template=sam_file) as bam_file:
            for alignment in sam_file:
                bam_file.write(alignment)
    
    @staticmethod
    def bam2sam(bam_filename, sam_filename):
        with pysam.AlignmentFile(bam_filename, "rb") as bam_file, \
                pysam.AlignmentFile(sam_filename, "w", template=bam_file) as sam_file:
            for alignment in bam_file:
                sam_file.write(alignment)


class RandomMockup:
    def __init__(self):
        pass
    
    @staticmethod
    def random():
        return 0.5
    
    @staticmethod
    def randint(a, b):
        return 0


if __name__ == '__main__':
    # python3 /usr/local/bin/nosetests -s /data/projects/varlock/scripts/mutator_test.py
    # python3 -m cProfile -s tottime /data/projects/varlock/scripts/mutator.py
    
    for i in range(2):
        MutatorTest.text2vac(MutatorTest.TEXT_VAC_FILENAMES[i], MutatorTest.VAC_FILENAMES[i])
        # vac_list = MutatorTest.vac2list(MutatorTest.VAC_FILENAME)
        # print(vac_list)
    
    MutatorTest.sam2bam(MutatorTest.IN_SAM_FILENAME, MutatorTest.IN_BAM_FILENAME)
    unittest.main()
