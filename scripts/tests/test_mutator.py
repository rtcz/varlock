import filecmp
import unittest
from random import Random

import pysam

from varlock import *
from .context import *
from .random_mockup import RandomMockup


class TestMutator(unittest.TestCase):
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
    
    def test_test(self):
        # TODO
        # Mutator.diff2bytes()
        pass
    
    def test_method(self):
        # EOF BAM case
        fai_list = parse_fai(FAI_FILENAME)
        mut = Mutator(fai_list, rnd=RandomMockup(), verbose=True)
        mut.mutate(
            in_vac_filename=IN_VAC_FILENAMES[0],
            in_bam_filename=IN_BAM_FILENAME,
            out_bam_filename=OUT_BAM_FILENAMES[0],
            out_diff_filename=OUT_DIFF_FILENAMES[0]
        )
        
        bam2sam(OUT_BAM_FILENAMES[0], OUT_SAM_FILENAMES[0])
        self.assertEqual(True, filecmp.cmp(OUT_SAM_FILENAMES[0], DESIRED_BAM_FILENAMES[0]))
        self.assertEqual(15, mut.alignment_counter)
        self.assertEqual(0, mut.unmapped_counter)
        self.assertEqual(12, mut.overlapping_counter)
        self.assertEqual(5, mut.max_snv_alignments)
        self.assertEqual(5, mut.snv_counter)
        self.assertEqual(12, mut.mut_counter)
        self.assertEqual(4, mut.diff_counter)
        
        # EOF VAC case
        fai_list = parse_fai(FAI_FILENAME)
        mut = Mutator(fai_list, rnd=RandomMockup(), verbose=False)
        mut.mutate(
            in_vac_filename=IN_VAC_FILENAMES[1],
            in_bam_filename=IN_BAM_FILENAME,
            out_bam_filename=OUT_BAM_FILENAMES[1],
            out_diff_filename=OUT_DIFF_FILENAMES[1]
        )
        bam2sam(OUT_BAM_FILENAMES[1], OUT_BAM_FILENAMES[1])
        self.assertEqual(True, filecmp.cmp(OUT_BAM_FILENAMES[1], DESIRED_BAM_FILENAMES[1]))
        self.assertEqual(15, mut.alignment_counter)
        self.assertEqual(0, mut.unmapped_counter)
        self.assertEqual(7, mut.overlapping_counter)
        self.assertEqual(3, mut.max_snv_alignments)
        self.assertEqual(4, mut.snv_counter)
        self.assertEqual(7, mut.mut_counter)
        self.assertEqual(3, mut.diff_counter)


if __name__ == '__main__':
    # python3 /usr/local/bin/nosetests -s /data/projects/varlock/scripts/mutator_test.py
    # python3 -m cProfile -s tottime /data/projects/varlock/scripts/mutator.py
    
    for i in range(2):
        Vac.text2vac(TEXT_VAC_FILENAMES[i], IN_VAC_FILENAMES[i])
        # vac_list = MutatorTest.vac2list(MutatorTest.IN_VAC_FILENAME)
        # print(vac_list)
    
    sam2bam(IN_SAM_FILENAME, IN_BAM_FILENAME)
    unittest.main()
