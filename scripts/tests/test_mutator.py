import filecmp
import unittest
from random import Random

import pysam

from varlock import *
from .context import *
from .random_mockup import RandomMockup


class TestMutator(unittest.TestCase):
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
        # Mutator.record2bytes()
        pass
    
    def test_mutate(self):
        # convert input from text to binary format
        for i in range(2):
            Vac.text2vac(TEXT_VAC_FILENAMES[i], VAC_FILENAMES[i])
        
        sam2bam(SAM_FILENAME, BAM_FILENAME)
        
        # EOF BAM case
        mut = Mutator(rnd=RandomMockup(), verbose=False)
        with pysam.AlignmentFile(BAM_FILENAME, "rb") as in_bam_file, \
                open(VAC_FILENAMES[0], "rb") as in_vac_file, \
                pysam.AlignmentFile(MUT_BAM_FILENAMES[0], "wb", template=in_bam_file) as out_bam_file, \
                open(DIFF_FILENAMES[0], "wb") as out_diff_file:
            mut.mutate(
                in_vac_file=in_vac_file,
                in_bam_file=in_bam_file,
                out_bam_file=out_bam_file,
                out_diff_file=out_diff_file
            )
        pysam.index(MUT_BAM_FILENAMES[0])
        bam2sam(MUT_BAM_FILENAMES[0], MUT_SAM_FILENAMES[0])
        
        self.assertEqual(True, filecmp.cmp(MUT_SAM_FILENAMES[0], DESIRED_MUT_FILENAMES[0]))
        self.assertEqual(15, mut.alignment_counter)
        self.assertEqual(0, mut.unmapped_counter)
        self.assertEqual(12, mut.overlapping_counter)
        self.assertEqual(5, mut.max_coverage)
        self.assertEqual(5, mut.variant_counter)
        self.assertEqual(12, mut.mut_counter)
        self.assertEqual(4, mut.diff_counter)
        
        # EOF VAC case
        mut = Mutator(rnd=RandomMockup(), verbose=False)
        with pysam.AlignmentFile(BAM_FILENAME, "rb") as in_bam_file, \
                open(VAC_FILENAMES[1], "rb") as in_vac_file, \
                pysam.AlignmentFile(MUT_BAM_FILENAMES[1], "wb", template=in_bam_file) as out_bam_file, \
                open(DIFF_FILENAMES[1], "wb") as out_diff_file:
            mut.mutate(
                in_vac_file=in_vac_file,
                in_bam_file=in_bam_file,
                out_bam_file=out_bam_file,
                out_diff_file=out_diff_file
            )
        bam2sam(MUT_BAM_FILENAMES[1], MUT_SAM_FILENAMES[1])
        pysam.index(MUT_BAM_FILENAMES[0])
        
        self.assertEqual(True, filecmp.cmp(MUT_SAM_FILENAMES[1], DESIRED_MUT_FILENAMES[1]))
        self.assertEqual(15, mut.alignment_counter)
        self.assertEqual(0, mut.unmapped_counter)
        self.assertEqual(7, mut.overlapping_counter)
        self.assertEqual(3, mut.max_coverage)
        self.assertEqual(4, mut.variant_counter)
        self.assertEqual(7, mut.mut_counter)
        self.assertEqual(3, mut.diff_counter)
        
        # TODO test case with more references (chromosomes)
    
    def test_diff(self):
        Diff.diff2text(DIFF_FILENAMES[0], TEXT_DIFF_FILENAMES[0])
        self.assertEqual(True, filecmp.cmp(TEXT_DIFF_FILENAMES[0], DESIRED_DIFF_FILENAMES[0]))
        
        Diff.diff2text(DIFF_FILENAMES[1], TEXT_DIFF_FILENAMES[1])
        self.assertEqual(True, filecmp.cmp(TEXT_DIFF_FILENAMES[1], DESIRED_DIFF_FILENAMES[1]))
        
        with open(DIFF_FILENAMES[0], "rb") as diff_file:
            # diff inside range
            self.assertEqual(0, Diff.seek_index(diff_file, 0, 9999999999))
            # diff after range
            self.assertEqual(20, Diff.seek_index(diff_file, 0, 1111111111))
            # diff before range
            self.assertEqual(20, Diff.seek_index(diff_file, 3333333333, 9999999999))
            # diff covers range
            self.assertEqual(0, Diff.seek_index(diff_file, 2830728741, 2830728781))
            self.assertEqual(5, Diff.seek_index(diff_file, 2830728781, 2830728786))
            self.assertEqual(10, Diff.seek_index(diff_file, 2830728786, 2830728806))
            # diff starts inside range and ends after range
            self.assertEqual(0, Diff.seek_index(diff_file, 2830728741, 9999999999))
            self.assertEqual(5, Diff.seek_index(diff_file, 2830728781, 9999999999))
            self.assertEqual(10, Diff.seek_index(diff_file, 2830728786, 9999999999))
            self.assertEqual(15, Diff.seek_index(diff_file, 2830728806, 9999999999))
            # diff starts before range and ends inside range
            self.assertEqual(0, Diff.seek_index(diff_file, 0, 2830728741))
            self.assertEqual(0, Diff.seek_index(diff_file, 0, 2830728781))
            self.assertEqual(0, Diff.seek_index(diff_file, 0, 2830728786))
            self.assertEqual(0, Diff.seek_index(diff_file, 0, 2830728806))
    
    def test_unmutate(self):
        mut = Mutator(rnd=RandomMockup(), verbose=False)
        with pysam.AlignmentFile(MUT_BAM_FILENAMES[0], "rb") as in_bam_file, \
                open(DIFF_FILENAMES[0], "rb") as in_diff_file, \
                pysam.AlignmentFile(UNMUT_BAM_FILENAMES[0], "wb", template=in_bam_file) as out_bam_file:
            mut.unmutate(
                in_bam_file=in_bam_file,
                in_diff_file=in_diff_file,
                ref_name='chr22',
                start_pos=1000021,
                end_pos=1000061,
                out_bam_file=out_bam_file
            )
        bam2sam(UNMUT_BAM_FILENAMES[0], UNMUT_SAM_FILENAMES[0])
        self.assertEqual(True, filecmp.cmp(UNMUT_SAM_FILENAMES[0], DESIRED_UNMUT_FILENAMES[0]))

        # TODO more tests (ranges etc.)
        # TODO test case with more references (chromosomes)

if __name__ == '__main__':
    # python3 /usr/local/bin/nosetests -s /data/projects/varlock/scripts/mutator_test.py
    # python3 -m cProfile -s tottime /data/projects/varlock/scripts/mutator.py
    unittest.main()
