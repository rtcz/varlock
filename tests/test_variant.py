import unittest

import pysam
import typing
from array import array

from varlock_src.variant import AlignedVariant

from varlock_src.variant import qual_str2array
from varlock_src.variant import qual_array2str

from varlock_src.cigar import Cigar


class TestVariant(unittest.TestCase):

    @staticmethod
    def build_alignment(seq: str, qual_seq: str, cigar: typing.Union[str, None] = None) -> pysam.AlignedSegment:
        """
        SAM/BAM files may include extra flanking bases that are not part of the alignment.
        These bases may be the result of the Smith-Waterman or other algorithms,
        which may not require alignments that begin at the first residue or end at the last.
        In addition, extra sequencing adapters, multiplex identifiers,
        and low-quality bases that were not considered for alignment may have been retained.
        :param seq: query sequence
        :param qual_seq: query quality sequence
        :param cigar: cigar sequence
        """
        alignment = pysam.AlignedSegment()
        alignment.query_sequence = seq
        alignment.reference_start = 1000
        alignment.query_qualities = qual_str2array(qual_seq)
        alignment.query_name = 'test'
        if cigar is None:
            cigar = '%dM' % len(seq)
        alignment.cigarstring = cigar

        return alignment

    def test_qual_convertors(self):
        qual_ar = array('B', [34, 34, 34, 34, 34, 37, 37, 38, 38, 38, 37, 38, 37, 37, 36, 37, 38, 34, 38, 38, 38, 38, 38, 38, 38, 37, 37, 37, 38, 38, 37, 38, 38, 38, 38])
        qual_str = 'CCCCCFFGGGFGFFEFGCGGGGGGGFFFGGFGGGG'

        self.assertEqual(qual_str, qual_array2str(qual_ar))
        self.assertEqual(qual_ar, qual_str2array(qual_str))

    def test_seq(self):
        # insertion:
        align = self.build_alignment("CCAGCTGTGGCAGGGGACAGACATTGTACGGTGGG", "CCCCCFFGGGFGFFEFGCGGGGGGGFFFGGFGGGG")
        av = AlignedVariant(align, 11, 12, 'A')
        av.seq = 'AGTTTT'

        self.assertEqual(av.alignment.query_sequence, "CCAGCTGTGGCAGTTTTGGGGACAGACATTGTACGGTGGG")
        self.assertEqual(qual_array2str(av.alignment.query_qualities), "CCCCCFFGGGFGFFFFFFFEFGCGGGGGGGFFFGGFGGGG")
        self.assertEqual(av.alignment.cigarstring, "12M5I23M")

        # deletion
        align = self.build_alignment("CCAGCTGTGGCAGTTTTGGGGACAGACATTGTACGGTGGG", "CCCCCFFGGGFGFFFFFFFEFGCGGGGGGGFFFGGFGGGG")
        av = AlignedVariant(align, 11, 17, 'AGTTTT')
        av.seq = 'A'

        self.assertEqual(av.alignment.query_sequence, "CCAGCTGTGGCAGGGGACAGACATTGTACGGTGGG")
        self.assertEqual(qual_array2str(av.alignment.query_qualities), "CCCCCFFGGGFGFFEFGCGGGGGGGFFFGGFGGGG")
        self.assertEqual(av.alignment.cigarstring, "12M5D23M")

        # another insertion
        align = self.build_alignment("TCATCTGAATCCACTGGGGAATGGGACGATTTTG", "GGGFGGGGFCGGGGGGGGDGGGGGGGGGGCCCCC")
        av = AlignedVariant(align, 4, 5, 'C')
        av.seq = 'CTTT'

        self.assertEqual(av.alignment.query_sequence, "TCATCTTTTGAATCCACTGGGGAATGGGACGATTTTG")
        self.assertEqual(qual_array2str(av.alignment.query_qualities), "GGGFGGGGGGGFCGGGGGGGGDGGGGGGGGGGCCCCC")
        self.assertEqual(av.alignment.cigarstring, "5M3I29M")


