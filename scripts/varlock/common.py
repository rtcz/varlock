import binascii
import hashlib

import pysam

BASES = ("A", "T", "G", "C")
UNKNOWN_BASE = "N"


def multi_random(p_dist, rnd):
    """
    Draw index from multinomial probability distribution.
    :param p_dist: Array of probabilities. When all probabilities are zero, each outcome has equal probability.
    Each value represents probability of one outcome.
    :param rnd: random number generator
    :return: Always returns index of outcome in p_dist array.
    """
    if len(p_dist) == 1:
        # only one choice
        return 0
    
    if sum(p_dist) == 0:
        # each outcome has equal probability
        return rnd.randint(0, len(p_dist) - 1)
    
    p_level = 0
    rnd_value = rnd.random()
    
    # make relative
    p_value = rnd_value * sum(p_dist)
    for i in range(len(p_dist)):
        p_level += p_dist[i]
        if p_level > p_value:
            # random outcome has been reached
            return i
    
    raise ValueError(
        "sum of probability distribution %d must be greater then the probability value %d" % (sum(p_dist), p_value))


def strip_chr(value):
    """
    Strip 'chr' prefix from string.
    :param value: string ot strip
    :return: stripped string
    """
    if value[:3] == 'chr':
        return value[3:]
    else:
        return value


def sam2bam(sam_filename, bam_filename):
    with pysam.AlignmentFile(sam_filename, "r") as sam_file, \
            pysam.AlignmentFile(bam_filename, "wb", template=sam_file) as bam_file:
        for alignment in sam_file:
            bam_file.write(alignment)


def bam2sam(bam_filename, sam_filename):
    with pysam.AlignmentFile(bam_filename, "rb") as bam_file, \
            pysam.AlignmentFile(sam_filename, "wh", template=bam_file) as sam_file:
        for alignment in bam_file:
            sam_file.write(alignment)


def bin2hex(byte_str):
    return binascii.hexlify(byte_str).decode()


def hex2bin(hex_str):
    return binascii.unhexlify(hex_str)


def calc_checksum(filename):
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    
    return hash_md5.digest()


def ref_pos2seq_pos(alignment, ref_pos):
    """
    Retrieve base position in sequence string at refence position.
    It is assumed that alignment and reference position are on the same reference.
    :param alignment: pysam.AlignedSegment
    :param ref_pos: reference position of base
    :return: 0-based position of base with specified reference position in sequence string
    None if alignment is not mapped at ref_pos (deletion)
    """
    # TODO check reference name, dont assume that reference position and alignment are on same reference
    # TODO optimalize
    seq_pos = None
    for current_seq_pos, current_ref_pos in alignment.get_aligned_pairs(matches_only=False, with_seq=False):
        # search for base in snv position
        if current_ref_pos == ref_pos:
            seq_pos = current_seq_pos
            break
    
    return seq_pos


def get_base_pileup(snv_alignments):
    pileup_col = []
    for snv_alignment in snv_alignments:
        if snv_alignment.pos is not None:
            # alignment is mapped at snv position
            snv_base = get_base(snv_alignment.alignment, snv_alignment.pos)
            pileup_col.append(snv_base)
    
    return pileup_col


def get_base(alignment, pos):
    """
    :param alignment: pysam.AlignedSegment
    :param pos: position in sequence
    :return: base at pos
    """
    return alignment.query_sequence[pos]


def set_base(alignment, pos, base):
    """
    Replace base at SNV position
    :param alignment: pysam.AlignedSegment
    :param pos: position in sequence
    :param base: mutated base letter
    :return: mutated sequence string
    """
    mut_seq = alignment.query_sequence[:pos]
    mut_seq += base
    mut_seq += alignment.query_sequence[pos + 1:]
    alignment.query_sequence = mut_seq
