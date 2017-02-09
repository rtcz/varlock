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


def count_bases(base_pileup):
    """
    Count allele count in base pileup column.
    :param base_pileup: list of base occurences
    :return: list of DNA base frequencies
    """
    alt_ac = [0] * 4
    for base in base_pileup:
        # skip unknown base
        if base != UNKNOWN_BASE:
            try:
                alt_ac[BASES.index(base)] += 1
            except KeyError:
                raise ValueError("Illegal DNA base %s" % base)
    
    return alt_ac


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


def calc_checksum(file_obj):
    hasher = hashlib.md5()
    for chunk in iter(lambda: file_obj.read(4096), b""):
        hasher.update(chunk)
    
    return hasher.hexdigest()
