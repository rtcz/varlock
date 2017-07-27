import numpy as np

import binascii
import gzip
import hashlib
import math
import os

import pysam

BASES = ("A", "T", "G", "C")
UNKNOWN_BASE = "N"

BASE_BITMASKS = [0b11000000, 0b00110000, 0b00001100, 0b00000011]
BASE2BITS = {'A': 0b00, 'T': 0b01, 'G': 0b10, 'C': 0b11}
BITS2BASE = {0b00: 'A', 0b01: 'T', 0b10: 'G', 0b11: 'C'}


def stream_cipher(seq: str, key: bytes):
    """
    :param seq: sequence of DNA bases
    :param key: key stream
    :return: cipher stream in form of DNA bases
    DNA base encryption using stream cipher.
    Input sequence of DNA bases is encrypted as another sequence of same length.
    Each 2 bits of the key stream are used to encrypt one base from the sequence.
    If there is not enough bits of the key stream to encrypt whole sequence, key stream repeats.
    """
    mut_seq = ''
    for i in range(len(seq)):
        # obtain key byte
        secret_byte = key[int(i / 4) % len(key)]
        # calculate padding
        rshift = 6 - (i % 4) * 2
        
        if seq[i] == UNKNOWN_BASE:
            mut_seq += UNKNOWN_BASE
        else:
            try:
                # obtain 2 bits of the key byte
                key_bits = (secret_byte & BASE_BITMASKS[i % 4]) >> rshift
                # apply xor on 2 key bits and 2 base bits
                bits = BASE2BITS[seq[i]] ^ key_bits
            except KeyError:
                raise ValueError("Illegal DNA base %s" % seq[i])
            mut_seq += BITS2BASE[bits]
    
    return mut_seq


def create_mut_map(alt_ac, ref_ac, rnd):
    """
    Creates mutation mapping for base pileup column.
    pileup base -> mutated base
    :param alt_ac: list of DNA bases frequencies
    :param ref_ac:
    :param rnd: random number generator
    :return: dict which is specific mutation mapping
    """
    ref_bases = list(BASES)
    alt_bases = list(BASES)
    ref_ac = list(ref_ac)
    
    # add random value to distinguish tied values
    alt_ac = [ac + rnd.random() for ac in alt_ac]
    
    # init mutation mapping
    mut_map = dict.fromkeys(BASES)
    # unknown base is always mapped to itself
    mut_map[UNKNOWN_BASE] = UNKNOWN_BASE
    # map bases but skip last unmapped base
    for i in range(len(BASES) - 1):
        # draw ref base with multinomial probability
        ref_base_id = multi_random(ref_ac, rnd)
        # draw most abundant base from alt alleles
        alt_base_id = np.argmax(alt_ac)  # type: int
        # add mapping
        mut_map[alt_bases[alt_base_id]] = ref_bases[ref_base_id]
        
        # delete processed items
        del ref_bases[ref_base_id]
        del ref_ac[ref_base_id]
        del alt_bases[alt_base_id]
        del alt_ac[alt_base_id]
    
    # last base mapping is obvious
    mut_map[alt_bases[0]] = ref_bases[0]
    
    return mut_map


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


def count_bases(base_pileup: list):
    """
    Count allele count in base pileup column.
    :param base_pileup: list of DNA bases
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


def bin2hex(byte_str: bytes):
    return binascii.hexlify(byte_str).decode()


def hex2bin(hex_str: str):
    return binascii.unhexlify(hex_str)


def file_size(file_obj):
    old_pos = file_obj.tell()
    file_obj.seek(0, os.SEEK_END)
    size = file_obj.tell()
    file_obj.seek(old_pos, os.SEEK_SET)
    return size


def filename_checksum(filename: str):
    hasher = hashlib.md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    
    return hasher.digest()


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
    # TODO optimalize (matches_only=True is possible)
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


def is_placed_alignment(alignment):
    """
    :param alignment: pysam.AlignedSegment
    :return: True if alignment is placed, it still can be unmapped.
    Unplaced alignment has no reference_name and reference_start.
    """
    return alignment.reference_id != -1


def set_base(alignment, pos, base):
    """
    Replace base at SNV position
    :param alignment: pysam.AlignedSegment
    :param pos: position in sequence
    :param base: mutated base letter
    :return: mutated sequence string
    """
    # TODO query_qualities
    mut_seq = alignment.query_sequence[:pos]
    mut_seq += base
    mut_seq += alignment.query_sequence[pos + 1:]
    alignment.query_sequence = mut_seq


def open_vcf(filename, mode):
    if filename[-3:] == '.gz':
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def seq2bytes(seq: str):
    """
    :param seq: DNA sequence
    :return: DNA sequence encoded as 2bits per BASE.
    Incomplete last byte is zero-filled.
    """
    seq_bytes = b''
    byte = 0
    for i in range(len(seq)):
        # calculate padding
        lshift = 6 - (i % 4) * 2
        try:
            byte |= BASE2BITS[seq[i]] << lshift
        except KeyError:
            raise ValueError("Illegal DNA base %s" % seq[i])
        
        if (i + 1) % 4 == 0 or len(seq) == i + 1:
            # end of byte or end of sequence
            seq_bytes += bytes([byte])
            byte = 0
    
    return seq_bytes


def bytes2seq(byte_list: bytes, seq_length: int):
    """
    :param byte_list: DNA sequence encoded as 2bits per BASE
    :param seq_length: length of result
    :return: DNA sequence
    """
    if math.ceil(seq_length / 4) > len(byte_list):
        raise ValueError('Not enough bytes for supplied length')
    
    if math.ceil(seq_length / 4) < len(byte_list):
        raise ValueError('Too much bytes for supplied length')
    
    seq = ''
    for i in range(seq_length):
        byte = byte_list[int(i / 4)]
        # calculate padding
        rshift = 6 - (i % 4) * 2
        # extract next 2 bits
        bits = (byte & BASE_BITMASKS[i % 4]) >> rshift
        seq += BITS2BASE[bits]
    
    return seq
