import numpy as np

import binascii
import gzip
import hashlib
import math
import os
import json
import varlock.po as po

import pysam

BASES = ("A", "T", "G", "C")
UNKNOWN_BASE = "N"

BASE_BITMASKS = [0b11000000, 0b00110000, 0b00001100, 0b00000011]
BASE2BITS = {'A': 0b00, 'T': 0b01, 'G': 0b10, 'C': 0b11}
BITS2BASE = {0b00: 'A', 0b01: 'T', 0b10: 'G', 0b11: 'C'}


def stream_cipher(seq: str, key: bytes):
    """
    :param seq: sequence of DNA base letters
    :param key: key stream
    :return: cipher stream in form of DNA base letters
    DNA base encryption using stream cipher.
    Input sequence of DNA bases is encrypted as another sequence of the same length.
    Each 2 bits of the key stream are used to encrypt one base from the sequence.
    If there is not enough bits of the key stream to encrypt whole sequence, key stream is repeats.
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


def snv_mut_map(alt_freqs: list, ref_freqs: list, rnd):
    """
    :param alt_freqs: list of alternative A,T,G,C DNA bases frequencies
    :param ref_freqs: list of reference A,T,G,C DNA bases frequencies
    :param rnd: random number generator
    :return: Map as bijection between bases and premuted bases.
    """
    assert len(alt_freqs) == len(ref_freqs) == len(BASES)
    # add random value to distinguish tied values
    alt_freqs = [ac + rnd.random() for ac in alt_freqs]
    mut_map = _mut_map(list(BASES), alt_freqs, list(ref_freqs), rnd)
    mut_map[UNKNOWN_BASE] = UNKNOWN_BASE
    return mut_map


def indel_mut_map(alt_freq_map: dict, ref_freq_map: dict, rnd):
    """
    Create indel mutation mapping from intersection
    of alternative indels with reference indels.
    :param alt_freq_map: {seq:count, ...}
    :param ref_freq_map: {seq:count, ...}
    :param rnd: random number generator
    :return: Map as bijection between indels and permuted indels.
    """
    # find intersection of indels
    seqs = list(set(alt_freq_map) & set(ref_freq_map))
    
    alt_freqs = [0] * len(seqs)
    ref_freqs = [0] * len(seqs)
    for i in range(len(seqs)):
        # add random value to distinguish tied values
        alt_freqs[i] = alt_freq_map[seqs[i]] + rnd.random()
        ref_freqs[i] = ref_freq_map[seqs[i]]
    
    return _mut_map(seqs, alt_freqs, ref_freqs, rnd)


def _mut_map(seqs: list, alt_freqs: list, ref_freqs: list, rnd):
    """
    Function tampers with parameters to avoid list copying.
    :param seqs: sequences to mutate
    :param alt_freqs: alternative sequences frequencies
    :param ref_freqs: reference sequences frequencies
    :return: Map as bijection between seqs and permuted seqs.
    """
    assert len(seqs) == len(alt_freqs) == len(ref_freqs)
    ref_seqs = seqs
    alt_seqs = seqs[:]
    
    mut_map = {}
    # map bases but skip last unmapped base
    
    for i in range(len(ref_freqs) - 1):
        # draw ref indel with multinomial probability
        ref_indel_id = multi_random(ref_freqs, rnd)
        # draw most abundant indel from alt alleles
        # TODO why not use multi random here?
        alt_indel_id = np.argmax(alt_freqs)  # type: int
        # add mapping
        mut_map[alt_seqs[alt_indel_id]] = ref_seqs[ref_indel_id]
        
        # delete processed items
        del ref_seqs[ref_indel_id]
        del ref_freqs[ref_indel_id]
        del alt_seqs[alt_indel_id]
        del alt_freqs[alt_indel_id]
    
    # last base mapping is obvious
    mut_map[alt_seqs[0]] = ref_seqs[0]
    
    return mut_map


def multi_random(p_dist: list, rnd):
    """
    Draw index from multinomial probability distribution.
    :param p_dist: Probability distribution. When all probabilities are zero, each outcome has equal probability.
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


def freq_map(values: list):
    result = {}
    for value in values:
        if value not in result:
            result[value] = 1
        else:
            result[value] += 1
    
    return result


def base_freqs(pileup: list):
    """
    Count allele count in base pileup column.
    :param pileup: list of DNA bases: A, T, G, C, N
    :return: list of DNA base frequencies [A_freq, T_freq, G_freq, C_freq]
    """
    freqs = [0] * 4
    for base in pileup:
        # skip unknown base
        if base != UNKNOWN_BASE:
            try:
                freqs[BASES.index(base)] += 1
            except KeyError:
                raise ValueError("Illegal DNA base %s" % base)
    
    return freqs


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


def bytes2hex(byte_str: bytes):
    return binascii.hexlify(byte_str).decode()


def hex2bytes(hex_str: str):
    return binascii.unhexlify(hex_str)


def file_size(file_obj):
    old_pos = file_obj.tell()
    file_obj.seek(0, os.SEEK_END)
    size = file_obj.tell()
    file_obj.seek(old_pos, os.SEEK_SET)
    return size


def checksum(filename: str):
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


def variant_seqs(variants: list):
    """
    :param variants: list of po.AlignedVariant
    :return: list of bases at specific position
    """
    pileup_col = []
    for variant in variants:  # type: po.AlignedVariant
        if variant.is_present():
            # alignment is mapped at snv position
            pileup_col.append(variant.seq)
    
    return pileup_col


def max_match_len(string: str, pos: int, words: list):
    """
    :param string:
    :param pos: position of matching
    :param words: words to match
    :return: length of longest matching word in a string at specific position
    """
    max_length = 0
    for word in words:
        if word == string[pos:len(word)]:
            # matching
            max_length = max(len(word), max_length)
    
    return max_length


def diff_aligned_variant(alignment: pysam.AlignedSegment, diff: po.GenomicPosition):
    """
    Factory that creates AlignedVariant from pysam alignment and DIFF record.
    :param alignment:
    :param diff:
    """
    # TODO optimize, skip when ref_pos is greater than alignment reference end ?
    pos = ref_pos2seq_pos(alignment, diff.ref_pos)
    if pos is None:
        # this should not occur
        raise ValueError("variant position not found")
    elif isinstance(diff, po.DiffSnvRecord):
        variant = po.AlignedVariant(alignment, pos)
    elif isinstance(diff, po.DiffIndelRecord):
        end_pos = pos + max_match_len(alignment.query_sequence, pos, diff.mut_map.values())
        if end_pos > pos:
            # there was a match
            variant = po.AlignedVariant(alignment, pos, end_pos)
        else:
            # match not found
            # this sould be rare
            variant = po.AlignedVariant(alignment)
    else:
        raise ValueError("%s is not DIFF record instance" % type(diff).__name__)
    
    return variant


def vac_aligned_variant(alignment: pysam.AlignedSegment, vac: po.GenomicPosition):
    """
    Factory that creates AlignedVariant from pysam alignment and VAC record.
    :param alignment:
    :param vac:
    """
    # TODO optimize, skip when ref_pos is greater than alignment reference end ?
    pos = ref_pos2seq_pos(alignment, vac.ref_pos)
    if pos is None:
        # this should not occur
        raise ValueError("variant position not found")
    elif isinstance(vac, po.VacSnvRecord):
        variant = po.AlignedVariant(alignment, pos)
    elif isinstance(vac, po.VacIndelRecord):
        end_pos = pos + max_match_len(alignment.query_sequence, pos, vac.seqs)
        if end_pos > pos:
            # there was a match
            variant = po.AlignedVariant(alignment, pos, end_pos)
        else:
            # match not found
            # this sould be rare
            variant = po.AlignedVariant(alignment)
    else:
        raise ValueError("%s is not VAC record instance" % type(vac).__name__)
    
    return variant


# def get_seq(alignment: pysam.AlignedSegment, pos: int, ref_seqs: list):
#     """
#     :param alignment:
#     :param pos: position
#     :param ref_seqs: reference sequences
#     :return: longest sequence in alignment at specific position matching reference sequence
#     """
#     max_seq = ''
#     for ref_seq in ref_seqs:
#         if ref_seq == alignment.query_sequence[pos, len(ref_seq)]:
#             # matching
#             if len(ref_seq) > len(max_seq):
#                 # longest
#                 max_seq = max_seq
#
#     return max_seq


# def get_base(alignment: pysam.AlignedSegment, pos: int):
#     """
#     :param alignment:
#     :param pos: position in sequence
#     :return: base at pos
#     """
#     return alignment.query_sequence[pos]


def is_placed_alignment(alignment: pysam.AlignedSegment):
    """
    :param alignment:
    :return: True if alignment is placed, it still can be unmapped.
    Unplaced alignment has no reference_name and reference_start.
    """
    return alignment.reference_id != -1


# def set_base(alignment, pos: int, base: str):
#     """
#     Replace base at SNV position
#     :param alignment: pysam.AlignedSegment
#     :param pos: position in sequence
#     :param base: mutated base letter
#     :return: mutated sequence string
#     """
#     # TODO query_qualities
#     mut_seq = alignment.query_sequence[:pos]
#     mut_seq += base
#     mut_seq += alignment.query_sequence[pos + 1:]
#     alignment.query_sequence = mut_seq


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
    :param byte_list: DNA sequence encoded as 2 bits per BASE
    :param seq_length: length of result sequence
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


def dict2bytes(value: dict):
    return json.dumps(value, separators=(',', ':')).encode()


def bytes2dict(value: bytes):
    return json.loads(value.decode())
