import binascii
import gzip
import hashlib
import json
import math
import os
import typing

import numpy as np
import pysam

from varlock_src import po
from varlock_src.random import VeryRandom
from varlock_src.variant import AlignedVariant

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


def snv_mut_map(alt_freqs: list, ref_freqs: list, rnd: VeryRandom):
    """
    :param alt_freqs: list of alternative A,T,G,C DNA bases frequencies
    :param ref_freqs: list of reference A,T,G,C DNA bases frequencies
    :param rnd: random number generator
    :return: Map as bijection between bases and premuted bases.
    """
    assert len(alt_freqs) == len(ref_freqs) == len(BASES)
    # add random value to distinguish tied values
    alt_freqs = [ac + rnd.random() for ac in alt_freqs]
    mut_map = _mut_map(
        seqs=list(BASES),
        alt_freqs=alt_freqs,
        ref_freqs=list(ref_freqs),
        rnd=rnd
    )
    mut_map[UNKNOWN_BASE] = UNKNOWN_BASE
    return mut_map


def indel_mut_map(alt_freq_map: dict, ref_freq_map: dict, rnd: VeryRandom):
    """
    Create indel mutation mapping of reference sequences.
    :param alt_freq_map: {seq:count, ...}
    :param ref_freq_map: {seq:count, ...}
    :param rnd: random number generator
    :return: Map as bijection between indels and permuted indels.
    """
    seqs = [None] * len(ref_freq_map)
    alt_freqs = [0] * len(ref_freq_map)
    ref_freqs = [0] * len(ref_freq_map)
    i = 0
    
    # dict needs to be sorted so that order of its items is deterministic
    for seq, freq in sorted(ref_freq_map.items(), key=lambda item: item[0]):
        seqs[i] = seq
        alt_freqs[i] = alt_freq_map.get(seq, 0) + rnd.random()
        ref_freqs[i] = freq
        i += 1
    
    return _mut_map(seqs, alt_freqs, ref_freqs, rnd)


# TODO rename alt_freqs to private_freqs, ref_freqs to public_freqs
def _mut_map(seqs: list, alt_freqs: list, ref_freqs: list, rnd: VeryRandom):
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
        # draw ref allele with multinomial probability
        ref_indel_id = rnd.multirand_index(ref_freqs)
        # draw most abundant alt allele
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


def checksum(filepath: str, as_bytes=False):
    """
    :param filepath:
    :param as_bytes: when true, function returns bytes instead of hex string
    :return: checksum in hex string or byte format
    """
    hasher = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hasher.update(chunk)
    
    checksum_bytes = hasher.digest()
    if as_bytes:
        return checksum_bytes
    else:
        return bytes2hex(checksum_bytes)


def ref_pos2seq_pos(alignment: pysam.AlignedSegment, ref_pos: int):
    """
    Retrieve base position in sequence string at refence position.
    Alignment and ref_pos are assumed to be of the same reference.
    :param alignment: pysam.AlignedSegment
    :param ref_pos: reference position of base
    :return: AlignedSegment.query_sequence position matched to ref_pos.
    None is returned if matching position is not found.
    """
    # TODO optimalize: (try matches_only=True)
    # TODO optimalize: case when alignment is full matched based on CIGAR (e.g. 30M)
    
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
    for variant in variants:  # type: AlignedVariant
        if variant.is_present():
            # alignment is mapped at snv position
            pileup_col.append(variant.seq)
            
            # if len(variant.seq) == 0:
            #     print(variant.alignment)
            #     print(variant.alignment.query_name)
            #     print(variant.alignment.query_sequence)
            #     print(variant._pos)
            #     print(variant._end_pos)
            #     print(variant._is_snv)
            #     exit(0)
    
    return pileup_col


# def max_match_len(sequence: str, pos: int, words: list):
#     """
#     :param sequence:
#     :param pos: position of matching
#     :param words: words to match
#     :return: Length of longest matching word in a string starting
#     at specified position.
#     When at least one word exceeds sequence length, zero is returned.
#     """
#     max_length = 0
#     for word in words:
#         end_pos = pos + len(word)
#         if end_pos > len(sequence):
#             max_length = 0
#             break
#
#         if word == sequence[pos:end_pos]:
#             # matching
#             max_length = max(len(word), max_length)
#
#     return max_length


def max_match_cigar(alignment: pysam.AlignedSegment, pos: int, ref_pos: int, words: list, reference: str):
    """
    Find longest matching sequence in the alignment with respect to the CIGAR string.
    :param alignment: aligned read
    :param pos: starting position in the read
    :param ref_pos: reference position of the starting position
    :param words: list of sequences to match
    :param reference: reference sequence
    :return: the length of longest matching word
    """
    max_length = 0
    
    longest = max(len(word) for word in words)
    
    # skip those that span more than the sequence TODO: upgrade
    if alignment.reference_end < ref_pos + longest:
        return 0
    
    # create how the sequence should look like: (this can be faster, but INDELS are uncommon, so maybe it does not matter?)
    full_cigar = ''
    ref_seq = ''
    writing = False
    for current_seq_pos, current_ref_pos in alignment.get_aligned_pairs():
        if not writing and current_seq_pos == pos:
            # position of the indel is reached
            writing = True
        if writing:
            if current_seq_pos is None:  # deletion
                full_cigar += 'D'
                ref_seq += 'N'
                break
            elif current_ref_pos is None:  # insertion
                full_cigar += 'I'
                ref_seq += alignment.query_sequence[current_seq_pos]
            else:  # match
                full_cigar += 'M'
                ref_seq += alignment.query_sequence[current_seq_pos]
            if len(full_cigar) == longest:
                break
    
    # check reference sequence
    if full_cigar[:len(reference)] == 'M' * len(reference) and ref_seq[:len(reference)] == reference:
        max_length = max(len(reference), max_length)
    
    # check alternative sequences
    for word in words:
        if word != reference:
            # check if this word is available in alignment (insert is inside or it ends with the deletion)
            if ref_seq[:len(word)] == word \
                    and ('I' in full_cigar[:len(word)] or 'D' == full_cigar[len(word):len(word) + 1]):
                max_length = max(len(word), max_length)
    
    return max_length


def is_placed_alignment(alignment: pysam.AlignedSegment):
    """
    :param alignment:
    :return: True if alignment is placed, it still can be unmapped.
    Unplaced alignment has no reference_name and reference_start.
    """
    return alignment.reference_id != -1


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


def create_aligned_variant(
        alignment: pysam.AlignedSegment,
        vac: typing.Union[po.VariantPosition, po.VariantDiff],
        vac_occurrence: bool,
        is_mutated: bool = False
):
    """
    Factory that creates AlignedVariant from pysam alignment and VAC record.
    Alignment and VAC are assumed to be of the same reference.
    :param alignment: aligned read
    :param vac: variant position of diff record if vac_occurrence is False
    :param vac_occurrence: if we compare to Snv/IndelOccurrence (encrypt) or to Snv/IndelDiff (decrpyt)
    :param is_mutated:
    """
    assert alignment.reference_name == vac.ref_name
    
    SnpCompare = po.SnvOccurrence
    IndelCompare = po.IndelOccurrence
    if not vac_occurrence:
        SnpCompare = po.SnvDiff
        IndelCompare = po.IndelDiff
    
    if alignment.is_unmapped or (vac.ref_pos >= alignment.reference_end or vac.ref_pos < alignment.reference_start):
        # variant is unmapped or either after or before the alignment
        variant = AlignedVariant(alignment, is_mutated=is_mutated)
    else:
        pos = ref_pos2seq_pos(alignment, vac.ref_pos)
        if pos is None:
            # message = "WARNING: reference position %d on alignment with range <%d,%d> not found (possible deletion in bam read)" \
            #           % (vac.ref_pos, alignment.reference_start, alignment.reference_end - 1)
            # print(message)
            variant = AlignedVariant(alignment, is_mutated=is_mutated)
        elif isinstance(vac, SnpCompare):
            variant = AlignedVariant(alignment, pos, is_mutated=is_mutated)
        elif isinstance(vac, IndelCompare):
            words = vac.seqs if vac_occurrence else list(vac.mut_map.values())
            end_pos = pos + max_match_cigar(alignment, pos, vac.ref_pos, words, vac.ref_seq)
            if end_pos > pos:
                # indel was found
                if vac_occurrence:
                    assert len(vac.seqs)
                
                if alignment.query_name == 'ERR015528.28046840' and alignment.reference_start == 356649:
                    print('YYY')
                    print('vac_pos: %s' % vac.ref_pos)
                    print('vac_seqs: %s' % vac.seqs)
                    print('ref_seq: %s' % vac.ref_seq)
                    # print('infer_query_length %d' % alignment.infer_query_length())
                    # print('query length %d' % len(alignment.query_sequence))
                    # print(alignment.cigarstring)
                    # print(alignment.query_sequence)
                    print('YYY')
                    # exit(0)
                
                variant = AlignedVariant(alignment, pos, end_pos, vac.ref_seq, is_mutated)
            else:
                # match not found or at least one variant exceeds alignment end
                variant = AlignedVariant(alignment, is_mutated=is_mutated)
        else:
            raise ValueError("%s is not %s/%s record instance" % (
                type(vac).__name__, type(SnpCompare).__name__, type(IndelCompare).__name__,))
    
    return variant
