import binascii
import gzip
import hashlib
import json
import math
import os

import numpy as np
import pysam

from varlock_src.random import VeryRandom

BASES = ("A", "T", "G", "C")
UNKNOWN_BASE = "N"

BASE_BITMASKS = [0b11000000, 0b00110000, 0b00001100, 0b00000011]
# BASE_BITMASKS = [0b00000011, 0b00110000, 0b00001100, 0b11000000]
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


def snv_mut_map(private_freqs: list, public_freqs: list, rnd: VeryRandom):
    """
    :param private_freqs: list of private allele frequencies
    :param public_freqs: list of public allele frequencies
    :param rnd: random number generator
    :return: Map as bijection between bases and premuted bases.
    """
    assert len(private_freqs) == len(public_freqs) == len(BASES)
    
    mut_map = _mut_map(
        seqs=list(BASES),
        private_freqs=[freq + rnd.random() for freq in private_freqs],
        public_freqs=public_freqs.copy(),
        rnd=rnd
    )
    return mut_map


# TODO use parameters: alleles, alt_freqs, ref_freqs
def indel_mut_map(private_freq_map: dict, public_freq_map: dict, rnd: VeryRandom):
    """
    Create indel mutation mapping of reference sequences.
    :param private_freq_map: {seq:count, ...}
    :param public_freq_map: {seq:count, ...}
    :param rnd: random number generator
    :return: Map as bijection between indels and permuted indels.
    """
    seqs = [None] * len(public_freq_map)
    alt_freqs = [0] * len(public_freq_map)
    ref_freqs = [0] * len(public_freq_map)
    i = 0
    
    # dict needs to be sorted so that order of its items is deterministic
    for seq, freq in sorted(public_freq_map.items(), key=lambda item: item[0]):
        seqs[i] = seq
        alt_freqs[i] = private_freq_map.get(seq, 0) + rnd.random()
        ref_freqs[i] = freq
        i += 1
    
    return _mut_map(seqs, alt_freqs, ref_freqs, rnd)


# TODO rename alt_freqs to private_freqs, ref_freqs to public_freqs
def _mut_map(seqs: list, private_freqs: list, public_freqs: list, rnd: VeryRandom):
    """
    Function tampers with parameters to avoid list copying.
    :param seqs: sequences to mutate
    :param private_freqs: private allele frequencies
    :param public_freqs: public allele frequencies
    :return: Map as bijection between seqs and permuted seqs.
    """
    assert len(seqs) == len(private_freqs) == len(public_freqs)
    public_seqs = seqs
    private_seqs = seqs[:]
    
    mut_map = {}
    # map bases but skip last unmapped base
    
    for i in range(len(public_freqs) - 1):
        # draw ref allele with multinomial probability
        public_allele_id = rnd.multirand_index(public_freqs)
        
        # draw most abundant alt allele
        # TODO why not use multi random here?
        private_allele_id = np.argmax(private_freqs)  # type: int
        # add mapping
        mut_map[private_seqs[private_allele_id]] = public_seqs[public_allele_id]
        
        # delete processed items
        del public_seqs[public_allele_id]
        del public_freqs[public_allele_id]
        del private_seqs[private_allele_id]
        del private_freqs[private_allele_id]
    
    # last base mapping is obvious
    mut_map[private_seqs[0]] = public_seqs[0]
    
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
    :param pileup: list of DNA bases: A, T, G, C
    :return: list of DNA base frequencies [A_freq, T_freq, G_freq, C_freq]
    """
    freqs = [0] * 4
    for base in pileup:
        try:
            freqs[BASES.index(base)] += 1
        except ValueError:
            raise ValueError("Illegal DNA base: %s" % base)
    
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


# TODO rewrite
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


# TODO rewrite
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


def byte2base_perm(perm_byte: int) -> list:
    """
    Convert single byte to permuted list of A,T,G,C bases in respective order.
    Each 2 bits represent one base.
    :return: permutation list
    """
    perm_list = [None] * 4
    for i in range(4):
        shift = 2 * i
        base_id = (perm_byte & (0b11 << shift)) >> shift
        perm_list[i] = BASES[base_id]
    
    return perm_list


def base_perm2byte(base_perm: list) -> int:
    """
    Convert permuted list of A,T,G,C bases to single byte in respective order.
    Each 2 bits represent one base.
    :return: permutation byte
    """
    assert len(base_perm) == 4
    mapping_byte = 0
    for i in range(4):
        shift = 2 * i
        mapping_byte |= BASES.index(base_perm[i]) << shift
    
    return mapping_byte


def dict2bytes(value: dict):
    return json.dumps(value, separators=(',', ':')).encode()


def bytes2dict(value: bytes):
    return json.loads(value.decode())
