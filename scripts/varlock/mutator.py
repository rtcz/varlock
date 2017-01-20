import random
import struct

import numpy as np
import pysam
from bitarray import bitarray

from .common import *


class Mutator:
    CIGAR_MATCH = 0  # M
    CIGAR_INS = 1  # I
    CIGAR_DEL = 2  # D
    CIGAR_REF_SKIP = 3  # N
    CIGAR_SOFT_CLIP = 4  # S
    CIGAR_HARD_CLIP = 5  # H
    CIGAR_PAD = 6  # P
    CIGAR_EQUAL = 7  # E
    CIGAR_DIFF = 8  # X
    NM_TAG = 9
    
    def __init__(self, fai_list, rnd=random.SystemRandom(), verbose=False):
        """
        :param fai_list: parsed FAI as list
        :param rnd: (secure) random generator
        :param verbose:
        """
        self.fai_list = fai_list
        # create dict for fast access
        self.fai_dict = fai_list2dict(self.fai_list)
        self.rnd = rnd
        self.verbose = verbose
    
    def __is_before_snv(self, alignment, snv):
        """
        Check if alignment is mapped before SNV.
        :return: True if alignment end is mapped before SNV
        """
        # reference_end points to one past the last aligned residue
        ref_end = alignment.reference_end - 1
        alignment_end = pos2index(alignment.reference_name, ref_end, self.fai_dict)
        return alignment_end < snv.index
    
    def __is_after_snv(self, alignment, snv):
        """
        Check if alignment is mapped after SNV.
        :return: True if alignment start is mapped after SNV
        """
        alignment_start = pos2index(alignment.reference_name, alignment.reference_start, self.fai_dict)
        return alignment_start > snv.index
    
    @staticmethod
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
            alt_base_id = np.argmax(alt_ac)
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
    
    def __read_alignment(self):
        try:
            return next(iter(self.in_bam_file))
        except StopIteration:
            return None
    
    def __read_snv(self):
        byte_string = self.in_vac_file.read(12)
        if len(byte_string) == 0:
            return None
        self.snv_counter += 1
        data_list = list(struct.unpack('<IHHHH', byte_string))
        ref_name, ref_pos = index2pos(data_list[0], self.fai_list)
        
        return Snv(index=data_list[0], ref_name=ref_name, ref_pos=ref_pos, ac=data_list[1:])
    
    @staticmethod
    def ref_pos2seq_pos(alignment, ref_pos):
        """
        Retrieve base position in sequence string at refence position.
        :param alignment: pysam.AlignedSegment
        :param ref_pos: reference position of base
        :return: 0-based position of base with specified reference position in sequence string
        None if alignment is not mapped at ref_pos (deletion)
        """
        seq_pos = None
        for current_seq_pos, current_ref_pos in alignment.get_aligned_pairs(matches_only=False, with_seq=False):
            # search for base in snv position
            if current_ref_pos == ref_pos:
                seq_pos = current_seq_pos
                break
        
        return seq_pos
    
    def __write_alignment(self, alignment):
        """
        :param alignment: pysam.AlignedSegment
        :return:
        """
        self.out_bam_file.write(alignment)
        self.alignment_counter += 1
        
        if self.verbose and self.alignment_counter % 10000 == 0:
            print("%d alignments processed" % self.alignment_counter)
    
    @classmethod
    def get_base_pileup(cls, snv_alignments):
        pileup_col = []
        for snv_alignment in snv_alignments:
            if snv_alignment.snv_pos is not None:
                # alignment is mapped at snv position
                snv_base = cls.get_base(snv_alignment.alignment, snv_alignment.snv_pos)
                pileup_col.append(snv_base)
        
        return pileup_col
    
    @staticmethod
    def get_base(alignment, pos):
        """
        :param alignment: pysam.AlignedSegment
        :param pos: position in sequence
        :return: base at pos
        """
        return alignment.query_sequence[pos]
    
    @staticmethod
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
    
    def __finish(self, snv, alignment, snv_alignments):
        """
        Write remaining alignments to BAM file.
        :param snv:
        :param alignment:
        :param snv_alignments:
        :return:
        """
        if self.verbose:
            if snv is None:
                print("EOF VAC")
            else:
                print("EOF BAM")
        
        # last mutation
        if snv is not None:
            self.__mutate_snv_alignments(snv_alignments, snv)
        
        for snv_alignment in snv_alignments:
            self.__write_alignment(snv_alignment.alignment)
        
        # write remaining alignments (in EOF VAC case only)
        while alignment is not None:
            self.__write_alignment(alignment)
            alignment = self.__read_alignment()
    
    def __init_counters(self):
        self.alignment_counter = 0  # counts written alignments
        self.diff_counter = 0
        self.mut_counter = 0
        self.snv_counter = 0  # counts read snvs
        
        self.unmapped_counter = 0
        self.overlapping_counter = 0
        self.max_snv_alignments = 0
        
        # counts relative mutation positions that exceed SNV_POS_BIT_SIZE
        self.sparse_mut_counter = 0
    
    def __init_files(self):
        self.in_bam_file = None
        self.in_vac_file = None
        self.out_bam_file = None
        self.out_diff_file = None
    
    def __init_last(self):
        # relative index to last mutation position
        self.last_mut_index = 0
        
        # reference name in last snv
        self.last_ref_name = None
    
    def mutate(self, in_vac_filename, in_bam_filename, out_bam_filename, out_diff_filename):
        """
        Mutate BAM file SNVs.
        :param in_vac_filename: input variant allele count filename
        :param in_bam_filename: input bam filename
        :param out_bam_filename: output bam filename
        :param out_diff_filename: output diff filename
        """
        self.__init_counters()
        self.__init_files()
        self.__init_last()
        
        with pysam.AlignmentFile(in_bam_filename, "rb") as self.in_bam_file, \
                open(in_vac_filename, "rb") as self.in_vac_file, \
                pysam.AlignmentFile(out_bam_filename, "wb", template=self.in_bam_file) as self.out_bam_file, \
                open(out_diff_filename, "wb") as self.out_diff_file:
            
            snv_alignments = []
            alignment = self.__read_alignment()
            snv = self.__read_snv()
            
            while True:
                if snv is None or alignment is None:
                    self.__finish(snv, alignment, snv_alignments)
                    break
                
                elif alignment.is_unmapped:
                    self.__write_alignment(alignment)
                    self.unmapped_counter += 1
                
                elif self.__is_before_snv(alignment, snv):
                    # print("alignment is before snv")
                    # time.sleep(0.1)
                    self.__write_alignment(alignment)
                    alignment = self.__read_alignment()
                
                elif self.__is_after_snv(alignment, snv):
                    self.__mutate_snv_alignments(snv_alignments, snv)
                    # done with this snv, read next
                    snv = self.__read_snv()
                    # process snv_alignments with new snv
                    snv_alignments = self.__process_snv_alignments(snv_alignments, snv)
                
                else:  # alignment is overlapping snv
                    # alignment dont have to be mapped to the snv (indel)
                    # find sequence position of snv
                    snv_pos = self.ref_pos2seq_pos(alignment, snv.ref_pos)
                    snv_alignments.append(SnvAlignment(alignment, snv_pos))
                    alignment = self.__read_alignment()
                    self.overlapping_counter += 1
                    # noinspection PyAttributeOutsideInit
                    self.max_snv_alignments = max(len(snv_alignments), self.max_snv_alignments)
            
            if self.verbose:
                print("written alignments %d" % self.alignment_counter)
                print("unmapped alignments %d" % self.unmapped_counter)
                print("overlapping alignments %d" % self.overlapping_counter)
                print("max snv alignments %d" % self.max_snv_alignments)
                print("read snvs %d" % self.snv_counter)
                print("mutations %d" % self.mut_counter)
                print("diffs (mutation mappings) %d" % self.diff_counter)
                print("sparse mutations %d" % self.sparse_mut_counter)
    
    def __mutate_snv_alignments(self, snv_alignments, snv):
        """
        Mutate alignments at current SNV position by allele frequencies of the SNV
        :param snv_alignments: alignments together with snv positions
        :param snv: SNV
        :return: True if NS mutation has occured
        """
        is_mutated = False
        # get pileup of bases from alignments at snv mapping position
        base_pileup = self.get_base_pileup(snv_alignments)
        alt_ac = count_bases(base_pileup)
        
        mut_map = self.create_mut_map(alt_ac=alt_ac, ref_ac=snv.ac, rnd=self.rnd)
        
        for snv_alignment in snv_alignments:
            if snv_alignment.snv_pos is not None:
                # alignment has snv to mutate
                is_mutated |= self.__mutate_alignment(snv_alignment.alignment, snv_alignment.snv_pos, mut_map)
                # done with current snv
                snv_alignment.snv_pos = None
        
        if is_mutated:
            # at least one alignment has been mutated
            self.__write_diff(
                snv=snv,
                mut_map=mut_map
            )
    
    def __mutate_alignment(self, alignment, snv_pos, mut_map):
        """
        Mutate alignment by mutation map at SNV position.
        :param alignment: pysam.AlignedSegment
        :param snv_pos: position of SNV in aligned sequence
        :param mut_map: mutation map
        :return: True if base has been mutated
        """
        is_mutated = False
        # alignment is mapped at snv position
        snv_base = self.get_base(alignment, snv_pos)
        mut_base = mut_map[snv_base]
        
        if snv_base != mut_base:
            # base has been mutated to another base
            self.mut_counter += 1
            is_mutated = True
            self.set_base(alignment, snv_pos, mut_base)
            self.check_cigar_str(alignment, snv_pos)
        
        return is_mutated
    
    def __next_rel_pos(self, snv):
        """
        :param snv:
        :return: Relative position to last NS SNV
        """
        if snv.ref_name != self.last_ref_name:
            self.last_ref_name = snv.ref_name
            self.last_mut_index = self.fai_dict[snv.ref_name].start
        
        rel_pos = snv.index - self.last_mut_index
        self.last_mut_index = snv.index
        return rel_pos
    
    BIN_2_MUT = {
        '00000': ('A', 'T', 'C', 'G'),
        '00001': ('A', 'T', 'G', 'C'),
        '00010': ('A', 'C', 'T', 'G'),
        '00011': ('A', 'C', 'G', 'T'),
        '00100': ('A', 'G', 'T', 'C'),
        '00101': ('A', 'G', 'C', 'T'),
        '00110': ('T', 'A', 'C', 'G'),
        '00111': ('T', 'A', 'G', 'C'),
        '01000': ('T', 'C', 'A', 'G'),
        '01001': ('T', 'C', 'G', 'A'),
        '01010': ('T', 'G', 'A', 'C'),
        '01011': ('T', 'G', 'C', 'A'),
        '01100': ('C', 'A', 'T', 'G'),
        '01101': ('C', 'A', 'G', 'T'),
        '01110': ('C', 'T', 'A', 'G'),
        '01111': ('C', 'T', 'G', 'A'),
        '10000': ('C', 'G', 'A', 'T'),
        '10001': ('C', 'G', 'T', 'A'),
        '10010': ('G', 'A', 'T', 'C'),
        '10011': ('G', 'A', 'C', 'T'),
        '10100': ('G', 'T', 'A', 'C'),
        '10101': ('G', 'T', 'C', 'A'),
        '10110': ('G', 'C', 'A', 'T'),
        '10111': ('G', 'C', 'T', 'A')
    }
    
    MUT_2_BIN = {
        ('A', 'T', 'C', 'G'): '00000',
        ('A', 'T', 'G', 'C'): '00001',
        ('A', 'C', 'T', 'G'): '00010',
        ('A', 'C', 'G', 'T'): '00011',
        ('A', 'G', 'T', 'C'): '00100',
        ('A', 'G', 'C', 'T'): '00101',
        ('T', 'A', 'C', 'G'): '00110',
        ('T', 'A', 'G', 'C'): '00111',
        ('T', 'C', 'A', 'G'): '01000',
        ('T', 'C', 'G', 'A'): '01001',
        ('T', 'G', 'A', 'C'): '01010',
        ('T', 'G', 'C', 'A'): '01011',
        ('C', 'A', 'T', 'G'): '01100',
        ('C', 'A', 'G', 'T'): '01101',
        ('C', 'T', 'A', 'G'): '01110',
        ('C', 'T', 'G', 'A'): '01111',
        ('C', 'G', 'A', 'T'): '10000',
        ('C', 'G', 'T', 'A'): '10001',
        ('G', 'A', 'T', 'C'): '10010',
        ('G', 'A', 'C', 'T'): '10011',
        ('G', 'T', 'A', 'C'): '10100',
        ('G', 'T', 'C', 'A'): '10101',
        ('G', 'C', 'A', 'T'): '10110',
        ('G', 'C', 'T', 'A'): '10111',
        
    }
    
    # variant bit size
    DIFF_BIT_SIZE = 24
    DIFF_POS_BIT_SIZE = 19
    MAX_DIFF_POS_VALUE = 2 ** DIFF_POS_BIT_SIZE - 1
    
    def __write_diff(self, snv, mut_map):
        rel_pos = self.__next_rel_pos(snv)
        self.diff_counter += 1
        
        while rel_pos > self.MAX_DIFF_POS_VALUE:
            self.sparse_mut_counter += 1
            # write dummy diff to keep maximum relative size
            
            self.out_diff_file.write(self.diff2bytes(self.MAX_DIFF_POS_VALUE, BASES))
            rel_pos -= self.MAX_DIFF_POS_VALUE
        
        mut_tuple = (mut_map['A'], mut_map['T'], mut_map['G'], mut_map['C'])
        self.out_diff_file.write(self.diff2bytes(rel_pos, mut_tuple))
    
    @classmethod
    def diff2bytes(cls, pos, mut_tuple):
        """
        :param pos: relative position
        :param mut_tuple: mapping for A, G, C, T bases
        :return:
        """
        bin_pos = bin(pos)[2:]
        bits = bitarray(bin_pos.zfill(cls.DIFF_POS_BIT_SIZE)) + bitarray(cls.MUT_2_BIN[mut_tuple])
        
        return bits.tobytes()
    
    @classmethod
    def check_cigar_str(cls, alignment, snv_pos):
        """
        Validate cigar operation at SNV position.
        TODO: Update cigar and quality strings.
        :param alignment: pysam.AlignedSegment
        :param snv_pos: position of SNV in alignment sequence
        """
        # start at -1 to work with 0-based position
        seq_pos = -1
        for cig_op_id, cig_op_len in alignment.cigartuples:
            if cig_op_id == cls.CIGAR_MATCH or cig_op_id == cls.CIGAR_INS:
                # match and insert are present in aligned sequence
                tmp_seq_pos = seq_pos + cig_op_len
            elif cig_op_id == cls.CIGAR_DEL:
                # deletion is not present in aligned sequence
                continue
            else:
                raise ValueError("Unsupported CIGAR operation %d" % cig_op_id)
            
            if tmp_seq_pos < snv_pos:
                # curr_seq_pos has not reached snv_seq_pos
                seq_pos = tmp_seq_pos
            else:
                # cigar segment covers snv position, it must be a match
                if cig_op_id != cls.CIGAR_MATCH:
                    # snv is not matched
                    # TODO update mapping quality
                    # TODO update cigar string
                    
                    print("snv_pos %d" % snv_pos)
                    print("seq_pos %d" % seq_pos)
                    print(alignment.cigarstring)
                    print(alignment.get_aligned_pairs(matches_only=False, with_seq=False))
                    
                    raise ValueError("Unsupported CIGAR operation %d" % cig_op_id)
                break
    
    def __process_snv_alignments(self, snv_alignments, new_snv):
        """
        Update ovelapping alignments with new SNV. Save alignments before new snv to file
        :param snv_alignments: SNV alignments overlapping previous SNVs
        :param new_snv: Next SNV from VAC file. This SNV must not be present in current snv_alignments yet.
        :return: list of remaining alignments
        """
        if new_snv is None:
            # could be end of VAC file
            return snv_alignments
        
        new_snv_alignments = snv_alignments[:]
        for snv_alignment in snv_alignments:
            # new_snv alignment is either before or overlapping next new_snv
            if self.__is_before_snv(snv_alignment.alignment, new_snv):
                self.__write_alignment(snv_alignment.alignment)
                # remove written alignment
                new_snv_alignments.remove(snv_alignment)
            else:  # is overlapping new_snv
                # find sequence position of new_snv
                snv_pos = self.ref_pos2seq_pos(snv_alignment.alignment, new_snv.ref_pos)
                if snv_pos is not None:
                    # alignment is mapped at another new_snv position
                    snv_alignment.snv_pos = snv_pos
        
        return new_snv_alignments
