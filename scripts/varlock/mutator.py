import random

import numpy as np

from .common import *
from .diff import Diff
from .vac import Vac


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
    
    def __init__(self, fai_filepath, rnd=random.SystemRandom(), verbose=False):
        """
        :param fai_filepath:
        :param rnd: (secure) random generator
        :param verbose:
        """
        self.fai_list = parse_fai(fai_filepath)
        # create dict for fast access
        self.fai_dict = fai_list2dict(self.fai_list)
        self.rnd = rnd
        self.verbose = verbose
    
    def __is_before_index(self, alignment, index):
        """
        Check if alignment is mapped before index.
        :return: True if alignment end is mapped before index
        """
        # reference_end points to one past the last aligned residue
        ref_end = alignment.reference_end - 1
        alignment_end = pos2index(alignment.reference_name, ref_end, self.fai_dict)
        return alignment_end < index
    
    def __is_after_index(self, alignment, index):
        """
        Check if alignment is mapped after index.
        :return: True if alignment start is mapped after index
        """
        alignment_start = pos2index(alignment.reference_name, alignment.reference_start, self.fai_dict)
        return alignment_start > index
    
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
    
    @staticmethod
    def __read_alignment(bam_file):
        try:
            return next(bam_file)
        except StopIteration:
            return None
    
    def __write_alignment(self, bam_file, alignment):
        """
        :param alignment: pysam.AlignedSegment
        :return:
        """
        bam_file.write(alignment)
        self.alignment_counter += 1
        
        if self.verbose and self.alignment_counter % 10000 == 0:
            print("%d alignments processed" % self.alignment_counter)
    
    @staticmethod
    def ref_pos2seq_pos(alignment, ref_pos):
        """
        Retrieve base position in sequence string at refence position.
        It is assumed that alignment and reference position are on the same reference.
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
    
    @classmethod
    def get_base_pileup(cls, snv_alignments):
        pileup_col = []
        for snv_alignment in snv_alignments:
            if snv_alignment.pos is not None:
                # alignment is mapped at snv position
                snv_base = cls.get_base(snv_alignment.alignment, snv_alignment.pos)
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
    
    def __init_counters(self):
        self.alignment_counter = 0  # counts written alignments
        self.diff_counter = 0
        self.mut_counter = 0
        self.variant_counter = 0  # counts read snvs
        
        self.unmapped_counter = 0
        self.overlapping_counter = 0
        self.max_coverage = 0
    
    def __read_next_diff(self, diff_file):
        try:
            index, mut_tuple = Diff.read_next(diff_file)
        except EOFError:
            return None
        
        self.variant_counter += 1
        ref_name, ref_pos = index2pos(index, self.fai_list)
        return DiffRecord(
            index=index,
            ref_name=ref_name,
            ref_pos=ref_pos,
            mut_map=dict(zip(mut_tuple, BASES))
        )
    
    def __read_next_vac(self, vac_file):
        try:
            index, ac_tuple = Vac.read_next(vac_file)
        except EOFError:
            return None
        
        self.variant_counter += 1
        ref_name, ref_pos = index2pos(index, self.fai_list)
        return VacRecord(
            index=index,
            ref_name=ref_name,
            ref_pos=ref_pos,
            ac=ac_tuple
        )
    
    def unmutate(
            self,
            in_bam_file,
            in_diff_file,
            ref_name,
            start_pos,
            end_pos,
            out_bam_file
    ):
        """
        :param in_bam_file: pysam.AlignmentFile
        :param in_diff_file: binary file
        :param ref_name:
        :param start_pos:
        :param end_pos:
        :param out_bam_file: pysam.AlignmentFile
        :return:
        """
        self.__init_counters()
        
        if start_pos >= end_pos:
            raise ValueError("End position must be greater than start position.")
        
        if not strip_chr(ref_name) in self.fai_dict:
            raise ValueError("Unknown reference name %s" % ref_name)
        
        snv_alignments = []
        bam_region = in_bam_file.fetch(ref_name, start_pos, end_pos)
        alignment = self.__read_alignment(bam_region)
        
        diff_start_index = pos2index(ref_name, start_pos, self.fai_dict)
        diff_end_index = pos2index(ref_name, end_pos, self.fai_dict)
        Diff.seek_index(in_diff_file, diff_start_index, diff_end_index)
        diff = self.__read_next_diff(in_diff_file)
        
        if self.verbose:
            print("Diff start position %d" % in_diff_file.tell())
        
        while True:
            if diff is None or alignment is None:
                # finish
                if self.verbose:
                    if diff is None:
                        print("EOF DIFF")
                    else:
                        print("EOF BAM")
                
                # last mutation
                if diff is not None:
                    # noinspection PyTypeChecker
                    self.__unmutate_overlap(snv_alignments, diff.mut_map)
                
                for snv_alignment in snv_alignments:
                    self.__write_alignment(out_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF VAC case only)
                while alignment is not None:
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = self.__read_alignment(bam_region)
                
                break
            
            elif alignment.is_unmapped:
                self.unmapped_counter += 1
                self.__write_alignment(out_bam_file, alignment)
                alignment = self.__read_alignment(bam_region)
            
            elif self.__is_before_index(alignment, diff.index):
                self.__write_alignment(out_bam_file, alignment)
                alignment = self.__read_alignment(bam_region)
            
            elif self.__is_after_index(alignment, diff.index):
                self.__unmutate_overlap(snv_alignments, diff.mut_map)
                # done with this vac, read next
                diff = self.__read_next_diff(in_diff_file)
                if diff is not None:
                    # could be end of DIFF file
                    self.__write_before_index(out_bam_file, snv_alignments, diff.index)
                    self.__set_seq_positions(snv_alignments, diff.ref_pos)
            
            else:  # alignment is overlapping vac
                # alignment dont have to be mapped to the vac (indel)
                self.overlapping_counter += 1
                # noinspection PyAttributeOutsideInit
                self.max_coverage = max(len(snv_alignments), self.max_coverage)
                # find sequence position of vac
                seq_pos = self.ref_pos2seq_pos(alignment, diff.ref_pos)
                snv_alignments.append(SnvAlignment(alignment, seq_pos))
                
                alignment = self.__read_alignment(bam_region)
        
        if self.verbose:
            self.__print_counters()
    
    def __unmutate_overlap(self, snv_alignments, mut_map):
        for snv_alignment in snv_alignments:
            if snv_alignment.pos is not None:
                # alignment has vac to mutate
                self.__mutate_alignment(snv_alignment.alignment, snv_alignment.pos, mut_map)
                # done with current vac
                snv_alignment.pos = None
    
    def mutate(self, in_vac_file, in_bam_file, out_bam_file, out_diff_file):
        """
        Mutate BAM file SNVs.
        :param in_vac_file: binary file
        input variant allele count file
        :param in_bam_file: pysam.AlignmentFile
        input bam file
        :param out_bam_file: pysam.AlignmentFile
        output bam file
        :param out_diff_file: binary file
        output diff file
        """
        self.__init_counters()
        
        snv_alignments = []
        alignment = self.__read_alignment(in_bam_file)
        vac = self.__read_next_vac(in_vac_file)
        
        while True:
            if vac is None or alignment is None:
                # finish
                if self.verbose:
                    if vac is None:
                        print("EOF VAC")
                    else:
                        print("EOF BAM")
                
                # last mutation
                if vac is not None:
                    # noinspection PyTypeChecker
                    self.__mutate_overlap(out_diff_file, snv_alignments, vac)
                
                for snv_alignment in snv_alignments:
                    self.__write_alignment(out_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF VAC case only)
                while alignment is not None:
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = self.__read_alignment(in_bam_file)
                
                break
            
            elif alignment.is_unmapped:
                self.unmapped_counter += 1
                self.__write_alignment(out_bam_file, alignment)
                alignment = self.__read_alignment(in_bam_file)
            
            elif self.__is_before_index(alignment, vac.index):
                self.__write_alignment(out_bam_file, alignment)
                alignment = self.__read_alignment(in_bam_file)
            
            elif self.__is_after_index(alignment, vac.index):
                self.__mutate_overlap(out_diff_file, snv_alignments, vac)
                # done with this vac, read next
                vac = self.__read_next_vac(in_vac_file)
                if vac is not None:
                    # could be end of VAC file
                    self.__write_before_index(out_bam_file, snv_alignments, vac.index)
                    self.__set_seq_positions(snv_alignments, vac.ref_pos)
            
            else:  # alignment is overlapping vac
                # alignment dont have to be mapped to the vac (indel)
                # find sequence position of vac
                seq_pos = self.ref_pos2seq_pos(alignment, vac.ref_pos)
                snv_alignments.append(SnvAlignment(alignment, seq_pos))
                
                alignment = self.__read_alignment(in_bam_file)
                
                self.overlapping_counter += 1
                # noinspection PyAttributeOutsideInit
                self.max_coverage = max(len(snv_alignments), self.max_coverage)
        
        if self.verbose:
            self.__print_counters()
    
    def __print_counters(self):
        print("written alignments %d" % self.alignment_counter)
        print("unmapped alignments %d" % self.unmapped_counter)
        print("overlapping alignments %d" % self.overlapping_counter)
        print("max coverage %d" % self.max_coverage)
        print("read snvs %d" % self.variant_counter)
        print("mutations %d" % self.mut_counter)
        print("diffs (mutation mappings) %d" % self.diff_counter)
    
    def __mutate_overlap(self, out_diff_file, snv_alignments, vac):
        """
        Mutate alignments at current SNV position by allele frequencies of the SNV
        :param snv_alignments: alignments together with vac positions
        :param vac: SNV
        :return: True if NS mutation has occured
        """
        is_mutated = False
        # get pileup of bases from alignments at vac mapping position
        base_pileup = self.get_base_pileup(snv_alignments)
        alt_ac = count_bases(base_pileup)
        
        mut_map = self.create_mut_map(alt_ac=alt_ac, ref_ac=vac.ac, rnd=self.rnd)
        
        for snv_alignment in snv_alignments:
            if snv_alignment.pos is not None:
                # alignment has vac to mutate
                is_mutated |= self.__mutate_alignment(snv_alignment.alignment, snv_alignment.pos, mut_map)
                # done with current vac
                snv_alignment.pos = None
        
        if is_mutated:
            # at least one alignment has been mutated
            self.__write_diff(
                out_diff_file=out_diff_file,
                index=vac.index,
                mut_map=mut_map
            )
    
    def __mutate_alignment(self, alignment, seq_pos, mut_map):
        """
        Mutate alignment by mutation map at SNV position.
        :param alignment: pysam.AlignedSegment
        :param seq_pos: position of SNV in aligned sequence
        :param mut_map: mutation map
        :return: True if base has been mutated
        """
        is_mutated = False
        # alignment is mapped at snv position
        snv_base = self.get_base(alignment, seq_pos)
        mut_base = mut_map[snv_base]
        
        if snv_base != mut_base:
            # base has been mutated to another base
            self.mut_counter += 1
            is_mutated = True
            self.set_base(alignment, seq_pos, mut_base)
            self.check_cigar_str(alignment, seq_pos)
        
        return is_mutated
    
    # def __next_rel_pos(self, snv):
    #     """
    #     :param snv:
    #     :return: Relative position from last NS SNV
    #     """
    #     if snv.ref_name != self.last_ref_name:
    #         self.last_ref_name = snv.ref_name
    #         self.last_mut_index = self.fai_dict[snv.ref_name].start
    #
    #     rel_pos = snv.index - self.last_mut_index
    #     self.last_mut_index = snv.index
    #     return rel_pos
    
    def __write_diff(self, out_diff_file, index, mut_map):
        self.diff_counter += 1
        mut_tuple = (mut_map['A'], mut_map['T'], mut_map['G'], mut_map['C'])
        Diff.write_next(out_diff_file, index, mut_tuple)
    
    # @classmethod
    # def record2bytes(cls, index, mut_tuple):
    #     """
    #     :param index: absolute genome position
    #     :param mut_tuple: mapping for A, G, C, T bases
    #     :return: bytes
    #     """
    #     struct.pack()
    #     data_list = list(struct.unpack('<IHHHH', byte_string))
    #
    #
    #     bin_pos = bin(pos)[2:]
    #     bits = bitarray(bin_pos.zfill(cls.DIFF_POS_BIT_SIZE)) + bitarray(cls.MUT_2_INDEX[mut_tuple])
    #
    #     return bits.tobytes()
    #
    # @classmethod
    # def bytes2record(cls, byte_string):
    #     """
    #     :param byte_string:
    #     :return: index, mutation map
    #     """
    #     bits = bitarray()
    #     bits.frombytes(byte_string)
    #     rel_pos = int(bits[:19].to01(), 2)
    #     mut_map = cls.BIN_2_MUT[bits[19:].to01()]
    #     return rel_pos, mut_map
    
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
    
    # def __process_snv_alignments(self, out_bam_file, snv_alignments, new_snv):
    #     """
    #     Update ovelapping alignments with new SNV. Save alignments before new snv to file
    #     :param snv_alignments: SNV alignments overlapping previous SNVs
    #     :param new_snv: Next SNV from VAC file. This SNV must not be present in current snv_alignments yet.
    #     :return: list of remaining alignments
    #     """
    #     if new_snv is None:
    #         # could be end of VAC file
    #         return snv_alignments
    #
    #     new_snv_alignments = snv_alignments[:]
    #     for snv_alignment in snv_alignments:
    #         # new_snv alignment is either before or overlapping next new_snv
    #         if self.__is_before_index(snv_alignment.alignment, new_snv.index):
    #             self.__write_alignment(out_bam_file, snv_alignment.alignment)
    #             # remove written alignment
    #             new_snv_alignments.remove(snv_alignment)
    #         else:  # is overlapping new_snv
    #             # find sequence position of new_snv
    #             snv_pos = self.ref_pos2seq_pos(snv_alignment.alignment, new_snv.ref_pos)
    #             if snv_pos is not None:
    #                 # alignment is mapped at another new_snv position
    #                 snv_alignment.snv_pos = snv_pos
    #
    #     return new_snv_alignments
    
    def __write_before_index(self, out_bam_file, snv_alignments, index):
        """
        Write snv alignments before index to file.
        :param out_bam_file:
        :param snv_alignments:
        :param index:
        :return:
        """
        tmp_snv_alignments = snv_alignments[:]
        for snv_alignment in tmp_snv_alignments:
            # new_snv alignment is either before or overlapping next new_snv
            if self.__is_before_index(snv_alignment.alignment, index):
                self.__write_alignment(out_bam_file, snv_alignment.alignment)
                # remove written alignment
                snv_alignments.remove(snv_alignment)
    
    def __set_seq_positions(self, snv_alignments, ref_pos):
        """
        Set SNV position derived from reference position for each alignment.
        :param snv_alignments:
        :param ref_pos:
        :return:
        """
        for snv_alignment in snv_alignments:
            snv_pos = self.ref_pos2seq_pos(snv_alignment.alignment, ref_pos)
            if snv_pos is not None:
                # alignment is mapped at another new_snv position
                snv_alignment.pos = snv_pos
