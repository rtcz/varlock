import io
import random

import numpy as np

from .common import *
from .diff import Diff
from .fasta_index import FastaIndex
from .iterator import BamIterator, DiffIterator, VacIterator
from .po import SnvAlignment


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
    
    def __init__(self, bam_file, rnd=random.SystemRandom(), verbose=False):
        """
        :param rnd: (secure) random generator
        :param verbose:
        """
        self.rnd = rnd
        self.verbose = verbose
        self.fai = FastaIndex(bam_file)
        self.bam_file = bam_file
        self._bam_checksum = None
    
    def __is_before_index(self, alignment, index):
        """
        Check if alignment is mapped before index.
        :return: True if alignment end is mapped before index
        """
        # reference_end points to one past the last aligned residue
        ref_end = alignment.reference_end - 1
        alignment_end = self.fai.pos2index(alignment.reference_name, ref_end)
        return alignment_end < index
    
    def __is_after_index(self, alignment, index):
        """
        Check if alignment is mapped after index.
        :return: True if alignment start is mapped after index
        """
        alignment_start = self.fai.pos2index(alignment.reference_name, alignment.reference_start)
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
        self.alignment_counter = 0  # written alignments
        self.diff_counter = 0  # processed diff records
        self.mut_counter = 0  # mutations pre alignment
        self.snv_counter = 0  # read vac records (SNVs)
        
        self.unmapped_counter = 0  # unmapped alignments
        self.overlapping_counter = 0  # overlapping alignments
        self.max_coverage = 0  # maximum alignments overlapping single SNV
    
    def __resolve_range(
            self,
            diff_file,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
    ):
        checksum, diff_start_index, diff_end_index = Diff.read_header(diff_file)
        self.validate_checksum(checksum)
        
        return self.resolve_range(
            diff_start_index,
            diff_end_index,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
        )
    
    def resolve_range(
            self,
            diff_start_index,
            diff_end_index,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
    ):
        """
        Resolve between DIFF range and user specified range.
        :param diff_start_index:
        :param diff_end_index:
        :param start_ref_name:
        :param start_ref_pos:
        :param end_ref_name:
        :param end_ref_pos:
        :return: tuple (start_index, end_index)
        """
        if start_ref_name is not None or end_ref_name is not None:
            # apply user range
            start_index = self.fai.resolve_start_index(start_ref_name, start_ref_pos)
            end_index = self.fai.resolve_end_index(end_ref_name, end_ref_pos)
            
            diff_range = self.fai.index2pos(diff_start_index) + self.fai.index2pos(diff_end_index)
            
            if start_index < diff_start_index:
                args = self.fai.index2pos(start_index) + diff_range
                # noinspection PyStringFormat
                raise ValueError("Start position [%s, %d] must be within DIFF range [%s, %d]:[%s, %d]." % args)
            if end_index > diff_end_index:
                args = self.fai.index2pos(end_index) + diff_range
                # noinspection PyStringFormat
                raise ValueError("End position [%s, %d] must be within DIFF range [%s, %d]:[%s, %d]." % args)
            return start_index, end_index
        
        else:
            # apply diff range
            return diff_start_index, diff_end_index
    
    def bam_checksum(self):
        if self._bam_checksum is None:
            self._bam_checksum = calc_checksum(self.bam_file.filename)
        
        return self._bam_checksum
    
    def validate_checksum(self, check_sum):
        if self.bam_checksum() != check_sum:
            print(bin2hex(self.bam_checksum()))
            print(bin2hex(check_sum))
            raise ValueError("Invalid checksum.")
    
    def unmutate(
            self,
            diff_file,
            out_bam_file,
            start_ref_name=None,
            start_ref_pos=None,
            end_ref_name=None,
            end_ref_pos=None,
    ):
        """
        Unmutate BAM file in range specified by DIFF file or by parameters.
        :param diff_file: binary file
        :param out_bam_file: pysam.AlignmentFile
        :param start_ref_name:
        :param start_ref_pos:
        :param end_ref_name:
        :param end_ref_pos:
        :return:
        """
        self.__init_counters()
        
        Diff.validate(diff_file)
        start_index, end_index = self.__resolve_range(
            diff_file,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
        )
        
        snv_alignments = []
        bam_iter = BamIterator(self.bam_file, start_index, end_index)
        alignment = next(bam_iter)
        
        diff_iter = DiffIterator(diff_file, self.fai, start_index, end_index)
        diff = next(diff_iter)
        
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
                
                # write remaining alignments (in EOF DIFF case only)
                while alignment is not None:
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = next(bam_iter)
                
                break
            
            elif alignment.is_unmapped:
                self.unmapped_counter += 1
                self.__write_alignment(out_bam_file, alignment)
                alignment = next(bam_iter)
            
            elif self.__is_before_index(alignment, diff.index):
                self.__write_alignment(out_bam_file, alignment)
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, diff.index):
                self.__unmutate_overlap(snv_alignments, diff.mut_map)
                # done with this diff, read next
                diff = next(diff_iter)
                if diff is not None:
                    # could be end of DIFF file
                    self.__write_before_index(out_bam_file, snv_alignments, diff.index)
                    self.__set_seq_positions(snv_alignments, diff.ref_pos)
            
            else:  # alignment is overlapping diff
                # find sequence position of vac
                seq_pos = self.ref_pos2seq_pos(alignment, diff.ref_pos)
                snv_alignments.append(SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.overlapping_counter += 1
                # noinspection PyAttributeOutsideInit
                self.max_coverage = max(len(snv_alignments), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.diff_counter = diff_iter.counter
    
    def __unmutate_overlap(self, snv_alignments, mut_map):
        for snv_alignment in snv_alignments:
            if snv_alignment.pos is not None:
                # alignment has vac to mutate
                self.__mutate_alignment(snv_alignment.alignment, snv_alignment.pos, mut_map)
                # done with current vac
                snv_alignment.pos = None
    
    def mutate(self, in_vac_file, out_bam_file, out_diff_file):
        """
        Mutate BAM file SNVs.
        :param in_vac_file: binary file
        input variant allele count file
        input bam file
        :param out_bam_file: pysam.AlignmentFile
        output bam file
        :param out_diff_file: binary file
        output diff file
        """
        self.__init_counters()
        
        snv_alignments = []
        
        bam_iter = BamIterator(self.bam_file)
        alignment = next(bam_iter)
        
        vac_iter = VacIterator(in_vac_file, self.fai)
        vac = next(vac_iter)
        
        Diff.write_header(
            diff_file=out_diff_file,
            bam_checksum=bytes(Diff.MD5_LENGTH),
            start_index=self.fai.first_index(),
            end_index=self.fai.last_index()
        )
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
                    alignment = next(bam_iter)
                
                break
            
            elif alignment.is_unmapped:
                self.unmapped_counter += 1
                self.__write_alignment(out_bam_file, alignment)
                alignment = next(bam_iter)
            
            elif self.__is_before_index(alignment, vac.index):
                self.__write_alignment(out_bam_file, alignment)
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, vac.index):
                self.__mutate_overlap(out_diff_file, snv_alignments, vac)
                # done with this vac, read next
                vac = next(vac_iter)
                if vac is not None:
                    # could be end of VAC file
                    self.__write_before_index(out_bam_file, snv_alignments, vac.index)
                    self.__set_seq_positions(snv_alignments, vac.ref_pos)
            
            else:  # alignment is overlapping vac
                # find sequence position of vac
                seq_pos = self.ref_pos2seq_pos(alignment, vac.ref_pos)
                snv_alignments.append(SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.overlapping_counter += 1
                # noinspection PyAttributeOutsideInit
                self.max_coverage = max(len(snv_alignments), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.snv_counter = vac_iter.counter
    
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
    
    def __write_diff(self, out_diff_file, index, mut_map):
        self.diff_counter += 1
        mut_tuple = (mut_map['A'], mut_map['T'], mut_map['G'], mut_map['C'])
        Diff.write_record(out_diff_file, index, mut_tuple)
    
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


class MutatorCaller:
    SAM_COMMENT_TAG = 'CO'
    MUT_COMMENT_PREFIX = 'MUT:'
    
    STAT_ALIGNMENT_COUNT = 0,
    STAT_UNMAPPED_COUNT = 1,
    STAT_OVERLAPPING_COUNT = 2,
    STAT_MAX_COVERAGE = 3,
    STAT_SNV_COUNT = 4,
    STAT_MUT_COUNT = 5,
    STAT_DIFF_COUNT = 6
    
    def __init__(self, rnd=random.SystemRandom(), verbose=False):
        self.rnd = rnd
        self.verbose = verbose
        self._stats = {}
    
    @classmethod
    def is_mutated(cls, bam_filename):
        with pysam.AlignmentFile(bam_filename, 'rb') as bam_file:
            return cls.__is_mutated(bam_file.header)
    
    @classmethod
    def __is_mutated(cls, header_map):
        if cls.SAM_COMMENT_TAG in header_map:
            for comment in header_map[cls.SAM_COMMENT_TAG]:
                if comment[:len(cls.MUT_COMMENT_PREFIX)] == cls.MUT_COMMENT_PREFIX:
                    return True
        return False
    
    @classmethod
    def __add_comment(cls, header_map, comment):
        if cls.SAM_COMMENT_TAG in header_map:
            header_map[cls.SAM_COMMENT_TAG].append(comment)
        else:
            header_map[cls.SAM_COMMENT_TAG] = [comment]
    
    @classmethod
    def __mut_header(cls, header_map, checksum):
        if cls.__is_mutated(header_map):
            raise ValueError("File appears to be already mutated.")
        else:
            comment = cls.MUT_COMMENT_PREFIX + bin2hex(checksum)
            cls.__add_comment(header_map, comment)
            return header_map
    
    @classmethod
    def __unmut_header(cls, header_map):
        if cls.__is_mutated(header_map):
            comment_list = []
            for i in range(len(header_map[cls.SAM_COMMENT_TAG])):
                comment = header_map[cls.SAM_COMMENT_TAG][i]
                if comment[:len(cls.MUT_COMMENT_PREFIX)] != cls.MUT_COMMENT_PREFIX:
                    comment_list.append(comment)
            
            header_map[cls.SAM_COMMENT_TAG] = comment_list
            return header_map
        else:
            raise ValueError('File does not appear to be mutated.')
    
    def stat(self, stat_id):
        if stat_id in self._stats:
            return self._stats[stat_id]
        else:
            raise ValueError('Stat not found.')
    
    def mutate(
            self,
            bam_filename,
            vac_filename,
            out_bam_filename
    ):
        self._stats = {}
        out_diff_file = io.BytesIO()
        with pysam.AlignmentFile(bam_filename, 'rb') as sam_file:
            mut = Mutator(sam_file, rnd=self.rnd, verbose=self.verbose)
            mut_header = self.__mut_header(sam_file.header, mut.bam_checksum())
            with pysam.AlignmentFile(out_bam_filename, 'wb', header=mut_header) as out_bam_file, \
                    open(vac_filename, 'rb') as vac_file:
                mut.mutate(
                    in_vac_file=vac_file,
                    out_bam_file=out_bam_file,
                    out_diff_file=out_diff_file
                )
            
            self._stats = {
                self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
                self.STAT_UNMAPPED_COUNT: mut.unmapped_counter,
                self.STAT_OVERLAPPING_COUNT: mut.overlapping_counter,
                self.STAT_MAX_COVERAGE: mut.max_coverage,
                self.STAT_SNV_COUNT: mut.snv_counter,
                self.STAT_MUT_COUNT: mut.mut_counter,
                self.STAT_DIFF_COUNT: mut.diff_counter
            }
        
        # TODO remove difference
        # print('before bam ' + bin2hex(calc_checksum(out_bam_file.filename)))
        # bam2sam(out_bam_file.filename, out_bam_file.filename + b'.sam')
        # print('before sam ' + bin2hex(calc_checksum(out_bam_file.filename + b'.sam')))
        # sam2bam(out_bam_file.filename + b'.sam', out_bam_file.filename)
        #
        # print('after bam ' + bin2hex(calc_checksum(out_bam_file.filename)))
        # bam2sam(out_bam_file.filename, out_bam_file.filename + b'.sam')
        # print('after sam ' + bin2hex(calc_checksum(out_bam_file.filename + b'.sam')))
        # exit(0)
        
        Diff.write_checksum(out_diff_file, calc_checksum(out_bam_file.filename))
        out_diff_file.seek(0)
        return out_diff_file
    
    def unmutate(
            self,
            bam_filename,
            diff_file,
            out_bam_filename):
        self._stats = {}
        
        with pysam.AlignmentFile(bam_filename, 'rb') as sam_file, \
                open(diff_file, 'rb') as diff_file:
            unmut_header = self.__unmut_header(sam_file.header)
            
            # print(unmut_header)
            # exit(0)
            
            mut = Mutator(sam_file, rnd=self.rnd, verbose=self.verbose)
            with pysam.AlignmentFile(out_bam_filename, 'wb', header=unmut_header) as out_sam_file:
                mut.unmutate(
                    diff_file=diff_file,
                    out_bam_file=out_sam_file,
                )
            
            self._stats = {
                self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
                self.STAT_UNMAPPED_COUNT: mut.unmapped_counter,
                self.STAT_OVERLAPPING_COUNT: mut.overlapping_counter,
                self.STAT_MAX_COVERAGE: mut.max_coverage,
                self.STAT_MUT_COUNT: mut.mut_counter,
                self.STAT_DIFF_COUNT: mut.diff_counter
            }
