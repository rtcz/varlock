import random

import numpy as np

from varlock.common import *
from varlock.diff import Diff
from varlock.fasta_index import FastaIndex
from varlock.iterator import DiffIterator, VacIterator, UnmappedOnlyBamIterator, UnmappedBamIterator, BaiBamIterator
from varlock.iterator import FullBamIterator
from varlock.po import SnvAlignment


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
        
        self.prev_alignment = None
    
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
    
    def __write_alignment(self, bam_file, alignment):
        """
        :param alignment: pysam.AlignedSegment
        :return:
        """
        are_placed = self.prev_alignment is not None and is_placed_alignment(
            alignment) and is_placed_alignment(self.prev_alignment)
        # FIXME is case of placed after unplaced alignment valid ?
        if are_placed and alignment.reference_start < self.prev_alignment.reference_start:
            # safety check
            raise IndexError('Unordered write: %s after %s' % (
                self.alignment2str(alignment),
                self.alignment2str(self.prev_alignment)
            ))
        
        bam_file.write(alignment)
        self.alignment_counter += 1
        self.prev_alignment = alignment
        
        if self.verbose and self.alignment_counter % 10000 == 0:
            print("%d alignments processed" % self.alignment_counter)
    
    def __init_counters(self):
        self.alignment_counter = 0  # written alignments
        self.diff_counter = 0  # processed diff records
        self.mut_counter = 0  # mutations pre alignment
        self.snv_counter = 0  # read vac records (SNVs)
        
        self.overlapping_counter = 0  # overlapping alignments
        self.max_coverage = 0  # maximum alignments overlapping single SNV
    
    def resolve_diff_range(
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
            
            diff_range = self.fai.index2pos(diff_start_index) + self.fai.index2pos(diff_end_index)  # type: tuple
            
            if start_index < diff_start_index:
                args = self.fai.index2pos(start_index) + diff_range  # type: tuple
                # noinspection PyStringFormat
                raise ValueError("Start position [%s, %d] must be within DIFF range [%s, %d]:[%s, %d]." % args)
            if end_index > diff_end_index:
                args = self.fai.index2pos(end_index) + diff_range  # type: tuple
                # noinspection PyStringFormat
                raise ValueError("End position [%s, %d] must be within DIFF range [%s, %d]:[%s, %d]." % args)
            return start_index, end_index
        
        else:
            # apply diff range
            return diff_start_index, diff_end_index
    
    def bam_checksum(self):
        if self._bam_checksum is None:
            if self.verbose:
                print("Calculating BAM's checksum.")
            
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
            start_ref_name: str = None,
            start_ref_pos: int = None,
            end_ref_name: str = None,
            end_ref_pos: int = None,
            include_unmapped: bool = False,
            unmapped_only: bool = False
    ):
        """
        Unmutate BAM file in range specified by DIFF file or by parameters.
        :param diff_file: binary file
        :param out_bam_file: pysam.AlignmentFile
        :param start_ref_name: inclusive
        :param start_ref_pos: 0-based, inclusive
        :param end_ref_name: inclusive
        :param end_ref_pos: 0-based, inclusive
        :param include_unmapped: Include all unplaced unmapped reads.
        :param unmapped_only: Only unmapped reads - both placed and unplaced.
         Overrides other parameters.
        """
        self.__init_counters()
        
        Diff.validate(diff_file)
        start_index, end_index = self.resolve_diff_range(
            diff_file,
            start_ref_name,
            start_ref_pos,
            end_ref_name,
            end_ref_pos
        )
        
        alignment_queue = []
        
        if unmapped_only:
            bam_iter = UnmappedOnlyBamIterator(self.bam_file)
        elif include_unmapped:
            bam_iter = UnmappedBamIterator(self.bam_file, start_index, end_index)
        else:
            bam_iter = BaiBamIterator(self.bam_file, start_index, end_index)
        
        alignment = next(bam_iter)
        diff_iter = DiffIterator(diff_file, self.fai, start_index, end_index)
        diff = next(diff_iter)
        
        if self.verbose:
            print('first diff: %s' % diff)
            print('first alignment: %s' % self.alignment2str(alignment))
        
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
                    self.__unmutate_overlap(alignment_queue, diff.mut_map)
                
                for snv_alignment in alignment_queue:
                    self.__write_alignment(out_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF DIFF case only)
                while alignment is not None:
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = next(bam_iter)
                
                if self.verbose:
                    print('last alignment: %s' % self.alignment2str(self.prev_alignment))
                    print("total of %d alignments processed" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, diff.index):
                if len(alignment_queue) == 0:
                    # good to go, all preceeding alignments are written
                    self.__write_alignment(out_bam_file, alignment)
                else:
                    # append to queue
                    alignment_queue.append(SnvAlignment(alignment, None))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, diff.index):
                self.__unmutate_overlap(alignment_queue, diff.mut_map)
                # done with this diff, read next
                prev_diff = diff
                diff = next(diff_iter)
                if diff is not None:
                    # could be end of DIFF file
                    self.__write_done(out_bam_file, alignment_queue, diff.index)
                    self.__set_seq_positions(alignment_queue, diff.ref_pos)
                elif self.verbose:
                    print('last diff: %s' % prev_diff)
            
            else:  # alignment is overlapping diff
                # find sequence position of vac
                seq_pos = ref_pos2seq_pos(alignment, diff.ref_pos)
                alignment_queue.append(SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.overlapping_counter += 1
                # noinspection PyAttributeOutsideInit
                # TODO this is not the real coverage
                self.max_coverage = max(len(alignment_queue), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.diff_counter = diff_iter.counter
    
    def __unmutate_overlap(self, alignment_queue, mut_map):
        for snv_alignment in alignment_queue:
            if snv_alignment.pos is not None:
                # alignment has vac to mutate
                self.__mutate_alignment(snv_alignment.alignment, snv_alignment.pos, mut_map)
                # done with current vac
                snv_alignment.pos = None
    
    def alignment2str(self, alignment):
        if is_placed_alignment(alignment):
            ref_name = alignment.reference_name
            ref_start = alignment.reference_start
            # unplaced alignment
            return '#%d %s:%d' % (self.fai.pos2index(ref_name, ref_start), ref_name, ref_start)
        else:
            return alignment.query_sequence
    
    def mutate(self, in_vac_file, out_bam_file, out_diff_file):
        """
        Mutate BAM file SNVs.
        :param in_vac_file: binary file
        input variant allele count file
        :param out_bam_file: pysam.AlignmentFile
        output bam file
        :param out_diff_file: binary file
        output diff file
        """
        self.__init_counters()
        
        alignment_queue = []
        
        bam_iter = FullBamIterator(self.bam_file)
        alignment = next(bam_iter)
        
        vac_iter = VacIterator(in_vac_file, self.fai)
        vac = next(vac_iter)
        
        if self.verbose:
            print('first vac: %s' % vac)
            print('first alignment: %s' % self.alignment2str(alignment))
        
        Diff.write_header(
            diff_file=out_diff_file,
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
                    self.__mutate_overlap(out_diff_file, alignment_queue, vac)
                
                for snv_alignment in alignment_queue:
                    self.__write_alignment(out_bam_file, snv_alignment.alignment)
                
                # write remaining alignments (in EOF VAC case only)
                while alignment is not None:
                    self.__write_alignment(out_bam_file, alignment)
                    alignment = next(bam_iter)
                
                if self.verbose:
                    print('last alignment: %s' % self.alignment2str(self.prev_alignment))
                    print("total of %d alignments processed" % self.alignment_counter)
                
                break
            
            elif alignment.is_unmapped or self.__is_before_index(alignment, vac.index):
                if len(alignment_queue) == 0:
                    # good to go, all preceeding alignments are written
                    self.__write_alignment(out_bam_file, alignment)
                else:
                    # append to queue
                    alignment_queue.append(SnvAlignment(alignment, None))
                
                alignment = next(bam_iter)
            
            elif self.__is_after_index(alignment, vac.index):
                self.__mutate_overlap(out_diff_file, alignment_queue, vac)
                # done with this vac, read next
                prev_vac = vac
                vac = next(vac_iter)
                if vac is not None:
                    # could be end of VAC file
                    self.__write_done(out_bam_file, alignment_queue, vac.index)
                    self.__set_seq_positions(alignment_queue, vac.ref_pos)
                elif self.verbose:
                    print('last vac: %s' % prev_vac)
            
            else:  # alignment is overlapping vac
                # find sequence position of vac
                seq_pos = ref_pos2seq_pos(alignment, vac.ref_pos)
                alignment_queue.append(SnvAlignment(alignment, seq_pos))
                alignment = next(bam_iter)
                
                self.overlapping_counter += 1
                # noinspection PyAttributeOutsideInit
                # TODO this is not the real coverage
                self.max_coverage = max(len(alignment_queue), self.max_coverage)
        
        # noinspection PyAttributeOutsideInit
        self.snv_counter = vac_iter.counter
    
    def __mutate_overlap(self, out_diff_file, alignment_queue, vac):
        """
        Mutate alignments at current SNV position by allele frequencies of the SNV
        :param alignment_queue: alignments together with vac positions
        :param vac: SNV
        :return: True if NS mutation has occured
        """
        is_mutated = False
        # get pileup of bases from alignments at vac mapping position
        base_pileup = get_base_pileup(alignment_queue)
        alt_ac = self.count_bases(base_pileup)
        
        mut_map = self.create_mut_map(alt_ac=alt_ac, ref_ac=vac.ac, rnd=self.rnd)
        
        for snv_alignment in alignment_queue:
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
        snv_base = get_base(alignment, seq_pos)
        mut_base = mut_map[snv_base]
        
        if snv_base != mut_base:
            # base has been mutated to another base
            self.mut_counter += 1
            is_mutated = True
            set_base(alignment, seq_pos, mut_base)
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
    
    def __write_done(self, out_bam_file, alignment_queue, index):
        """
        Write snv alignments while their end reference position is before index.
        :param out_bam_file:
        :param alignment_queue:
        :param index:
        :return:
        """
        tmp_alignment_queue = alignment_queue[:]
        for snv_alignment in tmp_alignment_queue:
            is_mapped = not snv_alignment.alignment.is_unmapped
            if is_mapped and not self.__is_before_index(snv_alignment.alignment, index):
                # break to keep alignments in order
                # otherwise shorter alignment could be written before longer alignment
                break
            
            self.__write_alignment(out_bam_file, snv_alignment.alignment)
            # remove written alignment
            alignment_queue.remove(snv_alignment)
    
    @staticmethod
    def __set_seq_positions(alignment_queue, ref_pos):
        """
        Set SNV position derived from reference position for each alignment.
        :param alignment_queue:
        :param ref_pos:
        :return:
        """
        for snv_alignment in alignment_queue:
            # TODO optimalize, skip when ref_pos is greater than alignment reference end
            snv_pos = ref_pos2seq_pos(snv_alignment.alignment, ref_pos)
            if snv_pos is not None:
                # alignment is mapped at another new_snv position
                snv_alignment.pos = snv_pos
    
    @classmethod
    def count_bases(cls, base_pileup):
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
    
    @classmethod
    def create_mut_map(cls, alt_ac, ref_ac, rnd):
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
