import io
import random

import pysam

from .common import bin2hex, calc_checksum
from .diff import Diff
from .mutator import Mutator


class BamMutator:
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
        
        # TODO resolve difference between mutated and converted (bam->sam->bam) bam
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
            out_bam_filename,
            start_ref_name=None,
            start_ref_pos=None,
            end_ref_name=None,
            end_ref_pos=None,
    ):
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
                    start_ref_name=start_ref_name,
                    start_ref_pos=start_ref_pos,
                    end_ref_name=end_ref_name,
                    end_ref_pos=end_ref_pos,
                )
            
            self._stats = {
                self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
                self.STAT_UNMAPPED_COUNT: mut.unmapped_counter,
                self.STAT_OVERLAPPING_COUNT: mut.overlapping_counter,
                self.STAT_MAX_COVERAGE: mut.max_coverage,
                self.STAT_MUT_COUNT: mut.mut_counter,
                self.STAT_DIFF_COUNT: mut.diff_counter
            }
