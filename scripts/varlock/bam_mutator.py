import random

import pysam
from pysam.calignmentfile import VALID_HEADER_TYPES, KNOWN_HEADER_FIELDS

from .common import calc_checksum, bin2hex
from .diff import Diff
from .mutator import Mutator


class BamMutator:
    MUT_TAG = 'mt'
    BAM_TAG = 'bm'
    VAC_TAG = 'vc'
    
    STAT_ALIGNMENT_COUNT = 'alignment_count'
    STAT_OVERLAPPING_COUNT = 'overlapping_count'
    STAT_MAX_COVERAGE = 'max_coverage'
    STAT_SNV_COUNT = 'snv_count'
    STAT_MUT_COUNT = 'mut_count'
    STAT_DIFF_COUNT = 'diff_count'
    
    def __init__(self, rnd=random.SystemRandom(), verbose=False):
        self.rnd = rnd
        self.verbose = verbose
        self._stats = {}
    
    @classmethod
    def __mut_header(cls, header_map, bam_checksum, vac_checksum):
        if cls.MUT_TAG in header_map:
            raise ValueError("File appears to be already mutated.")
        else:
            header_map[cls.MUT_TAG] = {cls.BAM_TAG: bam_checksum, cls.VAC_TAG: vac_checksum}
            return header_map
    
    @classmethod
    def __unmut_header(cls, header_map):
        if cls.MUT_TAG in header_map:
            del header_map[cls.MUT_TAG]
            return header_map
        else:
            raise ValueError('File does not appear to be mutated.')
    
    def stat(self, stat_id):
        if stat_id in self._stats:
            return self._stats[stat_id]
        else:
            raise ValueError('Stat not found.')
    
    def all_stats(self):
        return self._stats
    
    def mutate(
            self,
            bam_filename: str,
            vac_filename: str,
            mut_bam_filename: str,
            diff_file
    ):
        self._stats = {}
        with pysam.AlignmentFile(bam_filename, 'rb') as sam_file:
            mut = Mutator(sam_file, rnd=self.rnd, verbose=self.verbose)
            
            if self.verbose:
                print("Calculating VAC's checksum")
            
            vac_checksum = bin2hex(calc_checksum(vac_filename))
            bam_checksum = bin2hex(mut.bam_checksum())
            mut_header = self.__mut_header(sam_file.header, bam_checksum, vac_checksum)
            
            with pysam.AlignmentFile(mut_bam_filename, 'wb', header=mut_header) as out_bam_file, \
                    open(vac_filename, 'rb') as vac_file:
                mut.mutate(
                    in_vac_file=vac_file,
                    out_bam_file=out_bam_file,
                    out_diff_file=diff_file
                )
            
            self._stats = {
                self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
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
        
        Diff.write_checksum(diff_file, calc_checksum(out_bam_file.filename))
        diff_file.seek(0)
    
    def unmutate(
            self,
            bam_filename: str,
            diff_file: object,
            out_bam_filename: str,
            start_ref_name: str = None,
            start_ref_pos: int = None,
            end_ref_name: str = None,
            end_ref_pos: int = None
    ):
        self._stats = {}
        
        # modify pysam's header format
        VALID_HEADER_TYPES[self.MUT_TAG] = dict
        KNOWN_HEADER_FIELDS[self.MUT_TAG] = {self.BAM_TAG: str, self.VAC_TAG: str}
        
        with pysam.AlignmentFile(bam_filename, 'rb') as sam_file:
            unmut_header = self.__unmut_header(sam_file.header)
            mut = Mutator(sam_file, rnd=self.rnd, verbose=self.verbose)
            
            with pysam.AlignmentFile(out_bam_filename, 'wb', header=unmut_header) as out_sam_file:
                mut.unmutate(
                    diff_file=diff_file,
                    out_bam_file=out_sam_file,
                    start_ref_name=start_ref_name,
                    start_ref_pos=start_ref_pos,
                    end_ref_name=end_ref_name,
                    end_ref_pos=end_ref_pos
                )
            
            self._stats = {
                self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
                self.STAT_OVERLAPPING_COUNT: mut.overlapping_counter,
                self.STAT_MAX_COVERAGE: mut.max_coverage,
                self.STAT_MUT_COUNT: mut.mut_counter,
                self.STAT_DIFF_COUNT: mut.diff_counter
            }
