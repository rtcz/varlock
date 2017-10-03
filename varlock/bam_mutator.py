import random
import io

from varlock.fasta_index import FastaIndex
from varlock.mutator import Mutator

import varlock.common as cmn
import varlock.iters as iters
import varlock.bam as bam
import varlock.bdiff as bdiff


# TODO rename MutatorWrapper ?
class BamMutator:
    # written alignments
    STAT_ALIGNMENT_COUNT = 'alignment_count'
    
    # number of alignments covering VAC / DIFF position
    STAT_COVERING_COUNT = 'covering_count'
    
    # number of vac records read (variants)
    STAT_VAC_COUNT = 'vac_count'
    
    # total number of non-synonymously mutated variants (for each alignment)
    STAT_MUT_COUNT = 'mut_count'
    
    # number of DIFF records written / read
    STAT_DIFF_COUNT = 'diff_count'
    
    FROM_INDEX = 'from_index'
    TO_INDEX = 'to_index'
    SECRET = 'secret'
    MB_CHECKSUM = 'mb_checksum'
    
    def __init__(
            self,
            filename: str,
            rnd=random.SystemRandom(),
            verbose: bool = False
    ):
        """
        :param filename: BAM filename
        :param rnd:  random number generator
        :param verbose:
        """
        self._verbose = verbose
        self._rnd = rnd
        self._stats = {}
        self._bam_filename = filename
        
        with bam.open_bam(self._bam_filename, 'rb') as bam_file:
            self._bam_header = bam_file.header
        self._fai = FastaIndex(self._bam_header)
        
        if self._verbose:
            print("Calculating VAC's checksum")
        
        self._bam_checksum = cmn.checksum(self._bam_filename)
    
    def stat(self, stat_id):
        if stat_id in self._stats:
            return self._stats[stat_id]
        else:
            raise ValueError('Stat not found.')
    
    @property
    def stats(self):
        return self._stats
    
    def mutate(
            self,
            vac_filename: str,
            mut_bam_filename: str,
            secret: bytes
    ):
        """
        :param vac_filename:
        :param mut_bam_filename:
        :param secret: Secret key written into DIFF used for unmapped alignment encryption.
        Secret must be of size specified by DIFF format.
        :return diff_file:
        """
        self._stats = {}
        
        vac_checksum = cmn.bytes2hex(cmn.checksum(vac_filename))
        bam_checksum = cmn.bytes2hex(self._bam_checksum)
        header = bam.mut_header(self._bam_header, bam_checksum, vac_checksum)
        
        mut = Mutator(fai=self._fai, rnd=self._rnd, verbose=self._verbose)
        with bam.open_bam(mut_bam_filename, 'wb', header=header) as mut_bam_file, \
                iters.VacIterator(vac_filename, self._fai) as vac_iter, \
                iters.FullBamIterator(self._bam_filename) as bam_iter:
            bdiff_io = mut.mutate(
                mut_bam_file=mut_bam_file,
                vac_iter=vac_iter,
                bam_iter=bam_iter,
                secret=secret
            )
        
        self._stats = {
            self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
            self.STAT_COVERING_COUNT: mut.covering_counter,
            self.STAT_VAC_COUNT: mut.vac_counter,
            self.STAT_MUT_COUNT: mut.mut_counter,
            self.STAT_DIFF_COUNT: mut.diff_counter
        }
        
        # TODO resolve difference between mutated and converted (bam->sam->bam) and bam
        # print('before bam ' + bytes2hex(checksum(out_bam_file._filename)))
        # bam2sam(out_bam_file._filename, out_bam_file._filename + b'.sam')
        # print('before sam ' + bytes2hex(checksum(out_bam_file._filename + b'.sam')))
        # sam2bam(out_bam_file._filename + b'.sam', out_bam_file._filename)
        #
        # print('after bam ' + bytes2hex(checksum(out_bam_file._filename)))
        # bam2sam(out_bam_file._filename, out_bam_file._filename + b'.sam')
        # print('after sam ' + bytes2hex(checksum(out_bam_file._filename + b'.sam')))
        # exit(0)
        
        return bdiff_io.file(header={
            self.FROM_INDEX: self._fai.first_index(),
            self.TO_INDEX: self._fai.last_index(),
            self.MB_CHECKSUM: cmn.bytes2hex(cmn.checksum(mut_bam_filename)),
            self.SECRET: cmn.bytes2hex(secret)
        })
    
    def unmutate(
            self,
            bdiff_file: io.BytesIO,
            out_bam_filename: str,
            start_ref_name: str = None,
            start_ref_pos: int = None,
            end_ref_name: str = None,
            end_ref_pos: int = None,
            include_unmapped: bool = False,
            unmapped_only: bool = False
    ):
        """
        Unmutate BAM file in range specified by DIFF file or by parameters.
        :param bdiff_file:
        :param out_bam_filename:
        :param start_ref_name: inclusive
        :param start_ref_pos: 0-based, inclusive
        :param end_ref_name: inclusive
        :param end_ref_pos: 0-based, inclusive
        :param include_unmapped: Include all unplaced unmapped reads.
        :param unmapped_only: Only unmapped reads - both placed and unplaced.
         Overrides other parameters.
        
        When range is supplied partialy covered reads are also included,
        but only variants within range are unmutated.
        """
        self._stats = {}
        
        with bam.open_bam(self._bam_filename, 'rb') as bam_file:
            header = bam.unmut_header(bam_file.header)
            mut = Mutator(fai=self._fai, rnd=self._rnd, verbose=self._verbose)
            
            with bam.open_bam(out_bam_filename, 'wb', header=header) as out_bam_file:
                bdiff_io = bdiff.BdiffIO(bdiff_file)
                # validate checksum
                # TODO refactor
                if self._bam_checksum != cmn.hex2bytes(bdiff_io.header[self.MB_CHECKSUM]):
                    raise ValueError('BDIFF does not refer to this BAM')
                # TODO user friendly exception on missing bdiff_io header value
                start_index, end_index = self.resolve_range(
                    bdiff_from_index=bdiff_io.header[self.FROM_INDEX],
                    bdiff_to_index=bdiff_io.header[self.TO_INDEX],
                    start_ref_name=start_ref_name,
                    start_ref_pos=start_ref_pos,
                    end_ref_name=end_ref_name,
                    end_ref_pos=end_ref_pos
                )
                # TODO move iterators to with statement
                mut.unmutate(
                    bam_iter=iters.bam_iterator(
                        self._bam_filename,
                        start_index,
                        end_index,
                        unmapped_only,
                        include_unmapped
                    ),
                    bdiff_iter=iters.BdiffIterator(
                        bdiff_io=bdiff_io,
                        fai=self._fai,
                        start_index=start_index,
                        end_index=end_index
                    ),
                    out_bam_file=out_bam_file,
                    secret=cmn.hex2bytes(bdiff_io.header[self.SECRET])
                )
            
            self._stats = {
                self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
                self.STAT_COVERING_COUNT: mut.covering_counter,
                self.STAT_MUT_COUNT: mut.mut_counter,
                self.STAT_DIFF_COUNT: mut.diff_counter
            }
    
    def resolve_range(
            self,
            bdiff_from_index,
            bdiff_to_index,
            start_ref_name: str,
            start_ref_pos: int,
            end_ref_name: str,
            end_ref_pos: int
    ):
        """
        Resolve between DIFF range and user specified range.
        :param bdiff_from_index:
        :param bdiff_to_index:
        :param start_ref_name:
        :param start_ref_pos:
        :param end_ref_name:
        :param end_ref_pos:
        :return: tuple (from_index, to_index)
        """
        if start_ref_name is not None or end_ref_name is not None:
            # use user range
            from_index = self._fai.resolve_start_index(start_ref_name, start_ref_pos)
            to_index = self._fai.resolve_end_index(end_ref_name, end_ref_pos)
            # validate user range
            # TODO show 1-based positions
            if from_index < bdiff_from_index:
                args = (start_ref_name, start_ref_pos) \
                       + self._fai.index2pos(bdiff_from_index) \
                       + self._fai.index2pos(bdiff_to_index)  # type: tuple
                # noinspection PyStringFormat
                raise ValueError("Start position [%s, %d] must be within BDIFF range [%s, %d]:[%s, %d]." % args)
            if to_index > bdiff_to_index:
                args = (end_ref_name, end_ref_pos) \
                       + self._fai.index2pos(bdiff_from_index) \
                       + self._fai.index2pos(bdiff_to_index)  # type: tuple
                # noinspection PyStringFormat
                raise ValueError("End position [%s, %d] must be within BDIFF range [%s, %d]:[%s, %d]." % args)
            return from_index, to_index
        else:
            # use diff range
            return bdiff_from_index, bdiff_to_index
