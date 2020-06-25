import io

import src.bam as bam
import src.common as cmn
import src.iters as iters
from src.bdiff import BdiffIO
from src.fasta_index import FastaIndex
from src.mutator import Mutator
# TODO rename MutatorWrapper ?
from src.random import VeryRandom


class BamMutator:
    # written alignments
    STAT_ALIGNMENT_COUNT = 'alignment_count'
    
    # number of alignments covering VAC / DIFF position
    STAT_COVERING_COUNT = 'covering_count'
    
    # number of vac records read (variants)
    STAT_VAC_COUNT = 'vac_count'
    
    # total number of altered (non-synonymously mutated) variants (for each alignment)
    STAT_MUT_COUNT = 'mut_count'
    
    # number of altered alignments (not including unmapped)
    STAT_ALIGNMENT_MUT_COUNT = 'alignment_mut_count'
    
    # number of DIFF records written / read
    STAT_DIFF_COUNT = 'diff_count'
    
    # secret key used in unmapped alignment encryption
    BDIFF_SECRET_TAG = 'secret'
    
    # mutated BAM's checksum
    BDIFF_CHECKSUM_TAG = 'mb_checksum'
    
    def __init__(
            self,
            filename: str,
            verbose: bool = False
    ):
        """
        :param filename: BAM filename
        :param verbose:
        """
        self._verbose = verbose
        self._stats = {}
        self._bam_filename = filename
        
        with bam.open_bam(self._bam_filename, 'rb') as bam_file:
            self._bam_header = bam_file.header
        
        self._fai = FastaIndex(self._bam_header)
        
        self._checksum = None
    
    def stat(self, stat_id):
        if stat_id in self._stats:
            return self._stats[stat_id]
        else:
            raise ValueError('Stat not found.')
    
    # @property
    # def fai(self):
    #     return self._fai
    
    @property
    def stats(self):
        return self._stats
    
    @property
    def checksum(self):
        """
        :return: BAM's checksum as hex string
        """
        if self._checksum is None:
            if self._verbose:
                print("Calculating BAM's checksum")
            self._checksum = cmn.checksum(self._bam_filename)
        
        return self._checksum
    
    @staticmethod
    def extract_secret_bytes(bdiff_io: BdiffIO):
        """
        Extract secret bytes from BDIFF file.
        :param bdiff_io:
        :return:
        """
        secret_bytes = None
        secret_hex = bdiff_io.header.get(BamMutator.BDIFF_SECRET_TAG)
        if secret_hex is not None:
            secret_bytes = cmn.hex2bytes(secret_hex)
        
        return secret_bytes
    
    def mutate(
            self,
            vac_filename: str,
            mut_bam_filename: str,
            secret: bytes,
            mut_p: float,
            rnd: VeryRandom
    ):
        """
        :param vac_filename:
        :param mut_bam_filename:
        :param secret: Secret key written into DIFF used for unmapped alignment encryption.
        :param mut_p: random variant (mutation) probability per genome base
        :param rnd:  random number generator
        :return diff_file:
        """
        self._stats = {}
        
        header = bam.mut_header(self._bam_header, self.checksum, cmn.checksum(vac_filename))
        with bam.open_bam(mut_bam_filename, 'wb', header=header) as mut_bam_file, \
                iters.VariantIterator(vac_filename, self._fai, mut_p, rnd) as vac_iter, \
                iters.FullBamIterator(self._bam_filename) as bam_iter:
            mut = Mutator(fai=self._fai, verbose=self._verbose)
            bdiff_io = mut.mutate(
                mut_bam_file=mut_bam_file,
                variant_iter=vac_iter,
                bam_iter=bam_iter,
                secret=secret,
                rnd=rnd
            )
        
        self._stats = {
            self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
            self.STAT_COVERING_COUNT: mut.covering_counter,
            self.STAT_VAC_COUNT: mut.vac_counter,
            self.STAT_MUT_COUNT: mut.mut_counter,
            self.STAT_DIFF_COUNT: mut.diff_counter,
            self.STAT_ALIGNMENT_MUT_COUNT: mut.alignment_mut_counter
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
            BdiffIO.FROM_INDEX: self._fai.first_index(),
            BdiffIO.TO_INDEX: self._fai.last_index(),
            self.BDIFF_CHECKSUM_TAG: cmn.checksum(mut_bam_filename),
            self.BDIFF_SECRET_TAG: cmn.bytes2hex(secret)
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
            mut = Mutator(fai=self._fai, verbose=self._verbose)
            
            with bam.open_bam(out_bam_filename, 'wb', header=header) as out_bam_file:
                bdiff_io = BdiffIO(bdiff_file)
                if self._verbose:
                    print('SNV diff count %d' % bdiff_io.snv_count)
                    print('INDEL diff count %d' % bdiff_io.indel_count)
                
                secret = self.extract_secret_bytes(bdiff_io)
                if (include_unmapped or unmapped_only) and secret is None:
                    raise ValueError('BDIFF must contain secret to decrypt unmapped reads.')
                
                # validate checksum
                if self.checksum != bdiff_io.header[self.BDIFF_CHECKSUM_TAG]:
                    raise ValueError('BDIFF does not refer to this BAM')
                
                # TODO user friendly exception on missing bdiff_io header value
                start_index, end_index = self.resolve_range(
                    bdiff_from_index=bdiff_io.header[BdiffIO.FROM_INDEX],
                    bdiff_to_index=bdiff_io.header[BdiffIO.TO_INDEX],
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
                    secret=secret
                )
            
            self._stats = {
                self.STAT_ALIGNMENT_COUNT: mut.alignment_counter,
                self.STAT_COVERING_COUNT: mut.covering_counter,
                self.STAT_MUT_COUNT: mut.mut_counter,
                self.STAT_DIFF_COUNT: mut.diff_counter,
                self.STAT_ALIGNMENT_MUT_COUNT: mut.alignment_mut_counter
            }
    
    # TODO refactor
    def resolve_range(
            self,
            bdiff_from_index: int,
            bdiff_to_index: int,
            start_ref_name: str,
            start_ref_pos: int,
            end_ref_name: str,
            end_ref_pos: int
    ):
        """
        Resolve between the DIFF range and a range specified by user.
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
