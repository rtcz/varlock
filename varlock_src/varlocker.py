import io
import struct

import pysam
from Crypto.Cipher import PKCS1_OAEP
from Crypto.Hash import MD5
from Crypto.PublicKey import RSA
from Crypto.Signature import PKCS1_v1_5

from varlock_src.aes import FileAES
from varlock_src.bam import open_bam
from varlock_src.bam_mutator import BamMutator
from varlock_src.bdiff import BdiffIO
from varlock_src.common import open_vcf
from varlock_src.fasta_index import FastaIndex
from varlock_src.random import VeryRandom
from varlock_src.vac import Vac

import gzip
import pyfaidx

class Varlocker:
    """
    Main class of varlock package for reversible depersonalization of BAM file
    by altering contained SNVs and INDELs based on supplied VAC (VCF based) file.
    """
    AES_KEY_LENGTH = 32
    
    # key used in stream cipher encryption of unmapped alignments
    SECRET_KEY_LENGTH = 16
    
    def __init__(self, verbose=False):
        self._verbose = verbose
    
    def create_vac(self, bam_filename: str, vcf_filename: str, out_vac_filename: str, ref_fasta_filename: str):
        """
        BAM and VCF should use same reference genome.
        VCF must contain INFO column with sub-fields AC and AN.
        :param bam_filename: filename of the SAM/BAM file, from which the header is extracted
        :param vcf_filename: filename of the input VCF file
        :param out_vac_filename: filename of the output VAC file
        :param ref_fasta_filename: filename to reference FASTA file
        """
        # TODO verify that BAM and VCF share the same reference genome ?

        # load the reference FASTA
        ref_fasta = None
        if ref_fasta_filename is not None:
            if self._verbose:
                print('--- Loading Reference Fasta ---')
            ref_fasta = pyfaidx.Fasta(ref_fasta_filename)

        # is VCF gzipped?
        is_gzipped = vcf_filename.endswith('.gz')

        # open all files and create the VAC file
        if self._verbose:
            print('--- Processing VCF %s (gzipped: %s) ---' % (vcf_filename, is_gzipped))
        with io.BufferedReader(gzip.open(vcf_filename)) if is_gzipped else open(vcf_filename, 'rt') as vcf_file, \
                open_bam(bam_filename, 'rb') as sam_file, open(out_vac_filename, 'wb') as out_vac_file:
            vac = Vac(FastaIndex(sam_file.header), self._verbose)
            vac.vcf2vac(vcf_file, out_vac_file, ref_fasta)
    
    def encrypt(
            self,
            rsa_sign_key: RSA,
            rsa_enc_key: RSA,
            bam_filename: str,
            vac_filename: str,
            out_bam_filename: str,
            out_enc_diff_filename: str,
            mut_p: float,
            seed: str = None
    ):
        """
        Mutate BAM file and store it along with encrypted DIFF file.
        !!! BAM and VAC (VCF derived file) should be derived from the same reference genome !!!
        Output formats:
        .mut.bam
        .diff.enc
        :param rsa_sign_key: private key to sign DIFF
        :param rsa_enc_key: public key to encrypt AES key to encrypt DIFF
        :param bam_filename: input bam
        :param vac_filename: VAC file
        :param out_bam_filename: output bam
        :param out_enc_diff_filename: output binary file
        :param mut_p: random variant (mutation) probability per genome base
        :param seed: random number generator seed
        """
        if self._verbose:
            print('--- Mutating BAM ---')
        
        rnd = VeryRandom(seed)
        aes_key = rnd.rand_bytes(self.AES_KEY_LENGTH)
        mut = BamMutator(filename=bam_filename, verbose=self._verbose)
        
        with open(out_enc_diff_filename, 'wb') as enc_diff_file:
            
            diff_file = mut.mutate(
                vac_filename=vac_filename,
                mut_bam_filename=out_bam_filename,
                secret=rnd.rand_bytes(self.SECRET_KEY_LENGTH),
                mut_p=mut_p,
                rnd=rnd
            )
            
            if self._verbose:
                print('stats: %s' % mut.stats)
            
            signature = self._sign(diff_file, rsa_sign_key)
            self._write_aes_key(enc_diff_file, aes_key, rsa_enc_key)
            self._write_signature(enc_diff_file, signature)
            self._encrypt(diff_file, aes_key, enc_diff_file)
    
    def decrypt(
            self,
            rsa_key: RSA,
            bam_filename: str,
            enc_diff_filename: str,
            out_bam_filename: str,
            start_ref_name: str = None,
            start_ref_pos: int = None,
            end_ref_name: str = None,
            end_ref_pos: int = None,
            include_unmapped: bool = False,
            unmapped_only: bool = False,
            rsa_ver_key: RSA = None
    ):
        """
        Reverse of mutate operation. Restore original BAM file.
        Output formats:
        .bam
        :param rsa_key: private key to decrypt AES key to decrypt DIFF
        :param bam_filename: mutated BAM
        :param enc_diff_filename:
        :param out_bam_filename: output BAM
        :param start_ref_name: inclusive
        :param start_ref_pos: 0-based, inclusive
        :param end_ref_name: inclusive
        :param end_ref_pos: 0-based, inclusive
        :param include_unmapped: Include all unplaced unmapped reads.
        :param unmapped_only: Only unmapped reads - both placed and unplaced.
         Overrides other parameters.
        :param rsa_ver_key: RSA key with public key to verify DIFF
        """
        # TODO verify if bam_filename is mutated
        # TODO compare mutated BAM checksum to checksum stored in DIFF header
        
        if self._verbose:
            print('--- Unmutating BAM ---')
        
        # make sure that BAM is indexed
        pysam.index(bam_filename)
        with io.BytesIO() as diff_file, \
                open(enc_diff_filename, 'rb') as enc_diff_file:
            aes_key = self._read_aes_key(enc_diff_file, rsa_key)
            signature = self._read_signature(enc_diff_file)
            self._decrypt(enc_diff_file, aes_key, diff_file)
            self._verify(diff_file, signature, rsa_ver_key)
            
            # unmutate
            mut = BamMutator(bam_filename, verbose=self._verbose)
            mut.unmutate(
                diff_file,
                out_bam_filename,
                start_ref_name,
                start_ref_pos,
                end_ref_name,
                end_ref_pos,
                include_unmapped,
                unmapped_only
            )
    
    def reencrypt(
            self,
            rsa_key: RSA,
            rsa_enc_key: RSA,
            bam_filename: str,
            enc_diff_filename: str,
            out_enc_diff_filename: str,
            start_ref_name: str = None,
            start_ref_pos: int = None,
            end_ref_name: str = None,
            end_ref_pos: int = None,
            include_unmapped: bool = False,
            unmapped_only: bool = False,
            rsa_ver_key: RSA = None
    ):
        """
        Reencrypt DIFF file with supplied public_key.
        Output formats:
        .diff.enc
        :param rsa_key: private key to decrypt DIFF and sign new DIFF
        :param rsa_enc_key: public key to encrypt new DIFF
        :param bam_filename: mutated BAM
        :param enc_diff_filename: diff to reencrypt
        :param out_enc_diff_filename: reencrypted DIFF
        :param start_ref_name: inclusive
        :param start_ref_pos: 0-based, inclusive
        :param end_ref_name: inclusive
        :param end_ref_pos: 0-based, inclusive
        :param include_unmapped: Include all unplaced unmapped reads.
        :param unmapped_only: Only unmapped reads - both placed and unplaced.
         Overrides other parameters.
        :param rsa_ver_key: RSA key with public key to verify DIFF
        """
        # TODO verify if bam_filename is mutated
        # make sure that BAM is indexed
        pysam.index(bam_filename)
        with open(enc_diff_filename, 'rb') as enc_diff_file, \
                io.BytesIO() as diff_file:
            aes_key = self._read_aes_key(enc_diff_file, rsa_key)
            signature = self._read_signature(enc_diff_file)
            self._decrypt(enc_diff_file, aes_key, diff_file)
            self._verify(diff_file, signature, rsa_ver_key)
            bdiff = BdiffIO(diff_file)
            bam_mut = BamMutator(bam_filename)
            if bam_mut.checksum != bdiff.header.get(BamMutator.BDIFF_CHECKSUM_TAG):
                # checksum mismatch
                raise ValueError("Provided BDIFF is not associated with this BAM."
                                 " Reason: checksum mismatch.")
            
            from_index, to_index = bam_mut.resolve_range(
                bdiff.header[BdiffIO.FROM_INDEX],
                bdiff.header[BdiffIO.TO_INDEX],
                start_ref_name,
                start_ref_pos,
                end_ref_name,
                end_ref_pos
            )
            # use actual effective range
            bdiff.header[BdiffIO.FROM_INDEX] = from_index
            bdiff.header[BdiffIO.TO_INDEX] = to_index
            
            if (unmapped_only or include_unmapped) and BamMutator.BDIFF_SECRET_TAG not in bdiff.header:
                raise ValueError('BDIFF must contain secret to decrypt unmapped reads.')
            
            if unmapped_only:
                del bdiff.header[BdiffIO.FROM_INDEX]
                del bdiff.header[BdiffIO.TO_INDEX]
                out_diff = bdiff.file(bdiff.header)
            elif include_unmapped:
                out_diff = bdiff.file(bdiff.header)
            else:
                # mapped only
                del bdiff.header[BamMutator.BDIFF_SECRET_TAG]
                out_diff = bdiff.file(bdiff.header)
            
            with out_diff, open(out_enc_diff_filename, 'wb') as out_enc_diff_file:
                out_signature = self._sign(out_diff, rsa_key)
                self._write_aes_key(out_enc_diff_file, aes_key, rsa_enc_key)
                self._write_signature(out_enc_diff_file, out_signature)
                self._encrypt(out_diff, aes_key, out_enc_diff_file)
    
    def _encrypt(self, diff, aes_key, enc_diff):
        if self._verbose:
            print('--- Encrypting DIFF ---')
        
        diff.seek(0)
        out_aes = FileAES(aes_key)
        out_aes.encrypt(diff, enc_diff)
    
    def _decrypt(self, enc_diff, aes_key, diff):
        if self._verbose:
            print('--- Decrypting DIFF ---')
        # decrypt diff
        aes = FileAES(aes_key)
        aes.decrypt(enc_diff, diff)
    
    def _sign(self, diff, rsa_key: RSA):
        """
        :param diff:
        :param rsa_key: RSA key with private key for signing DIFF
        :return:
        """
        if self._verbose:
            print('--- Signing DIFF ---')
        
        # must use Crypto hash object
        hash_obj = MD5.new()
        diff.seek(0)
        hash_obj.update(diff.read())
        diff.seek(0)
        return PKCS1_v1_5.new(rsa_key).sign(hash_obj)
    
    def _verify(self, diff, signature, rsa_key: RSA):
        """
        :param diff:
        :param signature: DIFF signature
        :param rsa_key: RSA key with public key to verify DIFF
        """
        if rsa_key is not None:
            # verification step
            if self._verbose:
                print('--- Verifying DIFF ---')
            
            # must use Crypto hash object
            hash_obj = MD5.new()
            diff.seek(0)
            hash_obj.update(diff.read())
            diff.seek(0)
            if not PKCS1_v1_5.new(rsa_key).verify(hash_obj, signature):
                raise ValueError('DIFF has invalid signature')
    
    @staticmethod
    def _read_aes_key(enc_diff, rsa_key: RSA):
        """
        Retrieve AES key decrypted by RSA private key from header of encrypted DIFF
        :param enc_diff: encrypted DIFF file
        :param rsa_key: RSA key with private key for decryption of DIFF
        """
        aes_length = struct.unpack_from('<I', enc_diff.read(4))[0]
        enc_aes_key = enc_diff.read(aes_length)
        return PKCS1_OAEP.new(rsa_key).decrypt(enc_aes_key)
    
    @staticmethod
    def _write_aes_key(enc_diff, aes_key: bytes, rsa_key: RSA):
        """
        Write AES key encrypted by RSA public key into header of encrypted DIFF
        :param enc_diff: encrypted DIFF file
        :param aes_key: AES key to write
        :param rsa_key: RSA key with public key for encryption of AES key
        """
        enc_aes_key = PKCS1_OAEP.new(rsa_key).encrypt(aes_key)
        enc_diff.write(struct.pack('<I', len(enc_aes_key)))
        enc_diff.write(enc_aes_key)
    
    @staticmethod
    def _read_signature(enc_diff):
        """
        Retrieve decrypted signature (hash) from encrypted DIFF.
        :param enc_diff: encrypted DIFF file
        """
        signature_length = struct.unpack_from('<I', enc_diff.read(4))[0]
        return enc_diff.read(signature_length)
    
    @staticmethod
    def _write_signature(enc_diff, signature):
        """
        :param enc_diff:
        """
        enc_diff.write(struct.pack('<I', len(signature)))
        enc_diff.write(signature)
