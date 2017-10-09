import pysam
import io
import os
import struct

from Crypto.Cipher import PKCS1_OAEP
from Crypto.Signature import PKCS1_v1_5
from Crypto.PublicKey import RSA
from Crypto.Hash import MD5

from varlock.bam import open_bam
from varlock.common import open_vcf
from varlock.bdiff import BdiffIO
from varlock.fasta_index import FastaIndex
from varlock.vac import Vac
from varlock.bam_mutator import BamMutator
from varlock.mutator import Mutator
from varlock.aes import FileAES
from varlock import common as cmn


class Varlocker:
    """
    Main class of varlock package for reversible depersonalization of BAM file
    by altering contained SNVs and INDELs based on supplied VAC (VCF based) file.
    """
    AES_KEY_LENGTH = 32
    
    # key used in stream cipher encryption of unmapped alignments
    SECRET_KEY_LENGTH = 16
    
    @staticmethod
    def create_vac(bam_filename: str, vcf_filename: str, out_vac_filename: str, verbose=False):
        """
        BAM and VCF should use same reference genome.
        VCF must contain INFO column with sub-fields AC and AN.
        :param bam_filename:
        :param vcf_filename:
        :param out_vac_filename:
        :param verbose:
        :return:
        """
        # TODO verify that BAM and VCF share the same reference genome ?
        with open_vcf(vcf_filename, 'rt') as vcf_file, \
                open_bam(bam_filename, 'rb') as sam_file, \
                open(out_vac_filename, 'wb') as out_vac_file:
            vac = Vac(FastaIndex(sam_file.header), verbose)
            vac.vcf2vac(vcf_file, out_vac_file)
    
    @classmethod
    def encrypt(
            cls,
            rsa_sign_key: RSA,
            rsa_enc_key: RSA,
            bam_filename: str,
            vac_filename: str,
            out_bam_filename: str,
            out_enc_diff_filename: str,
            verbose=False
    ):
        """
        Mutate BAM file and store it along with encrypted DIFF file.
        !!! BAM and VAC (VCF derivate file) should use the same reference genome !!!
        Output formats:
        .mut.bam
        .diff.enc
        :param rsa_sign_key: private key to sign DIFF
        :param rsa_enc_key: public key to encrypt DIFF (contained AES key)
        :param bam_filename: input bam
        :param vac_filename: VAC file
        :param out_bam_filename: output bam
        :param out_enc_diff_filename: output binary file
        :param verbose:
        """
        aes_key = os.urandom(cls.AES_KEY_LENGTH)
        mut = BamMutator(bam_filename, verbose=verbose)
        
        with open(out_enc_diff_filename, 'wb') as enc_diff_file:
            if verbose:
                print('--- Mutating BAM ---')
            
            diff_file = mut.mutate(
                vac_filename,
                out_bam_filename,
                os.urandom(cls.SECRET_KEY_LENGTH)
            )
            
            if verbose:
                print('stats: %s' % mut.stats)
            
            signature = cls.__sign(diff_file, rsa_sign_key, verbose)
            cls.__write_aes_key(enc_diff_file, aes_key, rsa_enc_key)
            cls.__write_signature(enc_diff_file, signature)
            cls.__encrypt(diff_file, aes_key, enc_diff_file, verbose)
    
    @classmethod
    def decrypt(
            cls,
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
            rsa_ver_key: RSA = None,
            verbose: bool = False
    ):
        """
        Reverse of mutate operation. Restore original BAM file.
        Output formats:
        .bam
        :param rsa_key: An RSA key object (`RsaKey`) containing the private key
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
        :param verbose:
        """
        # TODO verify if bam_filename is mutated
        # TODO compare mutated BAM checksum to checksum stored in DIFF header
        
        # make sure that BAM is indexed
        pysam.index(bam_filename)
        with io.BytesIO() as diff_file, \
                open(enc_diff_filename, 'rb') as enc_diff_file:
            aes_key = cls.__read_aes_key(enc_diff_file, rsa_key)
            signature = cls.__read_signature(enc_diff_file)
            cls.__decrypt(enc_diff_file, aes_key, diff_file, verbose)
            cls.__verify(diff_file, signature, rsa_ver_key, verbose)
            
            secret = BdiffIO(diff_file).header.get(BamMutator.BDIFF_SECRET_TAG)
            if (unmapped_only or include_unmapped) and isinstance(secret, bytes):
                raise ValueError('BDIFF must contain secret to decrypt unmapped reads.')
            
            if verbose:
                print('--- Unmutating BAM ---')
            
            # unmutate
            mut = BamMutator(bam_filename, verbose=verbose)
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
    
    @classmethod
    def reencrypt(
            cls,
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
            rsa_ver_key: RSA = None,
            verbose: bool = False
    ):
        """
        Reencrypt DIFF file with supplied public_key.
        Output formats:
        .diff.enc
        :param rsa_key: RSA key with private key for decryption of DIFF and signing new DIFF
        :param rsa_enc_key: RSA key with public key for encryption of new DIFF
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
        :param verbose:
        """
        # TODO verify if bam_filename is mutated
        # make sure that BAM is indexed
        pysam.index(bam_filename)
        with open_bam(bam_filename, 'rb') as sam_file, \
                open(enc_diff_filename, 'rb') as enc_diff_file, \
                io.BytesIO() as diff_file:
            aes_key = cls.__read_aes_key(enc_diff_file, rsa_key)
            signature = cls.__read_signature(enc_diff_file)
            cls.__decrypt(enc_diff_file, aes_key, diff_file, verbose)
            cls.__verify(diff_file, signature, rsa_ver_key, verbose)
            
            bdiff = BdiffIO(diff_file)
            bam_mut = BamMutator(sam_file)
            if bam_mut.checksum != bdiff.header.get(BamMutator.BDIFF_CHECKSUM_TAG):
                # checksum mismatch
                raise ValueError("Provided BDIFF is not associated with this BAM."
                                 " Reason: checksum mismatch.")
            
            # TODO refactor, validate BDIFF header
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
            
            if (unmapped_only or include_unmapped) and BamMutator.BDIFF_SECRET_TAG in bdiff.header:
                raise ValueError('BDIFF must contain secret to decrypt unmapped reads.')
            
            if unmapped_only:
                # out_diff = Diff.truncate(diff_file)
                del bdiff.header[BdiffIO.FROM_INDEX]
                del bdiff.header[BdiffIO.TO_INDEX]
                out_diff = bdiff.file(bdiff.header)
            elif include_unmapped:
                # out_diff = Diff.slice(diff_file, start_index, end_index)
                out_diff = bdiff.file(bdiff.header)
            else:
                # mapped only
                # out_diff = Diff.slice(diff_file, start_index, end_index, False)
                del bdiff.header[BamMutator.BDIFF_SECRET_TAG]
                out_diff = bdiff.file(bdiff.header)
            
            with out_diff, open(out_enc_diff_filename, 'wb') as out_enc_diff_file:
                out_signature = cls.__sign(out_diff, rsa_key, verbose)
                cls.__write_aes_key(out_enc_diff_file, aes_key, rsa_enc_key)
                cls.__write_signature(out_enc_diff_file, out_signature)
                cls.__encrypt(out_diff, aes_key, out_enc_diff_file, verbose)
    
    @classmethod
    def __encrypt(cls, diff, aes_key, enc_diff, verbose):
        if verbose:
            print('--- Encrypting DIFF ---')
        
        diff.seek(0)
        out_aes = FileAES(aes_key)
        out_aes.encrypt(diff, enc_diff)
    
    @classmethod
    def __decrypt(cls, enc_diff, aes_key, diff, verbose):
        if verbose:
            print('--- Decrypting DIFF ---')
        # decrypt diff
        aes = FileAES(aes_key)
        aes.decrypt(enc_diff, diff)
    
    @classmethod
    def __sign(cls, diff, rsa_key: RSA, verbose):
        """
        :param diff:
        :param rsa_key: RSA key with private key for signing DIFF
        :return:
        """
        if verbose:
            print('--- Signing DIFF ---')
        
        # must use Crypto hash object
        hash_obj = MD5.new()
        diff.seek(0)
        hash_obj.update(diff.read())
        diff.seek(0)
        return PKCS1_v1_5.new(rsa_key).sign(hash_obj)
    
    @classmethod
    def __verify(cls, diff, signature, rsa_key: RSA, verbose):
        """
        :param diff:
        :param signature: DIFF signature
        :param rsa_key: RSA key with public key to verify DIFF
        :param verbose:
        """
        if rsa_key is not None:
            # verification step
            if verbose:
                print('--- Verifying DIFF ---')
            
            # must use Crypto hash object
            hash_obj = MD5.new()
            diff.seek(0)
            hash_obj.update(diff.read())
            diff.seek(0)
            if not PKCS1_v1_5.new(rsa_key).verify(hash_obj, signature):
                raise ValueError('DIFF has invalid signature')
    
    @classmethod
    def __read_aes_key(cls, enc_diff, rsa_key: RSA):
        """
        Retrieve AES key decrypted by RSA private key from header of encrypted DIFF
        :param enc_diff: encrypted DIFF file
        :param rsa_key: RSA key with private key for decryption of DIFF
        """
        aes_length = struct.unpack_from('<I', enc_diff.read(4))[0]
        enc_aes_key = enc_diff.read(aes_length)
        return PKCS1_OAEP.new(rsa_key).decrypt(enc_aes_key)
    
    @classmethod
    def __write_aes_key(cls, enc_diff, aes_key: bytes, rsa_key: RSA):
        """
        Write AES key encrypted by RSA public key into header of encrypted DIFF
        :param enc_diff: encrypted DIFF file
        :param aes_key: AES key to write
        :param rsa_key: RSA key with public key for encryption of AES key
        """
        enc_aes_key = PKCS1_OAEP.new(rsa_key).encrypt(aes_key)
        enc_diff.write(struct.pack('<I', len(enc_aes_key)))
        enc_diff.write(enc_aes_key)
    
    @classmethod
    def __read_signature(cls, enc_diff):
        """
        Retrieve decrypted signature (hash) from encrypted DIFF.
        :param enc_diff: encrypted DIFF file
        """
        signature_length = struct.unpack_from('<I', enc_diff.read(4))[0]
        return enc_diff.read(signature_length)
    
    @classmethod
    def __write_signature(cls, enc_diff, signature):
        """
        :param enc_diff:
        """
        enc_diff.write(struct.pack('<I', len(signature)))
        enc_diff.write(signature)

# @classmethod
# def verify(cls, rsa_dec_key: RSA, rsa_ver_key: RSA, enc_diff_filename: str, verbose: False):
#     """
#     :param rsa_dec_key: RSA key with private key to decrypt DIFF
#     :param rsa_ver_key: RSA key with public key to verify DIFF
#     :param enc_diff_filename: diff to reencrypt
#     :param verbose:
#     :return: True if the DIFF was created by authority with the public key
#     """
#     with open(enc_diff_filename, 'rb') as enc_diff_file, \
#             io.BytesIO() as diff_file:
#         # retrieve encrypted AES from DIFF
#         aes_length = struct.unpack_from('<I', enc_diff_file.read(4))[0]
#         enc_aes_key = enc_diff_file.read(aes_length)
#         aes_key = PKCS1_OAEP.new(rsa_dec_key).decrypt(enc_aes_key)
#
#         # retrieve signature from encrypted DIFF
#         signature_length = struct.unpack_from('<I', enc_diff_file.read(4))[0]
#         signature = enc_diff_file.read(signature_length)
#         signed_hash = PKCS1_OAEP.new(rsa_ver_key.publickey()).decrypt(signature)
#
#         if verbose:
#             print('Decrypting DIFF')
#
#         # decrypt diff
#         aes = vrl.FileAES(aes_key)
#         aes.decrypt(enc_diff_file, diff_file)
#
#         if verbose:
#             print('Calculating DIFF hash')
#
#         return cls.__verify(diff_file, signed_hash)

# @classmethod
# def __verify(cls, diff_file, signed_hash):
#     """
#     :param diff_file: (decrypted) DIFF
#     :param signed_hash: decrypted signature
#     :return: True if verified
#     """
#     diff_file.seek(0)
#     hasher = hashlib.md5()
#     hasher.update(diff_file.read())
#
#     return signed_hash == hasher.digest()
