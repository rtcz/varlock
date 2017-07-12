import io
import os
import struct

from Crypto.Cipher import PKCS1_OAEP
from Crypto.PublicKey import RSA

import varlock as vrl
from .bam import open_bam
from .common import open_vcf
from .diff import Diff
from .fasta_index import FastaIndex
from .vac import Vac


class Varlocker:
    RSA_KEY_LENGTH = 4096
    RSA_CIPHER_LENGTH = 512
    AES_KEY_LENGTH = 32
    
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
        with open_vcf(vcf_filename, 'rt') as vcf_file, \
                open_bam(bam_filename, 'rb') as sam_file, \
                open(out_vac_filename, 'wb') as out_vac_file:
            vac = Vac(FastaIndex(sam_file, keep_chr=False), verbose)
            vac.vcf2vac(vcf_file, out_vac_file)
    
    @classmethod
    def encrypt(
            cls,
            rsa_key: RSA,
            bam_filename: str,
            vac_filename: str,
            out_bam_filename: str,
            out_enc_diff_filename: str,
            verbose=False
    ):
        """
        Mutate BAM file and store it along with encrypted DIFF file.
        BAM and VAC (vcf derivate file) should use same reference genome.
        Output formats:
        .mut.bam
        .diff.enc
        :param rsa_key: An RSA key object (`RsaKey`) containing the public key
        :param bam_filename: input bam
        :param vac_filename: VAC file
        :param out_bam_filename: output bam
        :param out_enc_diff_filename: output binary file
        :param verbose:
        """
        aes_key = os.urandom(cls.AES_KEY_LENGTH)
        mut = vrl.BamMutator(verbose=verbose)
        
        if verbose:
            print('Mutating BAM')
        
        with io.BytesIO() as diff_file, open(out_enc_diff_filename, 'wb') as enc_diff_file:
            
            mut.mutate(bam_filename, vac_filename, out_bam_filename, diff_file)
            # store encrypted AES key in encrypted DIFF
            enc_aes_key = PKCS1_OAEP.new(rsa_key.publickey()).encrypt(aes_key)
            enc_diff_file.write(struct.pack('<I', len(enc_aes_key)))
            enc_diff_file.write(enc_aes_key)
            
            if verbose:
                print('Encrypting DIFF')
            
            # encrypt diff with new AES key
            # noinspection PyTypeChecker
            aes = vrl.FileAES(aes_key)
            aes.encrypt(diff_file, enc_diff_file)
    
    @staticmethod
    def decrypt(
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
        :param verbose:
        """
        with io.BytesIO() as diff_file, \
                open(enc_diff_filename, 'rb') as enc_diff_file:
            # retrieve encrypted AES key stored in encrypted DIFF
            aes_length = struct.unpack_from('<I', enc_diff_file.read(4))[0]
            enc_aes_key = enc_diff_file.read(aes_length)
            aes_key = PKCS1_OAEP.new(rsa_key).decrypt(enc_aes_key)
            
            if verbose:
                print('Decrypting DIFF')
            
            # decrypt diff
            aes = vrl.FileAES(aes_key)
            aes.decrypt(enc_diff_file, diff_file)
            
            secret = Diff.read_header(diff_file)[3]
            if (unmapped_only or include_unmapped) and secret == Diff.SECRET_PLACEHOLDER:
                raise ValueError('Attempt to decrypt unmapped reads without secret key')
            
            if verbose:
                print('Unmutating BAM')
            
            # unmutate
            mut = vrl.BamMutator(verbose=verbose)
            mut.unmutate(
                bam_filename,
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
            rsa_dec_key,
            rsa_enc_key,
            bam_filename: str,
            enc_diff_filename: str,
            out_enc_diff_filename: str,
            start_ref_name: str = None,
            start_ref_pos: int = None,
            end_ref_name: str = None,
            end_ref_pos: int = None,
            include_unmapped: bool = False,
            unmapped_only: bool = False,
            verbose: bool = True
    ):
        """
        Reencrypt DIFF file with supplied public_key.
        Output formats:
        .diff.enc
        :param rsa_dec_key: RSA key for decryption containing the private key
        :param rsa_enc_key: RSA key for encryption containing the public key
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
        :param verbose:
        """
        with open_bam(bam_filename, 'rb') as sam_file, \
                open(enc_diff_filename, 'rb') as enc_diff_file, \
                io.BytesIO() as diff_file:
            # retrieve encrypted AES key stored in encrypted DIFF
            aes_length = struct.unpack_from('<I', enc_diff_file.read(4))[0]
            enc_aes_key = enc_diff_file.read(aes_length)
            aes_key = PKCS1_OAEP.new(rsa_dec_key).decrypt(enc_aes_key)
            
            if verbose:
                print('Decrypting DIFF')
            
            # decrypt diff
            aes = vrl.FileAES(aes_key)
            aes.decrypt(enc_diff_file, diff_file)
            
            mut = vrl.Mutator(sam_file)
            
            # TODO refactor (diff file instance)
            Diff.validate(diff_file)
            checksum, diff_start_index, diff_end_index, secret = Diff.read_header(diff_file)
            mut.validate_checksum(checksum)
            start_index, end_index = mut.resolve_range(
                diff_start_index,
                diff_end_index,
                start_ref_name,
                start_ref_pos,
                end_ref_name,
                end_ref_pos
            )
            
            if (unmapped_only or include_unmapped) and secret == Diff.SECRET_PLACEHOLDER:
                raise ValueError('Attempt to reencrypt unmapped reads without secret key')
            
            if unmapped_only:
                new_diff = Diff.truncate(diff_file)
            elif include_unmapped:
                new_diff = Diff.slice(diff_file, start_index, end_index)
            else:
                # mapped only
                new_diff = Diff.slice(diff_file, start_index, end_index, False)
            
            with new_diff, open(out_enc_diff_filename, 'wb') as out_enc_diff_file:
                # generate new AES key
                aes_key = os.urandom(cls.AES_KEY_LENGTH)
                # store encrypted AES key in encrypted DIFF
                enc_aes_key = PKCS1_OAEP.new(rsa_enc_key.publickey()).encrypt(aes_key)
                out_enc_diff_file.write(struct.pack('<I', len(enc_aes_key)))
                out_enc_diff_file.write(enc_aes_key)
                if verbose:
                    print('Encrypting DIFF')
                
                # noinspection PyTypeChecker
                out_aes = vrl.FileAES(aes_key)
                out_aes.encrypt(new_diff, out_enc_diff_file)
