import io
import os
import struct

import pysam
from Crypto.Cipher import PKCS1_OAEP
from Crypto.PublicKey import RSA

import varlock as vrl
from varlock import Diff


class Varlocker:
    RSA_KEY_LENGTH = 4096
    RSA_CIPHER_LENGTH = 512
    AES_KEY_LENGTH = 32
    
    @staticmethod
    def create_vac(bam_filename: str, vcf_filename: str, out_vac_filename: str):
        """
        :param bam_filename:
        :param vcf_filename:
        :param out_vac_filename:
        :return:
        """
        with pysam.AlignmentFile(bam_filename, 'rb') as sam_file, \
                open(vcf_filename, 'rb') as vcf_file, \
                open(out_vac_filename, 'wb') as out_vac_file:
            vac = vrl.Vac(vrl.FastaIndex(sam_file))
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
        
        with mut.mutate(bam_filename, vac_filename, out_bam_filename) as diff_file, \
                open(out_enc_diff_filename, 'wb') as enc_diff_file:
            # store encrypted AES key in encrypted DIFF
            enc_aes_key = PKCS1_OAEP.new(rsa_key.publickey()).encrypt(aes_key)
            enc_diff_file.write(struct.pack('<I', len(enc_aes_key)))
            enc_diff_file.write(enc_aes_key)
            
            if verbose:
                print('Encrypting DIFF')
            
            # encrypt diff with new AES key
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
            
            if verbose:
                print('Unmutating BAM')
            
            # unmutate
            mut = vrl.BamMutator(verbose)
            mut.unmutate(
                bam_filename,
                diff_file,
                out_bam_filename,
                start_ref_name,
                start_ref_pos,
                end_ref_name,
                end_ref_pos
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
        :param verbose:
        """
        with pysam.AlignmentFile(bam_filename, 'rb') as sam_file, \
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
            # TODO refactor
            Diff.validate(diff_file)
            mut = vrl.Mutator(sam_file)
            start_index, end_index = mut.resolve_diff_range(
                diff_file,
                start_ref_name,
                start_ref_pos,
                end_ref_name,
                end_ref_pos
            )
            # slice diff
            with Diff.slice(diff_file, start_index, end_index) as sliced_diff, \
                    open(out_enc_diff_filename, 'rb') as out_enc_diff_file:
                # store encrypted AES key in encrypted DIFF
                enc_aes_key = PKCS1_OAEP.new(rsa_dec_key.publickey()).encrypt(aes_key)
                enc_diff_file.write(struct.pack('<I', len(enc_aes_key)))
                enc_diff_file.write(enc_aes_key)
                
                if verbose:
                    print('Encrypting DIFF')
                
                # encrypt diff with AES key
                out_aes_key = os.urandom(cls.AES_KEY_LENGTH)
                out_aes = vrl.FileAES(out_aes_key)
                out_aes.encrypt(sliced_diff, out_enc_diff_file)
    
    @classmethod
    def has_privilege(
            cls,
            private_key,
            enc_aes_key,
            enc_diff_file,
            mut_bam_file,
            ref_name=None,
            start_pos=None,
            end_pos=None
    ):
        """
        :param private_key:
        :param enc_aes_key:
        :param enc_diff_file:
        :param mut_bam_file:
        :param ref_name:
        :param start_pos:
        :param end_pos:
        :return:
        """
        pass
        # TODO is this method needed ???
        # if ref_name is not None and start_pos is not None and end_pos is not None:
        #     # range is specified
        #     if start_pos >= end_pos:
        #         # wrong range
        #         raise ValueError("End position must be greater than start position.")
        #
        # # TODO 1 try to decrypt aes
        # aes_key = private_key.decrypt(enc_aes_key)
        #
        # # TODO 2 try to decrypt diff, only header part
        #
        # # TODO 3 check the range specified in diff's header
        #
        # # TODO 4 compare checksum of BAM to checksum from diff's header
        #
        # # TODO 5 check if range is present in BAM
        #
        # return True
