import unittest

import pysam
from Crypto.PublicKey import RSA

import src.common as cmn
from src.vac import Vac
from src.varlocker import Varlocker


class TestVarlocker(unittest.TestCase):
    """
    Tests are order sensitive.
    """
    # TODO make tests isolated
    
    RESOURCE_PATH = 'tests/resources/varlocker/'
    KEY_PASS = 'password'
    
    @classmethod
    def setUpClass(cls):
        cls.locker = Varlocker()
    
    def test1_encrypt(self):
        cmn.sam2bam(self.RESOURCE_PATH + 'encrypt/input.sam', self.RESOURCE_PATH + 'encrypt/input.bam')
        pysam.index(self.RESOURCE_PATH + 'encrypt/input.bam')
        Vac.text2vac(self.RESOURCE_PATH + 'encrypt/input.vac.txt', self.RESOURCE_PATH + 'encrypt/input.vac')
        
        with open(self.RESOURCE_PATH + 'admin', 'r') as key_file, \
                open(self.RESOURCE_PATH + 'admin.pub', 'r') as pub_key_file:
            rsa_key = RSA.importKey(key_file.read(), passphrase=self.KEY_PASS)
            rsa_pub_key = RSA.importKey(pub_key_file.read())
            
            # creates DIFF with secret
            self.locker.encrypt(
                rsa_sign_key=rsa_key,
                rsa_enc_key=rsa_pub_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/input.bam',
                vac_filename=self.RESOURCE_PATH + 'encrypt/input.vac',
                out_bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                out_enc_diff_filename=self.RESOURCE_PATH + 'encrypt/output.diff.enc',
                mut_p=0
            )
        
        pysam.index(self.RESOURCE_PATH + 'encrypt/output.mut.bam')
        cmn.bam2sam(self.RESOURCE_PATH + 'encrypt/output.mut.bam', self.RESOURCE_PATH + 'encrypt/output.mut.sam')
    
    def test2_reencrypt(self):
        with open(self.RESOURCE_PATH + 'admin', 'r') as key_file, \
                open(self.RESOURCE_PATH + 'user.pub', 'r') as pub_key_file, \
                open(self.RESOURCE_PATH + 'admin.pub', 'r') as ver_key_file:
            rsa_key = RSA.importKey(key_file.read(), passphrase=self.KEY_PASS)
            rsa_enc_key = RSA.importKey(pub_key_file.read())
            rsa_ver_key = RSA.importKey(ver_key_file.read())
            
            # creates DIFF without secret
            self.locker.reencrypt(
                rsa_key=rsa_key,
                rsa_enc_key=rsa_enc_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                enc_diff_filename=self.RESOURCE_PATH + 'encrypt/output.diff.enc',
                out_enc_diff_filename=self.RESOURCE_PATH + 'reencrypt/output.diff.enc',
                include_unmapped=False,
                unmapped_only=False,
                rsa_ver_key=rsa_ver_key
            )
    
    def test3_decrypt(self):
        with open(self.RESOURCE_PATH + 'user', 'r') as key_file, \
                open(self.RESOURCE_PATH + 'admin.pub', 'r') as ver_key_file:
            rsa_key = RSA.importKey(key_file.read(), passphrase=self.KEY_PASS)
            rsa_ver_key = RSA.importKey(ver_key_file.read())
            
            self.locker.decrypt(
                rsa_key=rsa_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                enc_diff_filename=self.RESOURCE_PATH + 'reencrypt/output.diff.enc',
                out_bam_filename=self.RESOURCE_PATH + 'decrypt/output.bam',
                include_unmapped=False,
                unmapped_only=False,
                rsa_ver_key=rsa_ver_key,
            )
    
    def test4_reencrypt_unmapped(self):
        with open(self.RESOURCE_PATH + 'admin', 'r') as key_file, \
                open(self.RESOURCE_PATH + 'user.pub', 'r') as pub_key_file, \
                open(self.RESOURCE_PATH + 'admin.pub', 'r') as ver_key_file:
            rsa_key = RSA.importKey(key_file.read(), passphrase=self.KEY_PASS)
            rsa_enc_key = RSA.importKey(pub_key_file.read())
            rsa_ver_key = RSA.importKey(ver_key_file.read())
            
            # missing DIFF secret
            self.assertRaises(ValueError, lambda: self.locker.reencrypt(
                rsa_key=rsa_key,
                rsa_enc_key=rsa_enc_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                enc_diff_filename=self.RESOURCE_PATH + 'reencrypt/output.diff.enc',
                out_enc_diff_filename=self.RESOURCE_PATH + 'reencrypt/output2.diff.enc',
                include_unmapped=False,
                unmapped_only=True,
                rsa_ver_key=rsa_ver_key
            ))
            
            self.locker.reencrypt(
                rsa_key=rsa_key,
                rsa_enc_key=rsa_enc_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                enc_diff_filename=self.RESOURCE_PATH + 'encrypt/output.diff.enc',
                out_enc_diff_filename=self.RESOURCE_PATH + 'reencrypt/output2.diff.enc',
                include_unmapped=False,
                unmapped_only=True,
                rsa_ver_key=rsa_ver_key
            )
    
    def test5_decrypt_unmapped(self):
        with open(self.RESOURCE_PATH + 'admin', 'r') as key_file, \
                open(self.RESOURCE_PATH + 'admin.pub', 'r') as ver_key_file:
            rsa_key = RSA.importKey(key_file.read(), passphrase=self.KEY_PASS)
            rsa_ver_key = RSA.importKey(ver_key_file.read())
            
            # missing DIFF secret
            self.assertRaises(ValueError, lambda: self.locker.decrypt(
                rsa_key=rsa_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                enc_diff_filename=self.RESOURCE_PATH + 'reencrypt/output.diff.enc',
                out_bam_filename=self.RESOURCE_PATH + 'decrypt/output2.bam',
                include_unmapped=False,
                unmapped_only=True,
                rsa_ver_key=rsa_ver_key,
            ))
            
            self.locker.decrypt(
                rsa_key=rsa_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                enc_diff_filename=self.RESOURCE_PATH + 'encrypt/output.diff.enc',
                out_bam_filename=self.RESOURCE_PATH + 'decrypt/output2.bam',
                include_unmapped=False,
                unmapped_only=True,
                rsa_ver_key=rsa_ver_key,
            )
