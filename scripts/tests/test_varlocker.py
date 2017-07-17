import unittest

from Crypto.PublicKey import RSA

from varlock.varlocker import Varlocker


class TestVarlocker(unittest.TestCase):
    RESOURCE_PATH = 'tests/resources/varlocker/'
    KEY_PASS = 'password'
    
    @classmethod
    def setUpClass(cls):
        cls.locker = Varlocker()
    
    def test_encrypt(self):
        with open(self.RESOURCE_PATH + 'admin', 'r') as key_file, \
                open(self.RESOURCE_PATH + 'admin.pub', 'r') as pub_key_file:
            rsa_key = RSA.importKey(key_file.read(), passphrase=self.KEY_PASS)
            rsa_pub_key = RSA.importKey(pub_key_file.read())
            self.locker.encrypt(
                rsa_sign_key=rsa_key,
                rsa_enc_key=rsa_pub_key,
                bam_filename=self.RESOURCE_PATH + 'encrypt/input.bam',
                vac_filename=self.RESOURCE_PATH + 'encrypt/input.vac',
                out_bam_filename=self.RESOURCE_PATH + 'encrypt/output.mut.bam',
                out_enc_diff_filename=self.RESOURCE_PATH + 'encrypt/output.diff.enc'
            )
    
    def test_reencrypt(self):
        with open(self.RESOURCE_PATH + 'admin', 'r') as key_file, \
                open(self.RESOURCE_PATH + 'user.pub', 'r') as pub_key_file, \
                open(self.RESOURCE_PATH + 'admin.pub', 'r') as ver_key_file:
            rsa_key = RSA.importKey(key_file.read(), passphrase=self.KEY_PASS)
            rsa_enc_key = RSA.importKey(pub_key_file.read())
            rsa_ver_key = RSA.importKey(ver_key_file.read())
            
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
    
    def test_decrypt(self):
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
    
    def test_vac(self):
        # tested in vac_test
        pass
