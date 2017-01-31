import os
import sqlite3
from io import BytesIO

from Crypto.PublicKey import RSA

import varlock


class Varlock:
    RSA_KEY_LENGTH = 4096
    AES_KEY_LENGTH = 32
    
    PUBLIC_KEY_SUFFIX = '.pub'
    AES_KEY_SUFFIX = '.key.enc'
    
    DB_FILENAME = 'varlock.db'
    
    def __init__(self, repo_dir, fai_filepath):
        self.repo_dir = repo_dir
        self.fai_filepath = fai_filepath
        
        db_filepath = os.path.join(self.repo_dir, self.DB_FILENAME)
        self.db = sqlite3.connect(db_filepath)
    
    def insert_user(self, username, rsa_key, password):
        private_key = rsa_key.exportKey(format='PEM', passphrase=password)
        public_key = rsa_key.publickey().exportKey()
        
        self.db.execute("""
            INSERT INTO users(username, private_key, public_key) VALUES (?, ?, ?)
        """, (username, private_key, public_key))
    
    def delete_user(self, username):
        self.db.execute("""
            DELETE FROM users WHERE username = ?
        """, username)
    
    
    
    
    def __insert_user(self, key, username, password):
        """
        :param key:  An RSA key object (`RsaKey`)
        :param username:
        :param password:
        """
        byte_string = key.exportKey(format='PEM', passphrase=password)
        self.db.execute('INSERT INTO users (private_key) VALUES (?)', byte_string)
    
    def __insert_public_key(self, key, username):
        """
        :param key:  An RSA key object (`RsaKey`)
        :param username:
        """
        byte_string = key.publickey().exportKey()
        self.db.execute('INSERT INTO users (public_key) VALUES (?)', byte_string)
    
    def __select_private_key(self, username, password):
        """
        :param username:
        :param password:
        :return: An RSA key object (`RsaKey`)
        """
        key_filepath = os.path.join(self.key_dirpath, username)
        if os.path.isfile(key_filepath):
            raise FileNotFoundError(key_filepath)
        
        with open(key_filepath, 'rb') as key_file:
            return RSA.import_key(key_file.readall(), password)
    
    def __select_public_key(self, username):
        """
        :param username:
        :return: An RSA key object (`RsaKey`)
        """
        key_filepath = os.path.join(self.key_dirpath, username + self.PUBLIC_KEY_SUFFIX)
        if os.path.isfile(key_filepath):
            raise FileExistsError(key_filepath)
        
        with open(key_filepath, 'rb') as key_file:
            return RSA.import_key(key_file.readall())
    
    def __insert_aes_key(self, aes_key, username, sample_name):
        key_filepath = os.path.join(self.key_dirpath, self.aes_key_filepath(username, sample_name))
        if os.path.isfile(key_filepath):
            raise FileExistsError(key_filepath)
        
        with open(key_filepath, 'wb') as key_file:
            key_file.write(aes_key)
    
    def __select_aes_key(self, username, sample_name):
        key_filepath = os.path.join(self.key_dirpath, self.aes_key_filepath(username, sample_name))
        if not os.path.isfile(key_filepath):
            raise FileNotFoundError(key_filepath)
        
        with open(key_filepath, 'rb') as key_file:
            return key_file.read()
    
    def aes_key_filepath(self, username, sample_name):
        return username + '_' + sample_name + self.AES_KEY_SUFFIX
    
    def add_user(self, username, password):
        """
        Create asymetric key pair. Private key is encrypted by provided password.
        :return:
        """
        key = RSA.generate(self.RSA_KEY_LENGTH)
        self.__insert_private_key(key, username, password)
        self.__insert_public_key(key, username)
    
    def remove_user(self, username):
        """
        Remove asymetric key pair.
        """
        # TODO delete also user's aes key files (username_samplename.enc.key)
        private_key_filepath = os.path.join(self.key_dirpath, username)
        public_key_filepath = private_key_filepath + self.PUBLIC_KEY_SUFFIX
        if os.path.isfile(private_key_filepath):
            os.remove(os.path.join(self.key_dirpath, username))
        
        if os.path.isfile(public_key_filepath):
            os.remove(os.path.join(self.key_dirpath, username + self.PUBLIC_KEY_SUFFIX))
    
    def update_password(self, username, old_password, new_password):
        """
        Decrypt asymetric key pair by old password and encrypt it with new password.
        :param username:
        :param old_password:
        :param new_password:
        :return:
        """
        private_key_filepath = os.path.join(self.key_dirpath, username)
        if not os.path.isfile(private_key_filepath):
            raise IOError("User %s does not exists" % username)
        
        with open(private_key_filepath, 'rb') as in_file:
            try:
                key = RSA.import_key(in_file.read(), old_password)
                with open(private_key_filepath, 'wb') as out_file:
                    key_text = key.exportKey(format='PEM', passphrase=new_password)
                    out_file.write(key_text)
            except ValueError:
                raise ValueError('Incorrect password')
    
    # @staticmethod
    # def encrypt_file(raw_filepath, enc_filepath, key):
    #     """
    #     :param raw_filepath: binary file
    #     :param enc_filepath: binary file
    #     :param key:
    #     """
    #     aes = varlock.FileAES(key)
    #     with open(raw_filepath, 'rb') as raw_file, \
    #             open(enc_filepath, 'wb') as enc_file:
    #         aes.encrypt(raw_file, enc_file)
    
    # @staticmethod
    # def decrypt_file(enc_filepath, raw_filepath, key):
    #     aes = varlock.FileAES(key)
    #     with open(enc_filepath, 'rb') as enc_file, \
    #             open(raw_filepath, 'wb') as raw_file:
    #         aes.decrypt(enc_file, raw_file)
    
    # @staticmethod
    # def get_decrypted_file(enc_filepath, key):
    #     raw_file = BytesIO()
    #     aes = varlock.FileAES(key)
    #     with open(enc_filepath, 'rb') as enc_file:
    #         aes.decrypt(enc_file, raw_file)
    #
    #     return raw_file
    
    def mutate(self, username, sample_name, in_vac_filepath, in_bam_filepath):
        """
        Creates files:
        sample_name.mut.bam, sample_name.diff.enc, sample_name.key.enc
        Any existing user can mutate sample.
        :param username: diff is encrypted with public aes_key which belongs to this user
        :param sample_name:
        :param in_vac_filepath:
        :param in_bam_filepath:
        """
        mut = varlock.Mutator()
        
        # try to import public key before mutating
        # TODO check if users exists
        rsa_key = self.__select_public_key(username)
        aes_key = self.gen_password()
        
        if not os.path.isfile(in_vac_filepath):
            raise FileNotFoundError("File %s not found" % in_vac_filepath)
        
        if not os.path.isfile(in_bam_filepath):
            raise FileNotFoundError("File %s not found" % in_bam_filepath)
        
        # generate mutated bam and diff file
        diff_file = BytesIO()
        
        # mutate and create diff
        out_bam_filepath = os.path.join(self.sample_dirpath, sample_name + '.mut.bam')
        with pysam.AlignmentFile(in_bam_filepath, "rb") as in_bam_file, \
                open(in_vac_filepath, "rb") as in_vac_file, \
                pysam.AlignmentFile(out_bam_filepath, "wb", template=in_bam_file) as out_bam_file:
            mut.mutate(in_vac_file, in_bam_file, out_bam_file, diff_file)
        
        # encrypt diff with new AES key
        enc_diff_filepath = os.path.join(self.sample_dirpath, sample_name + '.diff.enc')
        with open(enc_diff_filepath, 'rb') as enc_diff_file:
            aes = varlock.FileAES(aes_key)
            aes.encrypt(diff_file, enc_diff_file)
        
        # encrypt AES key with public key
        enc_aes_key = rsa_key.encrypt(aes_key)
        self.__insert_aes_key(enc_aes_key, username, sample_name)
    
    def has_access(self, username, sample_name, ref_name=None, start_pos=None, end_pos=None):
        # TODO
        key_filepath = os.path.join(self.key_dirpath, self.aes_key_filepath(username, sample_name))
        return os.path.isfile(key_filepath)
    
    def add_access(self, username, password, to_username, ref_name=None, start_pos=None, end_pos=None):
        pass
    
    def remove_access(self, username, password, to_username, ref_name=None, start_pos=None, end_pos=None):
        pass
    
    def unmutate(self, username, password, sample_name, ref_name, start_pos, end_pos, out_bam_filepath):
        """
        :param username:
        :param password:
        :param sample_name:
        :param ref_name:
        :param start_pos:
        :param end_pos:
        :param out_bam_filepath:
        :return:
        """
        if not self.has_access(username, sample_name):
            raise PermissionError("User %s can not unmutate sample %s" % (username, sample_name))
        
        # try to import private key before unmutating
        # TODO check if users exists
        rsa_key = self.__select_private_key(username, password)
        aes_key = rsa_key.decrypt(self.__select_aes_key(username, sample_name))
        
        # decrypt diff
        diff_enc_filepath = os.path.join(self.sample_dirpath, sample_name + '.diff.enc')
        diff_file = BytesIO()
        aes = varlock.FileAES(aes_key)
        with open(diff_enc_filepath, 'rb') as diff_enc_file:
            aes.decrypt(diff_enc_file, diff_file)
        
        # unmutate
        in_bam_filepath = os.path.join(self.sample_dirpath, sample_name + '.mut.bam')
        mut = varlock.Mutator()
        with pysam.AlignmentFile(in_bam_filepath, "rb") as in_bam_file, \
                pysam.AlignmentFile(out_bam_filepath, "wb") as out_bam_file:
            mut.unmutate(in_bam_file, diff_file, ref_name, start_pos, end_pos, out_bam_file)
    
    @classmethod
    def gen_password(cls):
        return os.urandom(cls.AES_KEY_LENGTH)


if __name__ == "__main__":
    RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')
    
    vrl = Varlock(RESOURCES_DIR)
    vrl.remove_user('jozko')
    vrl.add_user('jozko', 'jozko123')
    vrl.update_password('jozko', 'jozko123', 'jozko321')
    vrl.update_password('jozko', 'jozko321', 'jozko123')
    # vrl.unmutate()
