#!/usr/bin/env python

import base64
import os
from hashlib import sha256

from Crypto import Random
from Crypto.Cipher import AES


class FileAES:
    """
    AES encryption for binary input
    """
    
    def __init__(self, key):
        """
        :param key: byte string
        The secret key to use in the symmetric cipher.
        It must be 16 (*AES-128*), 24 (*AES-192*), or 32 (*AES-256*) bytes long.
        """
        self.key = key
    
    def encrypt(self, in_file, out_file):
        iv = os.urandom(AES.block_size)
        out_file.write(iv)
        cipher = AES.new(self.key, AES.MODE_CBC, iv)
        while True:
            chunk = in_file.read(1024 * AES.block_size)
            if len(chunk) == 0 or len(chunk) % AES.block_size != 0:
                padding_length = AES.block_size - (len(chunk) % AES.block_size)
                chunk += padding_length * chr(padding_length)
                out_file.write(cipher.encrypt(chunk))
                break
            else:
                out_file.write(cipher.encrypt(chunk))
    
    def decrypt(self, in_file, out_file):
        iv = in_file.read(AES.block_size)
        cipher = AES.new(self.key, AES.MODE_CBC, iv)
        next_chunk = ''
        while True:
            chunk = next_chunk
            next_chunk = cipher.decrypt(in_file.read(1024 * AES.block_size))
            
            if len(next_chunk) == 0:
                padding_length = ord(chunk[-1])
                if padding_length < 1 or padding_length > AES.block_size:
                    raise ValueError("bad decrypt pad (%d)" % padding_length)
                
                # all the pad-bytes must be the same
                if chunk[-padding_length:] != (padding_length * chr(padding_length)):
                    # this is similar to the bad decrypt:evp_enc.c from openssl program
                    raise ValueError("bad decrypt")
                
                chunk = chunk[:-padding_length]
                out_file.write(chunk)
                break
            else:
                out_file.write(chunk)


class BinAES:
    """
    AES encryption for binary input
    """
    
    def __init__(self, key):
        """
        :param key: byte string
        The secret key to use in the symmetric cipher.
        It must be 16 (*AES-128*), 24 (*AES-192*), or 32 (*AES-256*) bytes long.
        """
        self.key = key
    
    def encrypt(self, raw_bytes):
        """
        :param raw_bytes: message bytes, length of message must be multiple of 16 bytes
        :return: cipher bytes
        """
        padded_bytes = self.pad(raw_bytes)
        iv = os.urandom(AES.block_size)
        cipher = AES.new(self.key, AES.MODE_CBC, iv)
        return iv + cipher.encrypt(padded_bytes)
    
    def decrypt(self, enc_bytes):
        """
        :param enc_bytes: cipher bytes
        :return: message bytes
        """
        iv = enc_bytes[:AES.block_size]
        cipher = AES.new(self.key, AES.MODE_CBC, iv)
        return self.unpad(cipher.decrypt(enc_bytes[AES.block_size:]))
    
    def key_hash(self):
        return self.key
    
    @staticmethod
    def pad(value):
        pad_count = AES.block_size - len(value) % AES.block_size
        return value + pad_count * b"\0"
    
    @staticmethod
    def unpad(value):
        return value.rstrip(b"\0")


class TextAES:
    """
    AES encryption for text input
    """
    
    def __init__(self, key):
        """
        transforms key to 32 bytes long hash
        :param key: text key
        """
        # transform key to 32 bytes long hash
        self.key = sha256(key.encode()).digest()
    
    def encrypt(self, raw):
        """
        :param raw: text message
        :return: cipher in base64 format
        """
        raw = self.pad(raw)
        iv = Random.new().read(AES.block_size)
        cipher = AES.new(self.key, AES.MODE_CBC, iv)
        return base64.b64encode(iv + cipher.encrypt(raw))
    
    def decrypt(self, enc):
        """
        :param enc: cipher in base64 format
        :return: text message
        """
        enc = base64.b64decode(enc)
        iv = enc[:16]
        cipher = AES.new(self.key, AES.MODE_CBC, iv)
        return self.unpad(cipher.decrypt(enc[16:]))
    
    def key_hash(self):
        return self.key
    
    @staticmethod
    def pad(value):
        pad_count = AES.block_size - len(value) % AES.block_size
        return value + pad_count * chr(pad_count)
    
    @staticmethod
    def unpad(value):
        return value[0:-ord(value[-1])]
