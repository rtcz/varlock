#!/usr/bin/env python

import base64
import hashlib
import os

from Crypto import Random
from Crypto.Cipher import AES


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
        self.key = hashlib.sha256(key.encode()).digest()
    
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
