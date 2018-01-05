import os

from Crypto.Cipher import AES

class FileAES:
    """
    AES encryption for binary input
    """
    
    def __init__(self, key: bytes):
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
                chunk += padding_length * bytes([padding_length])
                out_file.write(cipher.encrypt(chunk))
                break
            else:
                out_file.write(cipher.encrypt(chunk))
    
    def decrypt(self, in_file, out_file):
        iv = in_file.read(AES.block_size)
        cipher = AES.new(self.key, AES.MODE_CBC, iv)
        next_chunk = b''
        while True:
            chunk = next_chunk
            next_chunk = cipher.decrypt(in_file.read(1024 * AES.block_size))
            
            if len(next_chunk) == 0:
                padding_length = chunk[-1]
                if padding_length < 1 or padding_length > AES.block_size:
                    raise ValueError("bad decrypt pad (%d)" % padding_length)
                
                # all the pad-bytes must be the same
                if chunk[-padding_length:] != (padding_length * bytes([padding_length])):
                    # this is similar to the bad decrypt:evp_enc.c from openssl program
                    raise ValueError("bad decrypt")
                
                chunk = chunk[:-padding_length]
                out_file.write(chunk)
                break
            else:
                out_file.write(chunk)
