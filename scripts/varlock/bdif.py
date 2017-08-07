import json
import os
import struct
import io

import math

from varlock.common import bytes2seq, seq2bytes, dict2bytes, bytes2dict, file_size


class BdiffFile:
    """
    Diff is binary file, consisting of two blocks which are handled separately.
    The first block contains SNV records, the second block INDEL records.
    Each block uses temporary file while writing.
    These two temporary files are merged when writing is finished.

    Diff header:
    header_size, 4B
    header, <header_size>B
    number of SNV records, 4B
    number of INDEL records, 4B
    
    # md5 checksum of mutated BAM file, 16B
    # unmapped secret, 16B
    # start, first genomic index of DIFF range, 4B
    # end, last genomic index of DIFF range, 4B

    Diff SNV record: permutation of bases (A,T,G,C)
    index 4B, 0-based absolute genomic position
    mapping, 1B, index of base permutation

    DIFF INDEL record: permutation of sequences
    index 4B, 0-based absolute genomic position
    length 1B, index of INDEL value
    list:, sorted by sequence (ascending)
        INDEL allele count 2B
        INDEL base length 2B, length of INDEL sequence
        list:
            INDEL 4 base sequence, 1B
    """
    
    CHECKSUM_SIZE = 16
    SECRET_SIZE = 16
    INT_SIZE = 4
    SHORT_SIZE = 2
    SNV_RECORD_SIZE = 5
    
    INDEX_2_PERM = [
        ('A', 'T', 'G', 'C'),
        ('A', 'T', 'C', 'G'),
        ('A', 'C', 'T', 'G'),
        ('A', 'C', 'G', 'T'),
        ('A', 'G', 'T', 'C'),
        ('A', 'G', 'C', 'T'),
        ('T', 'A', 'C', 'G'),
        ('T', 'A', 'G', 'C'),
        ('T', 'C', 'A', 'G'),
        ('T', 'C', 'G', 'A'),
        ('T', 'G', 'A', 'C'),
        ('T', 'G', 'C', 'A'),
        ('C', 'A', 'T', 'G'),
        ('C', 'A', 'G', 'T'),
        ('C', 'T', 'A', 'G'),
        ('C', 'T', 'G', 'A'),
        ('C', 'G', 'A', 'T'),
        ('C', 'G', 'T', 'A'),
        ('G', 'A', 'T', 'C'),
        ('G', 'A', 'C', 'T'),
        ('G', 'T', 'A', 'C'),
        ('G', 'T', 'C', 'A'),
        ('G', 'C', 'A', 'T'),
        ('G', 'C', 'T', 'A')
    ]
    
    PERM_2_INDEX = {
        ('A', 'T', 'G', 'C'): 0,
        ('A', 'T', 'C', 'G'): 1,
        ('A', 'C', 'T', 'G'): 2,
        ('A', 'C', 'G', 'T'): 3,
        ('A', 'G', 'T', 'C'): 4,
        ('A', 'G', 'C', 'T'): 5,
        ('T', 'A', 'C', 'G'): 6,
        ('T', 'A', 'G', 'C'): 7,
        ('T', 'C', 'A', 'G'): 8,
        ('T', 'C', 'G', 'A'): 9,
        ('T', 'G', 'A', 'C'): 10,
        ('T', 'G', 'C', 'A'): 11,
        ('C', 'A', 'T', 'G'): 12,
        ('C', 'A', 'G', 'T'): 13,
        ('C', 'T', 'A', 'G'): 14,
        ('C', 'T', 'G', 'A'): 15,
        ('C', 'G', 'A', 'T'): 16,
        ('C', 'G', 'T', 'A'): 17,
        ('G', 'A', 'T', 'C'): 18,
        ('G', 'A', 'C', 'T'): 19,
        ('G', 'T', 'A', 'C'): 20,
        ('G', 'T', 'C', 'A'): 21,
        ('G', 'C', 'A', 'T'): 22,
        ('G', 'C', 'T', 'A'): 23,
    }
    
    # CHECKSUM = 'checksum'
    # START_INDEX = 'start_index'
    # END_INDEX = 'end_index'
    # SECRET = 'secret'
    
    HEADER_TEMP_EXT = '.header.temp'
    SNV_TEMP_EXT = '.snv.temp'
    INDEL_TEMP_EXT = '.indel.temp'
    
    def __init__(self, filename: str, mode: str, file_obj=None, header: dict = None):
        """
        :param filename:
        :param mode: 'r' for read, 'w' for write
        :param header: used when in write mode
        """
        assert mode in ['r', 'w']
        
        # import gzip
        # gzip.open()
        
        self._tell = 0
        self._mode = mode
        self._filename = filename
        if self._mode == 'r':
            self._init_read()
        elif self._mode == 'w':
            self._header = {} if header is None else header
            self._init_write()
    
    @staticmethod
    def open(header):
        pass
    
    @property
    def mode(self):
        return self._mode
    
    @property
    def header(self):
        """
        :return: header dict
        """
        return self._header
    
    @property
    def filename(self):
        return self._filename
    
    @property
    def indel_count(self):
        return self._indel_count
    
    @property
    def snv_count(self):
        return self._snv_count
    
    # def tell(self):
    #     if self._mode == 'w':
    #         return 0
    #     elif self._mode == 'r':
    #         pass
    
    # def _validate_header(self, header: dict):
    #     assert header is not None
    #     assert header[self.CHECKSUM_SIZE] is bytes \
    #            and len(header[self.CHECKSUM_SIZE]) == self.CHECKSUM_SIZE
    #     assert header[self.SECRET] is bytes \
    #            and len(header[self.SECRET]) == self.SECRET_SIZE
    #     assert header[self.START_INDEX] is int
    #     assert header[self.END_INDEX] is int
    
    def _snv_file_offset(self):
        if self._mode == 'w':
            return 0
        elif self._mode == 'r':
            return self.INT_SIZE + self._header_size
    
    def _indel_file_offset(self):
        if self._mode == 'w':
            return 0
        elif self._mode == 'r':
            return self.INT_SIZE + self._header_size + self._snv_count * self.SNV_RECORD_SIZE
    
    def _init_read(self):
        """
        Create file object for SNV and INDEL file block and seek its start.
        """
        self._snv_file = open(self._filename, 'rb')
        self._indel_file = open(self._filename, 'rb')
        
        # SNV file object is also used to read DIFF header
        self._header_size = struct.unpack('<I', self._snv_file.read(self.INT_SIZE))[0]
        self._header = bytes2dict(self._snv_file.read(self._header_size))
        self._snv_count, self._indel_count = struct.unpack('<II', self._snv_file.read(self.INT_SIZE * 4))
        # counters for tracking a number of read calls
        self._snv_counter = self._indel_counter = 0
        
        # seek start positions
        self._snv_file.seek(self._snv_file_offset())
        self._indel_file.seek(self._indel_file_offset())
    
    def _init_write(self):
        self._snv_file = open(self._filename + self.SNV_TEMP_EXT, 'wb')
        self._indel_file = open(self._filename + self.INDEL_TEMP_EXT, 'wb')
        self._snv_count = self._indel_count = 0
        # counters are not used in write mode
        self._snv_counter = self._indel_counter = None
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def read_snv(self):
        """
        :return: index, (base_1, base_2, base_3, base_4)
        """
        if self._snv_counter == self._snv_count:
            raise EOFError('end of SNVs')
        
        byte_string = self._snv_file.read(self.SNV_RECORD_SIZE)
        self._tell = self._snv_file.tell()
        index, mut_index = struct.unpack('<IB', byte_string)
        
        self._snv_counter += 1
        return index, self.INDEX_2_PERM[mut_index]
    
    def write_snv(self, index: int, perm_tuple: tuple):
        byte_string = struct.pack('<IB', index, self.PERM_2_INDEX[perm_tuple])
        
        self._snv_file.write(byte_string)
        self._snv_count += 1
    
    def read_indel(self):
        """
        :return: index, [SEQ_1, SEQ_2, ...]
        """
        if self._indel_counter == self._indel_count:
            raise EOFError('end of INDELs')
        
        byte_str = self._snv_file.read(self.INT_SIZE * 2)
        index, length = struct.unpack('<IB', byte_str)
        
        alt_list = []
        for i in range(length):
            base_length = struct.unpack('<H', self._snv_file.read(self.SHORT_SIZE))[0]
            seq_byte_size = math.ceil(base_length / 4)
            seq = bytes2seq(self._indel_file.read(seq_byte_size), base_length)
            alt_list.append(seq)
        
        self._indel_counter += 1
        return index, alt_list
    
    def write_indel(self, index: int, alt_list: list):
        assert len(alt_list) < 256
        
        record = struct.pack('<IB', index, len(alt_list))
        for alt in alt_list:
            record += seq2bytes(alt)
        
        self._indel_file.write(record)
        self._indel_count += 1
    
    def close(self):
        self._snv_file.close()
        self._indel_file.close()
        
        if self._mode == 'w':
            with open(self._filename, 'wb') as diff_file, \
                    open(self._filename + self.SNV_TEMP_EXT, 'rb') as snv_file, \
                    open(self._filename + self.INDEL_TEMP_EXT, 'rb') as indel_file:
                header_bytes = dict2bytes(self._header)
                diff_file.write(struct.pack('<I', len(header_bytes)))
                diff_file.write(header_bytes)
                diff_file.write(struct.pack('<II', self._snv_counter, self._indel_counter))
                diff_file.write(snv_file.read())
                diff_file.write(indel_file.read())
            
            os.remove(self._filename + self.HEADER_TEMP_EXT)
            os.remove(self._filename + self.SNV_TEMP_EXT)
            os.remove(self._filename + self.INDEL_TEMP_EXT)
    
    @staticmethod
    def from_text(text_filename: str, bdiff_filename: str):
        with open(text_filename, 'rt') as text_file:
            header = json.loads(text_file.readline().rstrip())
            snv_count = int(text_file.readline().rstrip())
            indel_count = int(text_file.readline().rstrip())
            with BdiffFile(bdiff_filename, 'wb', header) as bdiff_file:
                
                for i in range(snv_count):
                    index, permutation = text_file.readline().rstrip().split('\t')
                    bdiff_file.write_snv(index, permutation.split(','))
                
                for i in range(indel_count):
                    index, alternatives = text_file.readline().rstrip().split('\t')
                    bdiff_file.write_indel(index, alternatives.split(','))
    
    @staticmethod
    def to_text(bdiff_filename: str, text_filename: str):
        with BdiffFile(bdiff_filename, 'rb') as bdiff_file, \
                open(text_filename, 'wt') as text_file:
            header = json.dumps(bdiff_file.header)
            
            text_file.write('%s\n' % header)
            text_file.write('%d\n' % bdiff_file.snv_count)
            text_file.write('%d\n' % bdiff_file.indel_count)
            for i in range(bdiff_file.snv_count):
                index, perm_tuple = bdiff_file.read_snv()
                text_file.write('%d\t%s' % (index, ','.join(perm_tuple)))
            
            for i in range(bdiff_file.indel_count):
                index, alt_list = bdiff_file.read_indel()
                text_file.write('%d\t%s' % (index, ','.join(alt_list)))
    
    # def _tell_snv_index(self, index: int, up: bool):
    #     self._snv_file.seek(self._snv_file_offset())
    #     start_index = self.read_snv()[0]
    #
    #     self._snv_file.seek(self._snv_file_offset() + (self.snv_count - 1) * self.SNV_RECORD_SIZE)
    #     end_index = self.read_snv()[0]
    #
    #     if index > end_index:
    #         # index is after DIFF content
    #         if up:
    #             raise IndexError("Upper index not found")
    #         else:
    #             diff_file.seek(-self.SNV_RECORD_SIZE, os.SEEK_END)
    #             return diff_file.tell()
    #
    #     if index < start_index:
    #         # index is before DIFF content
    #         if up:
    #             diff_file.seek(.HEADER_SIZE)
    #             return diff_file.tell()
    #         else:
    #             raise IndexError("Lower index not found")
    #
    #     first = 0
    #     last = cls.body_size(diff_file) / cls.SNV_RECORD_SIZE - 1
    #     curr_index = index
    #     while first <= last:
    #         mid = int((first + last) / 2)
    #         diff_file.seek(cls.HEADER_SIZE + (mid * cls.SNV_RECORD_SIZE))
    #         curr_index = cls.__read_index(diff_file)
    #         if index == curr_index:
    #             break
    #         else:
    #             if index < curr_index:
    #                 last = mid - 1
    #             else:
    #                 first = mid + 1
    #
    #     diff_file.seek(-self.INT_SIZE, os.SEEK_CUR)
    #     offset = diff_file.tell()
    #
    #     # final move
    #     if index < curr_index:
    #         diff_file.seek(-self.SNV_RECORD_SIZE, os.SEEK_CUR)
    #         prev_offset = diff_file.tell()
    #         prev_index = cls.__read_index(diff_file)
    #         if not up or index == prev_index:
    #             offset = prev_offset
    #
    #     elif index > curr_index:
    #         diff_file.seek(self.SNV_RECORD_SIZE, os.SEEK_CUR)
    #         next_offset = diff_file.tell()
    #         next_index = cls.__read_index(diff_file)
    #         if up or index == next_index:
    #             offset = next_offset
    #
    #     diff_file.seek(offset)
    #
    #     # TODO search INDEL candidates
    #
    #     return offset, curr_index
    
    # def _tell_indel_index(self, index, up):
    #     self._indel_file.seek(self._indel_file_offset())
    #     # indel_index = self.read_indel()[0]
    #     # indel_offset = self._indel_file.tell()
    #     for i in range(self._indel_count):
    #         indel_index = self.read_indel()[0]
    #         if up:
    #             if index >=
    #         else:
    #             pass
            
    
    def tell_index(self, index: int, up=True):
        """
        Only write mode is supported.
        :param index: initial index
        :param up: direction of seeking
        :return: file position of closest index
        """
        # cls.validate(diff_file)
        
        # TODO check min length
        
        prev_snv_pos = self._snv_file.tell()
        prev_indel_pos = self._indel_file.tell()
        
        snv_index = None
        snv_index_pos = None
        offset = None
        
        if self._snv_count == self._indel_count == 0:
            return None
        
        if self._indel_count == 0:
            return self._tell_snv_index(index, up)
        elif self._snv_count == 0:
            return self._tell_indel_index(index, up)
        else:
            snv_offset, snv_index = self._tell_snv_index(index, up)
            indel_offset, indel_index = self._tell_indel_index(index, up)
            
            if abs(index - snv_index) < abs(index - indel_index):
                offset = snv_offset
            else:
                offset = indel_offset
        
        self._snv_file.seek(prev_snv_pos)
        self._indel_file.seek(prev_indel_pos)

