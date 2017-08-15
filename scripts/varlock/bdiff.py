import json
import struct
import io
import os

import math

from varlock.common import bytes2seq, seq2bytes, dict2bytes, bytes2dict


class BdiffIO:
    """
    Bdiff is binary file containing SNV and INDEL difference records.

    Bdiff header:
    header_size, 4B
    header, <header_size>B
    number of SNV records, 4B
    number of INDEL records, 4B
    
    Bdiff index:
    last_index 4B, index of last record
    index length 4B, number of index records
    index resolution 4B, distance between indexed values as number of records
    list: of index length
        index of a thousandth record, 4B
        file position of a thousandth record, 4B

    Diff SNV record: permutation of bases (A,T,G,C)
    index 4B, 0-based absolute genomic position
    length 1B, always zero to indicate SNV record (SNV has always 4 alternatives)
    mapping, 1B, index of base permutation

    DIFF INDEL record: permutation of sequences
    index 4B, 0-based absolute genomic position
    length 1B, number of alternatives
    list:, of length and sorted by sequence (ascending)
        INDEL allele count 2B
        INDEL base length 2B, length of INDEL sequence
        list:
            INDEL 4 base sequence, 1B
    """
    INT_SIZE = 4
    SHORT_SIZE = 2
    BYTE_SIZE = 1
    
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
    
    def __init__(self, diff_file: io.BytesIO = None, index_resolution=1000):
        """
        Instance is either in read or write mode.
        :param diff_file: provide for read mode
        """
        self._file_index = []
        self._header = {}
        self._snv_count = self._indel_count = 0
        self._last_index = 0
        self._header_size = 0
        self._index_resolution = index_resolution
        self._diff_file = None
        self._record_file = None
        
        if diff_file is not None:
            # read mode
            diff_file.seek(0)
            meta_size = struct.unpack('<I', diff_file.read(self.INT_SIZE))[0]
            self._header = bytes2dict(diff_file.read(meta_size))
            self._snv_count, self._indel_count = struct.unpack('<II', diff_file.read(self.INT_SIZE * 2))
            
            # read_record index
            self._last_index, \
            index_length, \
            self._index_resolution = struct.unpack('<III', diff_file.read(self.INT_SIZE * 3))
            
            for i in range(index_length):
                index, pos = struct.unpack('<II', diff_file.read(self.INT_SIZE * 2))
                self._file_index.append((index, pos))
            
            self._header_size = diff_file.tell()
            self._diff_file = diff_file
        else:
            # write mode
            self._record_file = io.BytesIO()
    
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            return self.read_record()
        except EOFError:
            raise StopIteration
    
    def is_empty(self):
        return self._snv_count + self.indel_count == 0
    
    @property
    def file_index(self):
        return self._file_index.copy()
    
    @property
    def index_resolution(self):
        return self._index_resolution
    
    @property
    def is_read_mode(self):
        return self._diff_file is not None
    
    @property
    def header_size(self):
        """
        :return: header size in bytes
        """
        return self._header_size
    
    @property
    def header(self):
        """
        :return: meta header dict
        """
        return self._header.copy()
    
    @property
    def indel_count(self):
        return self._indel_count
    
    @property
    def snv_count(self):
        return self._snv_count
    
    @property
    def first_index(self):
        """
        :return: genomic index of first record
        """
        if self.is_empty():
            return None
        
        if self.is_read_mode:
            curr_pos = self._diff_file.tell()
            self._diff_file.seek(self._header_size)
            first_index = struct.unpack('<I', self._diff_file.read(self.INT_SIZE))[0]
            self._diff_file.seek(curr_pos)
        else:
            curr_pos = self._record_file.tell()
            self._record_file.seek(0)
            first_index = struct.unpack('<I', self._record_file.read(self.INT_SIZE))[0]
            self._record_file.seek(curr_pos)
        
        return first_index
    
    @property
    def last_index(self):
        """
        :return: genomic index of last record
        """
        return None if self.is_empty() else self._last_index
    
    def read_record(self):
        """
        Reads next SNV or INDEL record
        :return: index, is_snv, alternatives
        SNV alternatives is a tuple
        INDEL alternatives is a list
        """
        assert self.is_read_mode
        byte_str = self._diff_file.read(self.INT_SIZE + self.BYTE_SIZE)
        if byte_str == b'':
            raise EOFError
        
        index, length = struct.unpack('<IB', byte_str)
        if length == 0:
            # SNV record
            return index, self._read_snv_alts()
        else:
            # INDEL record
            return index, self._read_indel_alts(length)
    
    def _read_snv_alts(self):
        """
        Read SNV alternatives from rest of the record
        :return: (base_1, base_2, base_3, base_4)
        """
        mut_index = struct.unpack('<B', self._diff_file.read(self.BYTE_SIZE))[0]
        return self.INDEX_2_PERM[mut_index]
    
    def _read_indel_alts(self, length):
        """
        Read INDEL alternatives from rest of the record
        :return: [SEQ_1, SEQ_2, ...]
        """
        alt_list = [None] * length
        for i in range(length):
            allele_count, base_length = struct.unpack('<HH', self._diff_file.read(self.SHORT_SIZE * 2))
            seq_byte_size = math.ceil(base_length / 4)
            seq = bytes2seq(self._diff_file.read(seq_byte_size), base_length)
            alt_list[i] = (allele_count, seq)
        
        return alt_list
    
    def write_snv(self, index: int, perm_tuple: tuple):
        assert not self.is_read_mode
        assert index > self._last_index
        
        record = struct.pack('<IBB', index, 0, self.PERM_2_INDEX[perm_tuple])
        self._record_file.write(record)
        
        self._last_index = index
        self._snv_count += 1
        if self._snv_count % self._index_resolution == 0:
            self._file_index.append((index, self._record_file.tell() - len(record)))
    
    def write_indel(self, index: int, alt_list: list):
        assert not self.is_read_mode
        assert index > self._last_index
        assert 1 < len(alt_list) < 256
        
        record = struct.pack('<IB', index, len(alt_list))
        for allele_count, seq in alt_list:
            record += struct.pack('<HH', allele_count, len(seq)) + seq2bytes(seq)
        
        self._record_file.write(record)
        
        self._last_index = index
        self._indel_count += 1
        if self._indel_count % self._index_resolution == 0:
            self._file_index.append((index, self._record_file.tell() - len(record)))
    
    def readcopy(self):
        return BdiffIO(self.file())
    
    def file(self, header: dict = None, close=True):
        """
        Creates Bdiff file using internal data.
        :param header: meta header dict, only in write mode
        :param close: closes internal stream in write mode or input stream in read mode
        :return: Bdiff file as BytesIO
        """
        # TODO is_closed ?
        diff_file = io.BytesIO()
        if self.is_read_mode:
            # copy of file content supplied in constructor
            diff_file.write(self._diff_file.read())
            if close:
                self._diff_file.close()
        else:
            # new stream from underlying record_file
            # meta
            meta_bytes = dict2bytes({} if header is None else header)
            diff_file.write(struct.pack('<I', len(meta_bytes)))
            diff_file.write(meta_bytes)
            
            # counts
            diff_file.write(struct.pack('<II', self._snv_count, self._indel_count))
            
            # index
            diff_file.write(struct.pack(
                '<III',
                self._last_index,
                len(self._file_index),
                self._index_resolution)
            )
            for index, pos in self._file_index:
                diff_file.write(struct.pack('<II', index, pos))
            
            # records
            self._record_file.seek(0)
            diff_file.write(self._record_file.read())
            if close:
                self._record_file.close()
            else:
                self._record_file.seek(0, os.SEEK_END)
        
        diff_file.seek(0)
        return diff_file
    
    def tell_range(self, lower_index, upper_index):
        """
        From file position of first record with index greater than or equal to lower_index
        to file position of last record with index lower than or equal to upper_index.
        :param lower_index: from index inclusive
        :param upper_index: to index inclusive
        :return: range tuple
        :raises: IndexError
        """
        assert self.is_read_mode
        assert lower_index <= upper_index
        
        start_pos = self.tell_index_gte(lower_index)
        end_pos = self.tell_index_lte(upper_index)
        
        if start_pos is None:
            raise IndexError("Empty range")
        
        if end_pos is None:
            raise IndexError("Empty range")
        
        # make end_pos inclusive
        # saved_pos = self._diff_file.tell()
        # self._diff_file.seek(end_pos)
        # self.read_record()
        # end_pos_inc = self._diff_file.tell()
        # self._diff_file.seek(saved_pos)
        
        return start_pos, end_pos
    
    def tell_index_gte(self, pivot_index: int):
        """
        :param pivot_index:
        :return: file position of index greater than or equal to pivot_index
        None if such index does not exists.
        """
        assert self.is_read_mode
        
        if self.last_index is None or pivot_index > self.last_index:
            # file is empty or desired index does not exists
            return None
        else:
            saved_pos = self._diff_file.tell()
            curr_pos = self.__indexed_pos(pivot_index)
            self._diff_file.seek(curr_pos)
            
            curr_index = self.read_record()[0]  # type: int
            while curr_index < pivot_index:
                curr_pos = self._diff_file.tell()
                curr_index = self.read_record()[0]  # type: int
            
            self._diff_file.seek(saved_pos)
            return curr_pos
    
    def tell_index_lte(self, pivot_index: int):
        """
        :param pivot_index:
        :return: file position of index lower than or equal to pivot_index
        None if such index does not exists.
        """
        assert self.is_read_mode
        
        if self.first_index is None or pivot_index < self.first_index:
            # desired index does not exists
            return None
        else:
            saved_pos = self._diff_file.tell()
            curr_pos = self.__indexed_pos(pivot_index)
            self._diff_file.seek(curr_pos)
            
            while True:
                next_pos = self._diff_file.tell()
                try:
                    next_index = self.read_record()[0]  # type: int
                except EOFError:
                    break
                if next_index <= pivot_index:
                    curr_pos = next_pos
                else:
                    break
            
            self._diff_file.seek(saved_pos)
            return curr_pos
    
    def __indexed_pos(self, pivot_index: int):
        """
        :param pivot_index:
        :return: indexed position of record with index lower than or equal to pivot_index.
        Zero position is implicitly indexed and returned if indexed position is not found.
        Returned position if file absolute.
        """
        indexed_pos = 0
        for index, pos in self._file_index:
            if pivot_index <= index:
                indexed_pos = pos
            else:
                break
        
        return indexed_pos + self._header_size


# def to_bytes_io(text_io: io.StringIO):
#     """
#     Convert Bdiff file from StringIO to BytesIO.
#     Indices are converted from 0-based to 1-based.
#     :param text_io:
#     :return:
#     """
#     text_io.seek(0)
#     header = json.loads(text_io.readline().rstrip())
#     # snv count line
#     text_io.readline()
#     # indel count line
#     text_io.readline()
#
#     bdiff = BdiffIO()
#     for line in text_io:
#         raw_index, raw_alts = line.rstrip().split('\t')
#         index = int(raw_index) - 1
#         alts = raw_alts.split(',')
#         assert len(alts) > 0
#         if len(alts[0]) == 1:
#             # is SNV
#             bdiff.write_snv(index, alts)
#         else:
#             # is indel
#             bdiff.write_indel(index, [alt.split(':') for alt in alts])
#
#     return bdiff.file(header)


def to_text_file(bytes_io: io.BytesIO(), filename: str):
    """
    Convert Bdiff file from BytesIO to text file.
    Indices are converted from 0-based to 1-based.
    :param filename:
    :param bytes_io:
    :return:
    """
    bdiff = BdiffIO(bytes_io)
    with open(filename, 'wt') as text_file:
        header = json.dumps(bdiff.header, sort_keys=True, indent=4)
        text_file.write('%s\n' % header)
        text_file.write('snvs\t%d\n' % bdiff.snv_count)
        text_file.write('indels\t%d\n' % bdiff.indel_count)
        
        text_file.write('last_index\t%d\n' % (bdiff.last_index + 1))
        text_file.write('index_res\t%d\n' % bdiff.index_resolution)
        
        text_file.write('index_size\t%d\n' % len(bdiff.file_index))
        for index, pos in bdiff.file_index:
            text_file.write('%d\t%d' % (index + 1, pos))
        
        for index, alts in bdiff:
            if isinstance(alts, tuple):
                text_file.write('%d\t%s\n' % (index + 1, ','.join(alts)))
            else:  # is indel
                str_alts = [':'.join(alt) for alt in alts]
                text_file.write('%d\t%s\n' % (index + 1, ','.join(str_alts)))

# def to_string_io(bytes_io: io.BytesIO()):
#     """
#     Convert Bdiff file from BytesIO to StringIO.
#     Indices are converted from 0-based to 1-based.
#     :param bytes_io:
#     :return:
#     """
#     bdiff = BdiffIO(bytes_io)
#     text_file = io.StringIO()
#
#     header = json.dumps(bdiff.header)
#
#     text_file.write('%s\n' % header)
#     text_file.write('%d\n' % bdiff.snv_count)
#     text_file.write('%d\n' % bdiff.indel_count)
#
#     for index, alts in bdiff:
#         if isinstance(alts, tuple):
#             text_file.write('%d\t%s\n' % (index + 1, ','.join(alts)))
#         else:  # is indel
#             str_alts = [':'.join(alt) for alt in alts]
#             text_file.write('%d\t%s\n' % (index + 1, ','.join(str_alts)))
#
#     text_file.seek(0)
#     return text_file
