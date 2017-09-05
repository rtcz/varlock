import json
import struct
import io
import os

import math

import varlock.common as cmn


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
    list:, of length, permutation of sorted (ascending) alternative sequences
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
            self._header = cmn.bytes2dict(diff_file.read(meta_size))
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
            return self._read_record()
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
        return self._header
    
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
    
    def _read_record(self):
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
    
    def read_record(self):
        """
        :return: genomic_index, is_indel, mutation_map
        """
        index, alts = self._read_record()
        if isinstance(alts, tuple):
            # is SNV
            return index, False, dict(zip(alts, cmn.BASES))
        else:
            # is INDEL
            return index, True, dict(zip(alts, sorted(alts)))
    
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
            base_length = struct.unpack('<H', self._diff_file.read(self.SHORT_SIZE))[0]
            seq_byte_size = math.ceil(base_length / 4)
            seq = cmn.bytes2seq(self._diff_file.read(seq_byte_size), base_length)
            alt_list[i] = seq
        
        return alt_list
    
    def _write_snv(self, index: int, base_perm: tuple):
        """
        :param index: genomic index
        :param base_perm: permutation of bases
        permutation represent values of map with keys A,T,G,C in this order
        :return:
        """
        assert not self.is_read_mode
        assert index > self._last_index
        
        record = struct.pack('<IBB', index, 0, self.PERM_2_INDEX[base_perm])
        self._record_file.write(record)
        
        self._last_index = index
        self._snv_count += 1
        if self._snv_count % self._index_resolution == 0:
            self._file_index.append((index, self._record_file.tell() - len(record)))
    
    def write_snv(self, index: int, mut_map: dict):
        self._write_snv(index, (mut_map['A'], mut_map['T'], mut_map['G'], mut_map['C']))
    
    def _write_indel(self, index: int, seq_perm: list):
        """
        :param index:
        :param seq_perm: permutation of sequences
        :return:
        """
        assert not self.is_read_mode
        assert index > self._last_index
        assert 1 < len(seq_perm) < 256
        
        record = struct.pack('<IB', index, len(seq_perm))
        for alt in seq_perm:
            record += struct.pack('<H', len(alt)) + cmn.seq2bytes(alt)
        
        self._record_file.write(record)
        
        self._last_index = index
        self._indel_count += 1
        if self._indel_count % self._index_resolution == 0:
            self._file_index.append((index, self._record_file.tell() - len(record)))
    
    def write_indel(self, index: int, mut_map: dict):
        """
        :param index: genomic index
        :param mut_map: original sequences -> permuted sequences
        """
        self._write_indel(index, self.seq_perm(mut_map))
    
    @staticmethod
    def seq_perm(mut_map: dict):
        """
        :param mut_map:
        :return: sequence permutation list
        """
        alts = [None] * len(mut_map)
        i = 0
        for key, value in sorted(mut_map.items(), key=lambda x: x[0]):
            alts[i] = value
            i += 1
        return alts
    
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
            meta_bytes = cmn.dict2bytes({} if header is None else header)
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
    
    def tell(self):
        if self.is_read_mode:
            return self._diff_file.tell()
        else:
            return self._record_file.tell()
    
    def seek(self, pos: int, mode: int = os.SEEK_SET):
        assert self.is_read_mode
        self._diff_file.seek(pos, mode)
    
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
            
            curr_index = self._read_record()[0]  # type: int
            while curr_index < pivot_index:
                curr_pos = self._diff_file.tell()
                curr_index = self._read_record()[0]  # type: int
            
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
                    next_index = self._read_record()[0]  # type: int
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
    
    @staticmethod
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
                is_indel = isinstance(alts, list)
                text_file.write('%d\t%d\t%s\n' % (index + 1, is_indel, ','.join(alts)))
    
    @staticmethod
    def from_text_file(filename: str):
        """
        Convert Bdiff text file to BytesIO.
        Indices are converted from 1-based to 0-based.
        :param filename:
        :return:
        """
        with open(filename, 'rt') as text_file:
            # header
            header_str = ''
            line = text_file.readline()
            header_str += line
            while line != '}\n':
                line = text_file.readline()
                header_str += line
            header = json.loads(header_str)
            # snvs
            text_file.readline()
            # indels
            text_file.readline()
            # last index
            text_file.readline()
            # index resolution
            index_res = int(text_file.readline().rstrip().split('\t')[1])
            # index size
            index_size = int(text_file.readline().rstrip().split('\t')[1])
            for i in range(index_size):
                text_file.readline()
            
            bdiff = BdiffIO(index_resolution=index_res)
            line = text_file.readline().rstrip()
            while line != '':
                raw_index, raw_is_indel, raw_alts = line.split('\t')
                index = int(raw_index) - 1
                is_snv = raw_is_indel == '0'
                alts = raw_alts.split(',')
                if is_snv:
                    # snvs
                    bdiff._write_snv(index, tuple(alts))
                else:
                    # indel
                    bdiff._write_indel(index, alts)
                
                line = text_file.readline().rstrip()
        
        return bdiff.file(header)
