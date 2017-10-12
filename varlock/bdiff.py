import json
import struct
import io
import os

import math
import shutil

import varlock.common as cmn


class BdiffIO:
    """
    Bdiff is binary file containing SNV and INDEL difference records.

    Bdiff meta header:
    header_size, 4B
    header, <header_size>B
    
    Bdiff counts:
    number of SNV records, 4B
    number of INDEL records, 4B
    
    Bdiff index:
    last_index 4B, index of last record
    index length 4B, number of index records
    index resolution 4B, distance between indexed values as number of records
    list: of index length
        index of a thousandth record, 4B
        file content position of N-th record, 4B

    Diff SNV record: permutation of DNA bases (A,T,G,C)
    index 4B, 0-based absolute genomic position
    length 1B, always zero to indicate SNV record (SNV has always 4 alternatives)
    reference index 1B, index of reference allele in the A,T,G,C list
    mapping, 1B, index of base permutation

    DIFF INDEL record: permutation of sequences
    index 4B, 0-based absolute genomic position
    length 1B, number of alternatives
    reference index 1B, index of the reference allele in the sorted list
    list:, of length, permutation of sorted (ascending) alternative sequences
        INDEL base length 2B, length of INDEL sequence
        list:
            INDEL 4 base sequence, 1B
    """
    # reserved header field, start of BAM's effective range, inclusive
    FROM_INDEX = '_from_index'
    
    # reserved header field, end of BAM's effective range, inclusive
    TO_INDEX = '_to_index'
    
    # # reserver header field
    # EFFECTIVE_RANGE = '_effective_range'
    
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
        CLass for handling BDIFF file. Instance is either in write or read mode.
        :param diff_file:  write_mode if None, read_mode otherwise
        :param index_resolution: relevant only for write_mode, read_mode uses value stored in the file
        """
        self._file_index = []
        self._header = {}
        self._snv_count = self._indel_count = 0
        self._last_index = 0
        self._data_offset = 0
        self._index_resolution = index_resolution
        self._bdiff_file = None
        self._record_file = None
        
        if diff_file is not None:
            # read mode
            diff_file.seek(0)
            header_size = struct.unpack('<I', diff_file.read(self.INT_SIZE))[0]
            self._header = cmn.bytes2dict(diff_file.read(header_size))
            self._snv_count, self._indel_count = struct.unpack('<II', diff_file.read(self.INT_SIZE * 2))
            
            # read_record index
            self._last_index, \
            index_length, \
            self._index_resolution = struct.unpack('<III', diff_file.read(self.INT_SIZE * 3))
            
            for i in range(index_length):
                index, pos = struct.unpack('<II', diff_file.read(self.INT_SIZE * 2))
                self._file_index.append((index, pos))
            
            self._data_offset = diff_file.tell()
            # file pointer is at start of the content
            self._bdiff_file = diff_file
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
        """
        :return: list of tuples (genomic_index, file_pos)
            file_pos is position relative to data_offset.
            Add data_offset to obtain absolute file position.
        """
        return self._file_index.copy()
    
    @property
    def index_resolution(self):
        return self._index_resolution
    
    @property
    def is_read_mode(self):
        return self._bdiff_file is not None
    
    @property
    def data_offset(self):
        """
        :return: header size in bytes
        """
        return self._data_offset
    
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
        :return: genomic index of the first record
        """
        if self.is_empty():
            return None
        
        if self.is_read_mode:
            curr_pos = self._bdiff_file.tell()
            self._bdiff_file.seek(self._data_offset)
            first_index = struct.unpack('<I', self._bdiff_file.read(self.INT_SIZE))[0]
            self._bdiff_file.seek(curr_pos)
        else:
            curr_pos = self._record_file.tell()
            self._record_file.seek(0)
            first_index = struct.unpack('<I', self._record_file.read(self.INT_SIZE))[0]
            self._record_file.seek(curr_pos)
        
        return first_index
    
    @property
    def last_index(self):
        """
        :return: genomic index of the last record
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
        byte_str = self._bdiff_file.read(self.INT_SIZE + self.BYTE_SIZE * 2)
        if byte_str == b'':
            raise EOFError
        
        index, length, ref_id = struct.unpack('<IBB', byte_str)
        if length == 0:
            # SNV record
            return index, ref_id, self._read_snv_alts()
        else:
            # INDEL record
            return index, ref_id, self._read_indel_alts(length)
    
    def read_record(self):
        """
        :return: genomic_index, is_indel, reference_sequence, mutation_map
        """
        index, ref_id, alts = self._read_record()
        if isinstance(alts, tuple):
            # is SNV
            return index, False, cmn.BASES[ref_id], dict(zip(alts, cmn.BASES))
        else:
            # is INDEL
            return index, True, alts[ref_id], dict(zip(alts, sorted(alts)))
    
    def _read_snv_alts(self):
        """
        Read SNV alternatives from rest of the record
        :return: (base_1, base_2, base_3, base_4)
        """
        mut_index = struct.unpack('<B', self._bdiff_file.read(self.BYTE_SIZE))[0]
        return self.INDEX_2_PERM[mut_index]
    
    def _read_indel_alts(self, length):
        """
        Read INDEL alternatives from rest of the record
        :return: [SEQ_1, SEQ_2, ...]
        """
        alt_list = [None] * length
        for i in range(length):
            base_length = struct.unpack('<H', self._bdiff_file.read(self.SHORT_SIZE))[0]
            seq_byte_size = math.ceil(base_length / 4)
            seq = cmn.bytes2seq(self._bdiff_file.read(seq_byte_size), base_length)
            alt_list[i] = seq
        
        return alt_list
    
    def _write_snv(self, index: int, ref_id: int, base_perm: tuple):
        """
        :param index: genomic index
        :param base_perm: permutation of bases
        permutation represent values of map with keys A,T,G,C in this order
        :return:
        """
        assert not self.is_read_mode
        # assert ref_id >= 0 < len(cmn.BASES)
        assert index > self._last_index
        
        record = struct.pack('<IBBB', index, 0, ref_id, self.PERM_2_INDEX[base_perm])
        self._record_file.write(record)
        
        self._last_index = index
        self._snv_count += 1
        if (self._snv_count + self._indel_count) % self._index_resolution == 0:
            self._file_index.append((index, self._record_file.tell() - len(record)))
    
    def write_snv(self, index: int, ref_id: int, mut_map: dict):
        """
        :param index:
        :param ref_id: index of reference allele in the A,T,G,C list
        :param mut_map:
        :return:
        """
        self._write_snv(index, ref_id, (mut_map['A'], mut_map['T'], mut_map['G'], mut_map['C']))
    
    def _write_indel(self, index: int, ref_id: int, seq_perm: list):
        """
        :param index:
        :param seq_perm:
        :return:
        """
        assert not self.is_read_mode
        assert index > self._last_index
        assert 1 < len(seq_perm) < 256
        
        record = struct.pack('<IBB', index, len(seq_perm), ref_id)
        for alt in seq_perm:
            record += struct.pack('<H', len(alt)) + cmn.seq2bytes(alt)
        
        self._record_file.write(record)
        
        self._last_index = index
        self._indel_count += 1
        if (self._snv_count + self._indel_count) % self._index_resolution == 0:
            self._file_index.append((index, self._record_file.tell() - len(record)))
    
    def write_indel(self, index: int, ref_seq: str, mut_map: dict):
        """
        :param index: genomic index
        :param ref_seq:
        :param mut_map: original sequences -> permuted sequences
        """
        seq_perm = self.seq_perm(mut_map)
        # noinspection PyTypeChecker
        self._write_indel(index, seq_perm.index(ref_seq), seq_perm)
    
    @staticmethod
    def seq_perm(mut_map: dict):
        """
        Sorts map by its keys and returns list of its values in respective order.
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
    
    @classmethod
    def _write_header(cls, diff_file: io.BytesIO, header: dict):
        meta_bytes = cmn.dict2bytes(header)
        diff_file.write(struct.pack('<I', len(meta_bytes)))
        diff_file.write(meta_bytes)
    
    @classmethod
    def _write_counts(cls, diff_file: io.BytesIO, snv_count: int, indel_count: int):
        diff_file.write(struct.pack('<II', snv_count, indel_count))
    
    @classmethod
    def _write_index(cls, diff_file: io.BytesIO(), file_index: list, last_index: int, index_resolution: int):
        # file index header
        diff_file.write(struct.pack(
            '<III',
            last_index,
            len(file_index),
            index_resolution
        ))
        # file index content
        for index, pos in file_index:
            diff_file.write(struct.pack('<II', index, pos))
    
    def _file_from_write(self, header: dict):
        assert not self.is_read_mode
        
        bdiff_file = io.BytesIO()
        self._write_header(bdiff_file, header)
        self._write_counts(bdiff_file, self._snv_count, self._indel_count)
        self._write_index(bdiff_file, self._file_index, self._last_index, self._index_resolution)
        
        self._record_file.seek(0)
        shutil.copyfileobj(self._record_file, bdiff_file)
        
        return bdiff_file
    
    def _file_from_copy(self):
        assert self.is_read_mode
        
        bdiff_file = io.BytesIO()
        self._bdiff_file.seek(0)
        shutil.copyfileobj(self._bdiff_file, bdiff_file)
        
        return bdiff_file
    
    def _build_index(self, from_index_pos: int, to_index: int):
        assert self.is_read_mode
        
        self._bdiff_file.seek(from_index_pos)
        
        file_index = []
        snv_count = 0
        indel_count = 0
        while True:
            try:
                index, ref_id, alts = self._read_record()
            except EOFError:
                break
            
            if index > to_index:
                break
            
            if isinstance(alts, tuple):
                # SNV
                snv_count += 1
            else:
                # INDEL
                indel_count += 1
            
            if (snv_count + indel_count) % self._index_resolution == 0:
                file_index.append((index, self._bdiff_file.tell()))
        
        return file_index, snv_count, indel_count
    
    def _file_from_slice(self, header: dict):
        assert self.is_read_mode
        
        bdiff_file = io.BytesIO()
        self._write_header(bdiff_file, header)
        
        if self.FROM_INDEX in header and self.TO_INDEX in header:
            from_index_pos, to_index_pos = self.tell_range(header[self.FROM_INDEX], header[self.TO_INDEX])
            file_index, snv_count, indel_count = self._build_index(from_index_pos, header[self.TO_INDEX])
            
            self._bdiff_file.seek(to_index_pos)
            last_index, ref_id, alts = self._read_record()
            end_pos = self._bdiff_file.tell()
            
            self._write_counts(bdiff_file, snv_count, indel_count)
            self._write_index(bdiff_file, file_index, last_index, self._index_resolution)
            
            self._bdiff_file.seek(from_index_pos)
            # TODO optimalize writing (use buffer)
            bdiff_file.write(self._bdiff_file.read(end_pos - from_index_pos))
            
            shutil.copyfileobj(self._bdiff_file, bdiff_file)
        
        else:
            # unmapped only
            self._write_counts(bdiff_file, 0, 0)
            self._write_index(bdiff_file, [], 0, self._index_resolution)
        
        return bdiff_file
    
    def file(self, header: dict = None, close=True):
        """
        Creates Bdiff file using internal data.
        :param header: meta header dict
        :param close: closes internal stream in write mode or input stream in read mode
        :return: Bdiff file as BytesIO
        """
        if self.is_read_mode:
            prev_pos = self._bdiff_file.tell()
            if header is None:
                # exact copy
                diff_file = self._file_from_copy()
            else:
                # header does affect the final content
                diff_file = self._file_from_slice(header)
            
            if close:
                self._bdiff_file.close()
            else:
                self._bdiff_file.seek(prev_pos, os.SEEK_SET)
        
        else:  # write mode
            prev_pos = self._record_file.tell()
            if header is None:
                diff_file = self._file_from_write({})
            else:
                # header does not affect the final content
                diff_file = self._file_from_write(header)
            
            if close:
                self._record_file.close()
            else:
                self._record_file.seek(prev_pos, os.SEEK_SET)
        
        diff_file.seek(0)
        return diff_file
    
    def tell(self):
        if self.is_read_mode:
            return self._bdiff_file.tell()
        else:
            return self._record_file.tell()
    
    def seek(self, pos: int, mode: int = os.SEEK_SET):
        assert self.is_read_mode
        self._bdiff_file.seek(pos, mode)
    
    def tell_range(self, from_index, to_index):
        """
        Genomic index range of file records between lower_index and upper_index.
        :param from_index: from index inclusive
        :param to_index: to index inclusive
        :return: (from_index_pos, to_index_pos)
            from_index_pos is file position of first record with index greater than or equal to from_index
            to_index_pos is position of last record with index lower than or equal to to_index
        :raises: IndexError
        """
        assert self.is_read_mode
        assert from_index <= to_index
        
        from_index_pos = self.tell_index_gte(from_index)
        to_index_pos = self.tell_index_lte(to_index)
        
        if from_index_pos is None:
            raise IndexError("Supplied range is empty")
        
        if to_index_pos is None:
            raise IndexError("Supplied range is empty")
        
        return from_index_pos, to_index_pos
    
    def tell_index_gte(self, pivot_index: int):
        """
        :param pivot_index:
        :return: the first file position of index greater than or equal to pivot_index
        None if such index does not exists.
        """
        assert self.is_read_mode
        
        if self.last_index is None or pivot_index > self.last_index:
            # file is empty or desired index does not exists
            return None
        else:
            saved_pos = self._bdiff_file.tell()
            curr_pos = self._indexed_pos(pivot_index)
            self._bdiff_file.seek(curr_pos)
            
            curr_index = self._read_record()[0]  # type: int
            while curr_index < pivot_index:
                curr_pos = self._bdiff_file.tell()
                curr_index = self._read_record()[0]  # type: int
            
            self._bdiff_file.seek(saved_pos)
            return curr_pos
    
    def tell_index_lte(self, pivot_index: int):
        """
        :param pivot_index:
        :return: the last file position of index lower than or equal to pivot_index
        None if such index does not exists.
        """
        assert self.is_read_mode
        
        if self.first_index is None or pivot_index < self.first_index:
            # desired index does not exists
            return None
        else:
            saved_pos = self._bdiff_file.tell()
            curr_pos = self._indexed_pos(pivot_index)
            self._bdiff_file.seek(curr_pos)
            
            while True:
                next_pos = self._bdiff_file.tell()
                try:
                    next_index = self._read_record()[0]  # type: int
                except EOFError:
                    break
                if next_index <= pivot_index:
                    curr_pos = next_pos
                else:
                    break
            
            self._bdiff_file.seek(saved_pos)
            return curr_pos
    
    def _indexed_pos(self, pivot_index: int):
        """
        :param pivot_index:
        :return: index position of record with index lower than or equal to pivot_index.
            Zero position is implicitly indexed and returned if indexed position is not found.
            Returned position is absolute within the file.
        """
        indexed_pos = 0
        for index, pos in self._file_index:
            if index <= pivot_index:
                indexed_pos = pos
            else:
                break
        
        return indexed_pos + self._data_offset
    
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
                text_file.write('%d\t%d\n' % (index + 1, pos))
            
            for index, ref_id, alts in bdiff:
                is_indel = isinstance(alts, list)
                text_file.write('%d\t%d\t%d\t%s\n' % (index + 1, is_indel, ref_id, ','.join(alts)))
    
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
            
            while True:
                line = text_file.readline().rstrip()
                if not line:
                    break
                raw_index, raw_is_indel, raw_ref_id, raw_alts = line.split('\t')
                index = int(raw_index) - 1
                is_snv = raw_is_indel == '0'
                ref_id = int(raw_ref_id)
                alts = raw_alts.split(',')
                if is_snv:
                    # snvs
                    bdiff._write_snv(index, ref_id, tuple(alts))
                else:
                    # indel
                    bdiff._write_indel(index, ref_id, alts)
        
        return bdiff.file(header)
