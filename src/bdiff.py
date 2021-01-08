import io
import json
import math
import os
import shutil
import struct

import src.common as cmn
from src.po import ZygosityChange


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

    Diff SNV record: single or two permutations of DNA bases (A,T,G,C)
    index, 4B, 0-based absolute genomic position
    length, 1B, always zero to indicate SNV record (always 4 possible alleles)
    reference index, 1B, index of reference allele within the A,T,G,C list
    zygosity, 1B, index from (homo2homo, homo2hetero, hetero2hetero, hetero2homo)
    mapping_a, 1B, permutation of alleles as indices of the A,T,G,C list
        data below is is only present when zygosity is homo2hetero or hetero2homo
    mapping_b, 1B, permutation of alleles as indiced of the A,T,G,C list
    seed, 4B, seed for RNG used to determine alleles with mapping_b

    DIFF INDEL record: permutation of sequences
    index, 4B, 0-based absolute genomic position
    length, 1B, number of alleles
    reference index, 1B, index of the reference allele within the list of alleles
    list (of length): alleles
        index, 1B, index of target allele
        INDEL sequence length, 2B
        # list: INDEL sequence
        #     INDEL base, 1B
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

    def __init__(self, diff_file: io.BytesIO = None, index_resolution=1000):
        """
        Class for handling BDIFF file. Instance is either in write or read mode.
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

    # TODO make public
    def _read_record(self) -> (int, bool, int, int, list, list, int):
        """
        Reads next SNV or INDEL record
        :return: index, is_snv, reference_id, zygosity_byte, permutation_a, permutation_b, rng_seed
        """
        assert self.is_read_mode
        byte_str = self._bdiff_file.read(self.INT_SIZE + self.BYTE_SIZE * 2)
        if byte_str == b'':
            raise EOFError

        index, length, ref_id = struct.unpack('<IBB', byte_str)
        is_snv = length == 0
        if is_snv:
            zygosity_byte, mapping_a_byte = struct.unpack('<BB', self._bdiff_file.read(self.BYTE_SIZE * 2))
            if zygosity_byte == ZygosityChange.HOMO_TO_HETERO.value \
                    or zygosity_byte == ZygosityChange.HETERO_TO_HOMO.value:
                # use double mapping
                mapping_b_byte, rng_seed = struct.unpack('<BI', self._bdiff_file.read(self.BYTE_SIZE + self.INT_SIZE))
                return index, is_snv, ref_id, zygosity_byte, \
                       cmn.byte2base_perm(mapping_a_byte), cmn.byte2base_perm(mapping_b_byte), rng_seed
            else:
                # use single mapping
                permutation = cmn.byte2base_perm(mapping_a_byte)
                return index, is_snv, ref_id, zygosity_byte, permutation, permutation, 0
        else:
            # TODO zygosity
            permutation = self._read_indel_alts(length)
            return index, is_snv, ref_id, 0, permutation, permutation, 0

    # def read_record(self) -> (int, VariantType, str, ZygosityChange, dict, dict):
    #     """
    #     :return: genomic_index, is_snv, reference_sequence, zygosity, mutation_map_a, mutation_map_b
    #     """
    #     index, is_snv, ref_id, zygosity, perm_a, perm_b, rng_seed = self._read_record()
    #     if is_snv:
    #         # SNV: stored alleles are the original ones
    #         return index, VariantType.SNV, cmn.BASES[ref_id], ZygosityChange(zygosity), \
    #                dict(zip(cmn.BASES, perm_a)), dict(zip(cmn.BASES, perm_b)), rng_seed
    #     else:
    #         # INDEL it is assumed that record was saved as:
    #         # sorted(alleles) -> permuted(alleles)
    #         # TODO zygosity
    #         return index, VariantType.INDEL, perm_a[ref_id], zygosity, \
    #                dict(zip(perm_a, sorted(perm_a))), dict(zip(perm_b, sorted(perm_b))), rng_seed

    def _read_indel_alts(self, length):
        """
        Read mapped alleles from rest of an INDEL record
        :return: [SEQ_1, SEQ_2, ...]
        """
        alt_list = [None] * length
        for i in range(length):
            base_length = struct.unpack('<H', self._bdiff_file.read(self.SHORT_SIZE))[0]
            seq_byte_size = math.ceil(base_length / 4)
            seq = cmn.bytes2seq(self._bdiff_file.read(seq_byte_size), base_length)
            alt_list[i] = seq

        return alt_list

    def _write_record(self, record: bytes, index: int):
        self._record_file.write(record)
        self._last_index = index
        if (self._snv_count + self._indel_count) % self._index_resolution == 0:
            self._file_index.append((index, self._record_file.tell() - len(record)))

    # TODO make public
    def _write_snv(self, index: int, ref_id: int, zygosity: ZygosityChange, perm_a: list, perm_b: list, rng_seed: int):
        """
        :param index: genomic index
        :param base_perm: permutation of bases
        permutation represent values of map with keys A,T,G,C in respective order
        :return:
        """
        assert not self.is_read_mode
        assert ref_id >= 0 < len(cmn.BASES)
        assert index > self._last_index
        record = struct.pack('<IBBBB', index, 0, ref_id, zygosity.value, cmn.base_perm2byte(perm_a))

        if zygosity.is_changed():
            # add second mapping
            record += struct.pack('<BI', cmn.base_perm2byte(perm_b), rng_seed)
        else:
            # single mapping
            assert perm_a == perm_b

        self._snv_count += 1
        self._write_record(record, index)

    # def write_snv(
    #         self,
    #         index: int,
    #         ref_id: int,
    #         zygosity: ZygosityChange,
    #         mask_map_a: dict,
    #         mask_map_b: dict,
    #         rng_seed: int
    # ):
    #     """
    #     :param rng_seed:
    #     :param mask_map_b:
    #     :param mask_map_a:
    #     :param zygosity:
    #     :param index:
    #     :param ref_id: index of reference allele in the A,T,G,C list
    #     :param mut_map:
    #     :return:
    #     """
    #     self._write_snv(
    #         index=index,
    #         ref_id=ref_id,
    #         zygosity=zygosity,
    #         perm_a=[mask_map_a['A'], mask_map_a['T'], mask_map_a['G'], mask_map_a['C']],
    #         perm_b=[mask_map_b['A'], mask_map_b['T'], mask_map_b['G'], mask_map_b['C']],
    #         rng_seed=rng_seed
    #     )

    def _write_indel(self, index: int, ref_id: int, seq_perm: list):
        """
        :param index: genomic position
        :param ref_id: index of reference allele within the alleles
        # :param alleles: list of allele sequences
        # :param perm_indices: permutation indices
        :return:
        """
        assert not self.is_read_mode
        assert index > self._last_index
        # assert 1 < len(alleles) < 256
        # assert len(alleles) == len(perm_indices)
        # assert 0 <= ref_id < len(alleles)

        # record = struct.pack('<IBB', index, len(alleles), ref_id)
        # for i in range(len(alleles)):
        #     assert 1 <= len(alleles[i]) < 256
        #     assert 0 <= perm_indices[i] < len(alleles)
        #     record += struct.pack('<BB', perm_indices[i], len(alleles[i])) + cmn.seq2bytes(alleles[i])

        record = struct.pack('<IBB', index, len(seq_perm), ref_id)
        for alt in seq_perm:
            record += struct.pack('<H', len(alt)) + cmn.seq2bytes(alt)

        self._indel_count += 1
        self._write_record(record, index)

    # TODO replace ???
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
        # same items of dict can be written in different order between calls
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
                record = self._read_record()
                index = record[0]
                is_snv = record[1]
            except EOFError:
                break

            if index > to_index:
                break

            if is_snv:
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
            last_index = self._read_record()[0]
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
        Indices are always 0-based.
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
            text_file.write('last_index\t%d\n' % bdiff.last_index)

            text_file.write('index_res\t%d\n' % bdiff.index_resolution)
            text_file.write('index_size\t%d\n' % len(bdiff.file_index))
            for index, pos in bdiff.file_index:
                text_file.write('%d\t%d\n' % (index, pos))

            for record in bdiff:
                if ZygosityChange(int(record[3])).is_changed():
                    record_list = list(record)
                    record_list[1] = int(record[1])
                    record_list[4] = ','.join(record[4])
                    record_list[5] = ','.join(record[5])
                else:
                    record_list = list(record[:5])
                    record_list[1] = int(record[1])
                    record_list[4] = ','.join(record[4])

                text_file.write('\t'.join([str(value) for value in record_list]) + '\n')

    @staticmethod
    def from_text_file(filename: str):
        """
        Convert Bdiff text file to BytesIO.
        Indices are always 0-based.
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
                segments = line.split('\t')
                if len(segments) == 7:
                    raw_index, raw_is_snv, raw_ref_id, raw_zygosity, raw_perm_a, raw_perm_b, raw_rng_seed = segments
                elif len(segments) == 5:
                    raw_index, raw_is_snv, raw_ref_id, raw_zygosity, raw_perm_a = segments
                    raw_perm_b = raw_perm_a
                    raw_rng_seed = 0
                else:
                    raise ValueError

                index = int(raw_index)
                is_snv = raw_is_snv == '1'
                ref_id = int(raw_ref_id)
                zygosity = ZygosityChange(int(raw_zygosity))
                perm_a = raw_perm_a.split(',')
                perm_b = raw_perm_b.split(',')
                rng_seed = int(raw_rng_seed)
                if is_snv:
                    # snvs
                    bdiff._write_snv(index, ref_id, zygosity, perm_a, perm_b, rng_seed)
                else:
                    # indel
                    bdiff._write_indel(index, ref_id, perm_a)

        return bdiff.file(header)
