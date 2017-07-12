import io
import struct

from .common import bin2hex, hex2bin


class Diff:
    """
    Diff is binary file, where each record represent one SNV mapping from original BAM to mutated BAM.
    
    Diff record:
    index, 4B, absolute position of SNV in genome
    mapping, 1B, id of bases permutation

    Diff header:
    md5 checksum of mutated BAM file, 16B
    start, first genomic index of DIFF range, 4B
    end, last genomic index of DIFF range, 4B
    unmapped secret, 64B
    """
    RECORD_FORMAT = "<IB"  # int, char
    RECORD_SIZE = 5
    
    CHECKSUM_SIZE = 16
    INT_SIZE = 4
    SECRET_SIZE = 16
    SECRET_PLACEHOLDER = bytes(SECRET_SIZE)
    
    HEADER_SIZE = CHECKSUM_SIZE + INT_SIZE + INT_SIZE + SECRET_SIZE
    
    INDEX_2_MUT = [
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
    
    MUT_2_INDEX = {
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
    
    @classmethod
    def read_header(cls, diff_file):
        """
        :param diff_file:
        :return:
        """
        diff_file.seek(0)
        checksum = diff_file.read(cls.CHECKSUM_SIZE)
        start_index, end_index = struct.unpack('<II', diff_file.read(cls.INT_SIZE * 2))
        secret = diff_file.read(cls.SECRET_SIZE)
        return checksum, start_index, end_index, secret
    
    @classmethod
    def validate_header_range(cls, diff_file):
        checksum, header_start_index, header_end_index, secret = cls.read_header(diff_file)
        content_start_index = cls.__read_index(diff_file)
        diff_file.seek(-cls.RECORD_SIZE, 2)
        content_end_index = cls.__read_index(diff_file)
        
        if header_start_index > content_start_index:
            raise ValueError("Invalid header start index")
        
        if header_end_index < content_end_index:
            raise ValueError("Invalid header end index")
    
    @classmethod
    def write_header(
            cls,
            diff_file,
            start_index: int,
            end_index: int,
            bam_checksum: bytes = bytes(CHECKSUM_SIZE),
            secret: bytes = bytes(SECRET_SIZE)
    ):
        assert start_index <= end_index
        cls.write_checksum(diff_file, bam_checksum)
        diff_file.write(struct.pack('<II', start_index, end_index))
        cls.write_secret(diff_file, secret)
    
    @classmethod
    def read_record(cls, diff_file):
        byte_string = diff_file.read(cls.RECORD_SIZE)
        if len(byte_string) == 0:
            raise EOFError()
        index, mut_index = struct.unpack(cls.RECORD_FORMAT, byte_string)
        return index, cls.INDEX_2_MUT[mut_index]
    
    @classmethod
    def write_record(cls, diff_file, index: int, mut_tuple: tuple):
        byte_string = struct.pack(cls.RECORD_FORMAT, index, cls.MUT_2_INDEX[mut_tuple])
        diff_file.write(byte_string)
    
    @classmethod
    def write_checksum(cls, diff_file, bam_checksum: bytes):
        if len(bam_checksum) != cls.CHECKSUM_SIZE:
            raise ValueError('Invalid checksum size')
        
        diff_file.write(bam_checksum)
    
    @classmethod
    def write_secret(cls, diff_file, secret: bytes):
        if len(secret) != cls.SECRET_SIZE:
            raise ValueError('Invalid secret size')
        
        diff_file.write(secret)
    
    @classmethod
    def get_start_pos(cls, diff_file):
        diff_file.seek(cls.HEADER_SIZE, 0)
        return diff_file.tell()
    
    @classmethod
    def get_end_pos(cls, diff_file):
        diff_file.seek(0, 2)
        return diff_file.tell()
    
    @classmethod
    def seek_subrange(cls, diff_file, start_index: int, end_index: int):
        assert start_index <= end_index
        
        start_offset = cls.seek_closest_index(diff_file, start_index, True)
        end_offset = cls.seek_closest_index(diff_file, end_index, False)
        
        return start_offset, end_offset + cls.RECORD_SIZE
    
    @classmethod
    def __read_index(cls, diff_file):
        index = struct.unpack('<I', diff_file.read(cls.INT_SIZE))[0]
        # diff_file.seek(-cls.INT_LENGTH, 1)
        return index
    
    @classmethod
    def body_size(cls, diff_file):
        diff_file.seek(0, 2)
        return diff_file.tell() - cls.HEADER_SIZE
    
    @classmethod
    def validate(cls, diff_file):
        body_size = cls.body_size(diff_file)
        if body_size < 0:
            raise IOError("Diff too short")
        
        if body_size == 0:
            # file is empty - return position after header
            raise EOFError("Diff is empty")
        
        if body_size % cls.RECORD_SIZE != 0:
            raise IOError("Diff file has invalid number of bytes")
        
        cls.validate_header_range(diff_file)
    
    @classmethod
    def seek_closest_index(cls, diff_file, index: int, up=True):
        """
        :param diff_file: file to seek
        :param index: initial index
        :param up: direction of seeking
        :return: file position of closest index
        """
        cls.validate(diff_file)
        
        diff_file.seek(cls.HEADER_SIZE, 0)
        start_index = cls.__read_index(diff_file)
        diff_file.seek(-Diff.RECORD_SIZE, 2)
        end_index = cls.__read_index(diff_file)
        
        if index > end_index:
            # index is after DIFF content
            if up:
                raise IndexError("Upper index not found")
            else:
                diff_file.seek(-cls.RECORD_SIZE, 2)
                return diff_file.tell()
        
        if index < start_index:
            # index is before DIFF content
            if up:
                diff_file.seek(cls.HEADER_SIZE, 0)
                return diff_file.tell()
            else:
                raise IndexError("Lower index not found")
        
        first = 0
        last = cls.body_size(diff_file) / cls.RECORD_SIZE - 1
        curr_index = index
        while first <= last:
            mid = int((first + last) / 2)
            diff_file.seek(cls.HEADER_SIZE + (mid * cls.RECORD_SIZE))
            curr_index = cls.__read_index(diff_file)
            if index == curr_index:
                break
            else:
                if index < curr_index:
                    last = mid - 1
                else:
                    first = mid + 1
        
        diff_file.seek(-cls.INT_SIZE, 1)
        offset = diff_file.tell()
        
        # final move
        if index < curr_index:
            diff_file.seek(-cls.RECORD_SIZE, 1)
            prev_offset = diff_file.tell()
            prev_index = cls.__read_index(diff_file)
            if not up or index == prev_index:
                offset = prev_offset
        
        elif index > curr_index:
            diff_file.seek(cls.RECORD_SIZE, 1)
            next_offset = diff_file.tell()
            next_index = cls.__read_index(diff_file)
            if up or index == next_index:
                offset = next_offset
        
        diff_file.seek(offset)
        return offset
    
    @classmethod
    def truncate(cls, diff_file):
        """
        :param diff_file: diff to slice from
        :return: sliced diff
        """
        diff_file.seek(0)
        sliced_diff = io.BytesIO()
        sliced_diff.write(diff_file.read(cls.HEADER_SIZE))
        sliced_diff.seek(0)
        return sliced_diff
    
    @classmethod
    def slice(cls, diff_file, start_index: int, end_index: int, keep_secret=True):
        """
        :param diff_file: diff to slice from
        :param start_index: start genome pos
        :param end_index: end genome pos
        :param keep_secret: True to copy secret from source diff.
        :return: sliced diff
        """
        assert start_index <= end_index
        
        start_offset, end_offset = cls.seek_subrange(diff_file, start_index, end_index)
        sliced_diff = io.BytesIO()
        # header
        checksum, _, _, secret = cls.read_header(diff_file)
        if keep_secret:
            cls.write_header(sliced_diff, start_index, end_index, checksum, secret)
        else:
            cls.write_header(sliced_diff, start_index, end_index, checksum)
        # body
        diff_file.seek(start_offset)
        sliced_diff.write(diff_file.read(end_offset - start_offset))
        sliced_diff.seek(0)
        return sliced_diff
    
    @classmethod
    def diff2text(cls, diff_filepath, text_filepath):
        with open(diff_filepath, 'rb') as diff_file, \
                open(text_filepath, 'wt') as text_file:
            checksum, start_index, end_index, unmapped_secret = cls.read_header(diff_file)
            text_file.write('%s\n%d-%d\n%s\n' % (bin2hex(checksum), start_index, end_index, bin2hex(unmapped_secret)))
            while True:
                try:
                    index, mut_tuple = cls.read_record(diff_file)
                    text_file.write('%d\t%s\n' % (index, ','.join(mut_tuple)))
                except EOFError:
                    break
    
    @classmethod
    def text2diff(cls, text_filepath, diff_filepath):
        with open(diff_filepath, 'wb') as diff_file, \
                open(text_filepath, 'rt') as text_file:
            checksum_hex = text_file.readline().rstrip()
            start_index, end_index = text_file.readline().rstrip().split('-')
            secret = text_file.readline().rstrip()
            cls.write_header(diff_file, int(start_index), int(end_index), hex2bin(checksum_hex),
                             hex2bin(secret))
            for line in text_file:
                index, mut = line.rstrip().split('\t')
                Diff.write_record(diff_file, int(index), tuple(mut.split(',')))
