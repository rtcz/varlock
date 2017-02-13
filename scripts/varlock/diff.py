import io
import struct


class Diff:
    """
    Diff is binary file, where each record represent one SNV mapping from original BAM to mutated BAM.
    
    Diff record:
    index, 4B, absolute position of SNV in genome
    mapping, 1B, id of bases permutation

    Diff header:
    md5 checksum of mutated BAM file, 16B
    start, first genome index of DIFF range, 4B
    end, last genome index of DIFF range, 4B
    """
    RECORD_FORMAT = "<IB"  # int, char
    RECORD_LENGTH = 5  # bytes
    INT_LENGTH = 4  # bytes
    
    HEADER_LENGTH = 24
    MD5_LENGTH = 16
    
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
        diff_file.seek(0)
        checksum = diff_file.read(cls.MD5_LENGTH)
        start_index, end_index = struct.unpack('<II', diff_file.read(cls.INT_LENGTH * 2))
        return checksum, start_index, end_index
    
    @classmethod
    def validate_header_range(cls, diff_file):
        checksum, header_start_index, header_end_index = cls.read_header(diff_file)
        content_start_index = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
        diff_file.seek(-cls.RECORD_LENGTH, 2)
        content_end_index = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
        
        if header_start_index > content_start_index:
            raise ValueError("Invalid header start index")
        
        if header_end_index < content_end_index:
            raise ValueError("Invalid header end index")
    
    @classmethod
    def write_header(cls, diff_file, bam_checksum, start_index, end_index):
        assert start_index <= end_index
        
        if len(bam_checksum) != cls.MD5_LENGTH:
            raise ValueError('Invalid checksum length')
        
        diff_file.write(bam_checksum)
        diff_file.write(struct.pack('<II', start_index, end_index))
    
    @classmethod
    def read_record(cls, diff_file):
        byte_string = diff_file.read(cls.RECORD_LENGTH)
        if len(byte_string) == 0:
            raise EOFError()
        index, mut_index = struct.unpack(cls.RECORD_FORMAT, byte_string)
        return index, cls.INDEX_2_MUT[mut_index]
    
    @classmethod
    def write_record(cls, diff_file, index, mut_tuple):
        byte_string = struct.pack(cls.RECORD_FORMAT, index, cls.MUT_2_INDEX[mut_tuple])
        diff_file.write(byte_string)
    
    @classmethod
    def get_start_pos(cls, diff_file):
        diff_file.seek(cls.HEADER_LENGTH, 0)
        return diff_file.tell()
    
    @classmethod
    def get_end_pos(cls, diff_file):
        diff_file.seek(0, 2)
        return diff_file.tell()
    
    @classmethod
    def seek_range(cls, diff_file, start_index, end_index):
        assert start_index <= end_index
        
        start_pos = cls.seek(diff_file, start_index)
        end_pos = cls.seek(diff_file, end_index)
        
        if start_pos == end_pos == cls.get_start_pos(diff_file):
            raise IndexError("Range is out of DIFF")
        
        if start_pos == end_pos == cls.get_end_pos(diff_file):
            raise IndexError("Range is out of DIFF")
        
        return start_pos, end_pos
    
    @classmethod
    def seek(cls, diff_file, index):
        """
        Seek position in DIFF file at or right after specified index.
        :param diff_file:
        :param index:
        :return: seeked position
        :raises: IndexError
        """
        # TODO move to diff validation
        diff_file.seek(0, 2)
        full_size = diff_file.tell() - cls.HEADER_LENGTH
        if full_size < 0:
            raise IOError("Diff too short")
        
        if full_size == 0:
            # file is empty - return position after header
            return cls.HEADER_LENGTH
        
        if full_size % cls.RECORD_LENGTH != 0:
            raise IOError("Diff file has invalid number of bytes")
        
        diff_file.seek(cls.HEADER_LENGTH, 0)
        start_diff = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
        diff_file.seek(-Diff.RECORD_LENGTH, 2)
        end_diff = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
        
        if index > end_diff:
            # index is after DIFF content
            # go to EOF
            diff_file.seek(0, 2)
            return diff_file.tell()
        
        if index < start_diff:
            # index is before DIFF content
            # go to SOF
            diff_file.seek(cls.HEADER_LENGTH, 0)
            return diff_file.tell()
        
        # mind diff with one record
        if full_size == cls.RECORD_LENGTH:
            diff_file.seek(cls.HEADER_LENGTH, 0)
            return diff_file.tell()
        
        # range start is inside diff content
        curr_index = index
        size = int(full_size / 2)
        while size >= cls.RECORD_LENGTH:
            diff_file.seek(cls.HEADER_LENGTH + size, 0)
            curr_index = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
            if curr_index > index:
                size -= int(size / 2)
            elif curr_index < index:
                size += int(size / 2)
            else:
                break
        
        diff_file.seek(-cls.INT_LENGTH, 1)
        
        # final move
        if curr_index > index:
            diff_file.seek(-cls.RECORD_LENGTH, 1)
        elif curr_index < index:
            diff_file.seek(+cls.RECORD_LENGTH, 1)
        
        return diff_file.tell()
    
    @classmethod
    def slice(cls, diff_file, start_index, end_index):
        """
        :param diff_file: diff to slice from
        :param start_index: start genome pos
        :param end_index: end genome pos
        :return:
        """
        assert start_index <= end_index
        
        start_pos, end_pos = cls.seek_range(diff_file, start_index, end_index)
        
        sliced_diff = io.BytesIO()
        
        # header
        checksum = cls.read_header(diff_file)[0]
        cls.write_header(diff_file, checksum, start_index, end_index)
        # content
        diff_file.seek(start_pos)
        sliced_diff.write(diff_file.read(end_pos - start_pos))
        return sliced_diff
        
        # @classmethod
        # def diff2text(cls, diff_filepath, text_filepath):
        #     with open(diff_filepath, "rb") as diff_file, \
        #             open(text_filepath, "wt") as text_file:
        #         while True:
        #             try:
        #                 index, mut_tuple = cls.read_record(diff_file)
        #                 text_file.write("%d\t%s\n" % (index, mut_tuple))
        #             except EOFError:
        #                 break
