import struct


class Diff:
    STRUCT_FORMAT = "<IB"  # int, char
    STRUCT_LENGTH = 5  # bytes
    INT_LENGTH = 4  # bytes
    
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
    def read_next(cls, diff_file):
        byte_string = diff_file.read(cls.STRUCT_LENGTH)
        if len(byte_string) == 0:
            raise EOFError()
        index, mut_index = struct.unpack(cls.STRUCT_FORMAT, byte_string)
        return index, cls.INDEX_2_MUT[mut_index]
    
    @classmethod
    def write_next(cls, diff_file, index, mut_tuple):
        byte_string = struct.pack(cls.STRUCT_FORMAT, index, cls.MUT_2_INDEX[mut_tuple])
        diff_file.write(byte_string)
    
    @classmethod
    def seek_index(cls, diff_file, start_range, end_range):
        diff_file.seek(0, 0)
        start_diff = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
        diff_file.seek(-Diff.STRUCT_LENGTH, 2)
        end_diff = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
        
        if start_range > end_diff or end_range < start_diff:
            # diff content is either before or after range
            # go to EOF
            diff_file.seek(0, 2)
            return diff_file.tell()
        
        if start_range < start_diff:
            # first index is inside diff content
            # go to SOF
            diff_file.seek(0, 0)
            return diff_file.tell()
        
        # range start is inside diff content
        
        diff_file.seek(0, 2)
        size = int(diff_file.tell() / 2)
        
        while size >= cls.STRUCT_LENGTH:
            diff_file.seek(size)
            curr_index = struct.unpack('<I', diff_file.read(cls.INT_LENGTH))[0]
            if curr_index > start_range:
                size -= int(size / 2)
            elif curr_index < start_range:
                size += int(size / 2)
            else:
                break
        
        # final move
        if curr_index > start_range:
            diff_file.seek(-cls.STRUCT_LENGTH, 1)
        elif curr_index < start_range:
            diff_file.seek(+cls.STRUCT_LENGTH, 1)
        
        diff_file.seek(-cls.INT_LENGTH, 1)
        return diff_file.tell()
    
    @classmethod
    def diff2text(cls, diff_filepath, text_filepath):
        with open(diff_filepath, "rb") as diff_file, \
                open(text_filepath, "wt") as text_file:
            while True:
                try:
                    index, mut_tuple = cls.read_next(diff_file)
                    text_file.write("%d\t%s\n" % (index, mut_tuple))
                except EOFError:
                    break
                    
                    # @classmethod
                    # def file2list(cls, diff_file):
                    #     diff_list = []
                    #     while True:
                    #         byte_string = diff_file.read(cls.STRUCT_LENGTH)
                    #         if byte_string == '':
                    #             break
                    #         record = Diff.bytes2record(byte_string)
                    #         diff_list.append(record)
                    #
                    #     return diff_list
