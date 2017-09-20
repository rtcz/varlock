import struct

from .common import *
import os


# TODO refactor like file object interface
class Vac:
    # TODO add header with VCF's checksum ?
    """
    Class for handling VAC file. VAC stands for Variant Allele Count.
    VAC is binary file, consisting of two blocks.
    The first block contains SNV records, the second one INDEL records.
    VAC uses two temporary files when writing. These files are deleted after writing is finished.
    
    VAC header:
    number of SNVs, 4B
    number of INDELs, 4B
    
    VAC SNV record:
    index 4B, 0-based absolute genomic position
    A allele count 2B
    T allele count 2B
    C allele count 2B
    G allele count 2B
    
    VAC INDEL record:
    index 4B, 0-based absolute genomic position
    length 1B, number of alternatives
    list:, sorted by allele count and sequence (descending)
        INDEL allele count 2B
        INDEL base length 2B, length of INDEL sequence
        list:
            INDEL 4 base sequence, 1B
    """
    VCF_CHROM_ID = 0
    VCF_POS_ID = 1
    VCF_REF_ID = 3
    VCF_ALT_ID = 4
    VCF_INFO_ID = 7
    
    VCF_COL_SEP = "\t"
    VCF_LIST_SEP = ","
    VCF_INFO_SEP = ";"
    
    SHORT_SIZE = 2
    INT_SIZE = 4
    # LONG_SIZE = 8
    
    HEADER_SIZE = INT_SIZE * 2
    
    MAX_ALLELE_COUNT = 2 ** 16 - 1
    """
    Max value of 2 bytes
    """
    
    SNV_FORMAT = "<IHHHH"  # int, short, short, short
    SNV_RECORD_SIZE = 12  # bytes
    
    # INDEL_FORMAT = "<IH"  # int, short
    
    SNV_TEMP_EXT = '.snv.temp'
    INDEL_TEMP_EXT = '.indel.temp'
    
    def __init__(self, fai, verbose=False):
        """
        :param fai: FastaIndex
        :param verbose:
        """
        self.fai = fai
        self.current_pos = None
        self.current_chrom = None
        self.verbose = verbose
    
    @classmethod
    def snvs_size(cls, count: int):
        return cls.SNV_RECORD_SIZE * count
    
    @staticmethod
    def is_snv(allele_list: list):
        """
        :param allele_list: variant reference and alternatives
        :return: True if allele_list represents SNP
        """
        # must contain at least one reference and one allele
        assert len(allele_list) >= 2
        
        return all(allele in BASES for allele in allele_list)
    
    @staticmethod
    def is_indel(allele_list: list):
        """
        :param allele_list: variant reference and alternatives
        :return: True if allele_list represents INDEL
        """
        # must contain at least one reference and one allele
        assert len(allele_list) >= 2
        
        seq_found = False
        for allele in allele_list:
            if len(allele) >= 2:
                # at least one allele must have multiple bases
                seq_found = True
            if any([base not in BASES for base in allele]):
                return False
        
        return seq_found
    
    @classmethod
    def parse_ac(cls, info_list: list):
        """
        Parse allele count from VCF INFO column.
        :param info_list: list of VCF INFO key=value strings
        :return: list of allele counts
        """
        for item in info_list:
            # find allele count item
            if item[:3] == 'AC=':
                # AC=value1,value2
                return list(map(int, item[3:].split(cls.VCF_LIST_SEP)))
        
        raise ValueError("VCF record is missing INFO AC attribute")
    
    @staticmethod
    def parse_an(info_list: list):
        """
        Parse allele number from VCF INFO column
        :param info_list: list of VCF INFO key=value strings
        :return: allele number
        """
        for item in info_list:
            # find total allele number
            if item[:3] == 'AN=':
                # AN=value
                return int(item[3:])
        
        raise ValueError("VCF record is missing INFO AN attribute")
    
    def compact_base_count(self, count_list: list):
        """
        Compact base counts so each of them is lower than MAX_BASE_COUNT
        :param count_list:
        :return: compacted ac_list
        """
        max_count = max(count_list)
        if max_count > self.MAX_ALLELE_COUNT:
            ratio = self.MAX_ALLELE_COUNT / float(max_count)
            new_count_list = np.round(np.array(count_list) * ratio).astype(int).tolist()
            if self.verbose:
                print("compacting base counts from %s to %s" % (count_list, new_count_list))
            return new_count_list
        else:
            return count_list
    
    def snv_count_map(self, allele_list: list, count_list: list):
        """
        :return: SNV count map
        """
        assert len(allele_list) == len(count_list)
        
        count_map = dict(zip(allele_list, self.compact_base_count(count_list)))
        
        for base in BASES:
            if base not in count_map:
                # add missing base
                count_map[base] = 0
        
        return count_map
    
    def indel_count_map(self, allele_list: list, count_list: list):
        """
        :return: sorted INDEL count map as list of tuples
        """
        assert len(allele_list) == len(count_list)
        return sorted(
            zip(self.compact_base_count(count_list), allele_list),
            key=lambda x: (x[0], x[1]),
            reverse=True
        )
    
    @classmethod
    def read_header(cls, vac_file):
        return struct.unpack('<II', vac_file.read(cls.HEADER_SIZE))
    
    @classmethod
    def write_header(cls, vac_file, snv_count: int, indel_count: int):
        vac_file.write(struct.pack('<II', snv_count, indel_count))
    
    @classmethod
    def read_snv_record(cls, vac_file):
        """
        :param vac_file:
        :return: index, allele count tuple
        """
        byte_str = vac_file.read(cls.SNV_RECORD_SIZE)
        if len(byte_str) == 0:
            raise EOFError()
        index, a_ac, t_ac, c_ac, g_ac = struct.unpack(cls.SNV_FORMAT, byte_str)
        return index, (a_ac, t_ac, c_ac, g_ac)
    
    @classmethod
    def write_snv_record(cls, vac_file, index: int, ac_tuple: tuple):
        """
        :param vac_file:
        :param index: genomic index
        :param ac_tuple: counts of bases in format (A,T,G,C)
        :return:
        """
        # if sum(ac_tuple) > 0:
        vac_file.write(struct.pack(cls.SNV_FORMAT, index, *ac_tuple))
    
    @classmethod
    def read_indel_record(cls, vac_file):
        """
        :param vac_file:
        :return:
        """
        byte_str = vac_file.read(cls.INT_SIZE)
        if len(byte_str) == 0:
            raise EOFError()
        
        index = struct.unpack('<I', byte_str)[0]
        length = int.from_bytes(vac_file.read(1), byteorder='little')
        
        counts = [0] * length
        seqs = [''] * length
        for i in range(length):
            count, base_length = struct.unpack('<HH', vac_file.read(cls.SHORT_SIZE * 2))
            seq_byte_size = math.ceil(base_length / 4)
            seq = bytes2seq(vac_file.read(seq_byte_size), base_length)
            counts[i] = count
            seqs[i] = seq
        
        return index, counts, seqs
    
    # TODO parameters: index, counts, seqs
    @staticmethod
    def write_indel_record(indel_file, index: int, indel_map: list):
        """
        :param indel_file:
        :param index: genomic index
        :param indel_map: sorted list of tuples [(SEQ_A,COUNT_A), (SEQ_B:COUNT_B), ...]
        """
        assert len(indel_map) < 256
        # single allele can be replaced only by itself
        # if len(indel_map) > 1 and sum(indel_map.values()) > 0:
        # noinspection PyTypeChecker
        record = struct.pack('<I', index) + bytes([len(indel_map)])
        
        # sort indel_map by count and sequence
        # for sequence, count in sorted(indel_map.items(), key=lambda x: (x[1], x[0])):
        for allele_count, sequence in indel_map:
            record += struct.pack('<HH', allele_count, len(sequence))
            record += seq2bytes(sequence)
        
        indel_file.write(record)
    
    def vcf2vac(self, vcf_file, vac_file):
        """
        Converts VCF to binary VAC file.
        Can preserve compatibility between different BAM and VCF reference genomes
        by stripping chr prefix via FastaIndex.
        :param vcf_file: input VCF file
        :param vac_file: output VAC file
        """
        variant_cnt = 0
        snv_count = 0
        indel_count = 0
        
        snv_filename = vac_file.name + self.SNV_TEMP_EXT
        indel_filename = vac_file.name + self.INDEL_TEMP_EXT
        
        with open(snv_filename, 'wb') as snv_file, \
                open(indel_filename, 'wb') as indel_file:
            for line in vcf_file:
                if line[0] == "#":
                    # skip header line
                    continue
                variant_cnt += 1
                
                # split until last used column only
                data = line.split(self.VCF_COL_SEP, maxsplit=8)
                allele_list = [data[self.VCF_REF_ID]] + data[self.VCF_ALT_ID].split(self.VCF_LIST_SEP)
                
                is_snv = self.is_snv(allele_list)
                if is_snv or self.is_indel(allele_list):
                    chrom = data[self.VCF_CHROM_ID]
                    
                    pos = int(data[self.VCF_POS_ID]) - 1  # vcf has 1-based index, convert it to 0-based index
                    index = self.fai.pos2index(chrom, pos)
                    
                    count_list = self.parse_allele_counts(data[self.VCF_INFO_ID])
                    if is_snv:
                        snv_count += 1
                        # allele count map
                        count_map = self.snv_count_map(allele_list, count_list)
                        count_tuple = (count_map['A'], count_map['T'], count_map['G'], count_map['C'])
                        self.write_snv_record(snv_file, index, count_tuple)
                    
                    else:  # is indel
                        indel_count += 1
                        indel_map = self.indel_count_map(allele_list, count_list)
                        self.write_indel_record(indel_file, index, indel_map)
                
                if self.verbose and variant_cnt % 10000 == 0:
                    print("variant %d" % variant_cnt)
        
        self.write_header(vac_file, snv_count, indel_count)
        
        if self.verbose:
            print("total variants %d" % variant_cnt)
            print("total SNVs %d" % snv_count)
        
        self.__final_merge(vac_file, snv_filename, indel_filename)
    
    def __final_merge(self, vac_file, snv_filename: str, indel_filename: str):
        if self.verbose:
            print("merging temporary files")
        
        with open(snv_filename, 'rb') as snv_file:
            # for chunk in iter(lambda: snv_file.read(4096), b""):
            #     vac_file.write(chunk)
            
            vac_file.write(snv_file.read())
        
        with open(indel_filename, 'rb') as indel_file:
            # for chunk in iter(lambda: indel_file.read(4096), b""):
            #     vac_file.write(chunk)
            
            vac_file.write(indel_file.read())
        
        os.remove(snv_filename)
        os.remove(indel_filename)
    
    def parse_allele_counts(self, info_value: str):
        """
        :param info_value: VCF INFO value
        :return: list of allele counts
        [ref_count, allele_1_count, allele_2_count, ...]
        """
        info_list = info_value.split(self.VCF_INFO_SEP)
        # alternative allele count list
        info_ac = self.parse_ac(info_list)
        # total allele number
        info_an = self.parse_an(info_list)
        # sanity check
        ref_count = info_an - sum(info_ac)
        assert ref_count >= 0
        return [ref_count] + info_ac
    
    @classmethod
    def text2vac(cls, text_filepath, vac_filepath):
        """
        Text genomic indices are 1-based. After conversion to binary they change to 1-based.
        :param text_filepath:
        :param vac_filepath:
        :return:
        """
        with open(vac_filepath, 'wb') as vac_file, \
                open(text_filepath, 'rt') as text_file:
            
            snv_count = int(text_file.readline().rstrip())
            vac_file.write(struct.pack('<I', snv_count))
            
            indel_count = int(text_file.readline().rstrip())
            vac_file.write(struct.pack('<I', indel_count))
            
            # write SNVs
            for i in range(snv_count):
                index_str, ac_str = text_file.readline().rstrip().split(cls.VCF_COL_SEP)
                ac_tuple = tuple(map(int, ac_str.split(',')))
                cls.write_snv_record(vac_file, int(index_str) - 1, ac_tuple)
            
            # write INDELs
            for line in text_file:
                index_str, indel_str = line.rstrip().split('\t')
                indel_map = []
                for segment in indel_str.split(','):
                    allele_count, sequence = segment.split(':')
                    indel_map.append((allele_count, sequence))
                
                cls.write_indel_record(vac_file, int(index_str) - 1, indel_map)
    
    @classmethod
    def vac2text(cls, vac_filepath, text_filepath):
        """
        Binary genomic indices are 0-based. After conversion to text format they change to 1-based.
        :param vac_filepath:
        :param text_filepath:
        :return:
        """
        with open(vac_filepath, 'rb') as vac_file, \
                open(text_filepath, 'wt') as text_file:
            snv_count, indel_count = cls.read_header(vac_file)
            
            text_file.write('%d\n' % snv_count)
            text_file.write('%d\n' % indel_count)
            for i in range(snv_count):
                index, count_list = cls.read_snv_record(vac_file)
                ac_string = ','.join(map(str, count_list))
                text_file.write('%d\t%s\n' % (index + 1, ac_string))
            
            for i in range(indel_count):
                index, counts, seqs = cls.read_indel_record(vac_file)
                
                text_file.write('%d\t' % (index + 1))
                indels = []
                for count, seq in zip(counts, seqs):
                    indels.append('%d:%s' % (count, seq))
                
                text_file.write(','.join(indels))
                text_file.write('\n')
