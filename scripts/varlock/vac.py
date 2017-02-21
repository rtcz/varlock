import struct

import numpy as np

from .common import *


class Vac:
    """
    Class for handling VAC file. VAC stands for Variant Allele Count.
    VAC is binary file, where each record represent one SNV from VCF.
    
    Vac record:
    index, 4B, absolute position of SNV in genome
    A allele count, 2B
    T allele count, 2B
    C allele count, 2B
    G allele count, 2B
    """
    VCF_CHROM_ID = 0
    VCF_POS_ID = 1
    VCF_REF_ID = 3
    VCF_ALT_ID = 4
    VCF_INFO_ID = 7
    
    COL_SEP = "\t"
    LIST_SEP = ","
    INFO_SEP = ";"
    
    MAX_BASE_COUNT = 2 ** 16 - 1
    """
    Max value of 2 bytes
    """
    
    STRUCT_FORMAT = "<IHHHH"  # int, short, short, short
    STRUCT_LENGTH = 12  # bytes
    
    def __init__(self, fai, verbose=False):
        """
        :param fai: FastaIndex
        :param verbose:
        """
        self.fai = fai
        self.current_pos = None
        self.current_chrom = None
        self.verbose = verbose
    
    @staticmethod
    def is_snp(ref, alt):
        """
        :param ref: reference string
        :param alt: array of alleles
        :return: True if variant is SNP
        """
        return len(ref) == 1 and all(len(value) == 1 for value in alt)
    
    @classmethod
    def parse_ac(cls, info_list):
        """
        Parse allele count from VCF INFO column.
        :param info_list: list of VCF INFO key=value strings
        :return: list of allele counts
        """
        for item in info_list:
            # find allele count item
            if item[:3] == 'AC=':
                # AC=value1,value2
                return list(map(int, item[3:].split(cls.LIST_SEP)))
        
        raise ValueError("VCF record is missing INFO AC attribute")
    
    @staticmethod
    def parse_an(info_list):
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
    
    def compact_base_count(self, ac_list):
        """
        Compact base counts so each of them is lower than MAX_BASE_COUNT
        :param ac_list:
        :return: compacted ac_list
        """
        max_count = max(ac_list)
        if max_count > self.MAX_BASE_COUNT:
            ratio = float(max_count) / self.MAX_BASE_COUNT
            new_ac_list = np.round(np.array(ac_list) * ratio).astype(int)
            if self.verbose:
                print("compacting base counts from %s to %s" % (ac_list, new_ac_list))
            return new_ac_list
        else:
            return ac_list
    
    def __snv2vac(self, ref_name, pos, ref, alt_list, info_list, vac_file):
        """
        :param ref_name: reference sequence name
        :param pos: 0-based position
        :param ref: reference base
        :param alt_list: list of alternative alleles
        :param info_list: list of VCF INFO key=value strings
        :return: byte string
        """
        # alternative allele count list
        info_ac = self.parse_ac(info_list)
        # total allele count
        info_an = self.parse_an(info_list)
        # reference allele count
        ref_ac = info_an - sum(info_ac)
        
        base_list = [ref] + alt_list
        count_list = self.compact_base_count([ref_ac] + info_ac)
        
        ref_name = ref_name if self.fai.keep_chr else strip_chr(ref_name)
        
        index = self.fai.pos2index(ref_name, pos)
        ac_map = dict(zip(base_list, count_list))
        
        for base in BASES:
            if base not in ac_map:
                ac_map[base] = 0
        
        self.write_record(vac_file, index, (ac_map['A'], ac_map['T'], ac_map['G'], ac_map['C']))
    
    @classmethod
    def write_record(cls, vac_file, index, ac_tuple):
        vac_file.write(struct.pack(cls.STRUCT_FORMAT, index, *ac_tuple))
    
    def vcf2vac(self, vcf_file, vac_file):
        """
        Converts VCF to binary VAC file.
        Preserves compatibility between different BAM and VCF reference genomes by stripping chr prefix.
        :param vcf_file: input VCF file
        :param vac_file: output VAC file
        """
        variant_counter = 0
        snp_counter = 0
        
        for line in vcf_file:
            if line[0] == "#":
                # skip header line
                continue
            
            variant_counter += 1
            
            # split until last used column only
            data = line.split(self.COL_SEP, maxsplit=8)
            ref = data[self.VCF_REF_ID]
            alt_list = data[self.VCF_ALT_ID].split(self.LIST_SEP)
            
            if self.is_snp(ref, alt_list):
                snp_counter += 1
                chrom = data[self.VCF_CHROM_ID]
                pos = int(data[self.VCF_POS_ID]) - 1  # vcf has 1-based index, convert to 0-based index
                info_list = data[self.VCF_INFO_ID].split(self.INFO_SEP)
                self.__snv2vac(
                    ref_name=chrom,
                    pos=pos,
                    ref=ref,
                    alt_list=alt_list,
                    info_list=info_list,
                    vac_file=vac_file
                )
            
            if self.verbose and variant_counter % 10000 == 0:
                print("variant %d" % variant_counter)
        
        if self.verbose:
            print("total variants %d" % variant_counter)
            print("total SNVs %d" % snp_counter)
    
    @classmethod
    def read_record(cls, vac_file):
        byte_string = vac_file.read(cls.STRUCT_LENGTH)
        if len(byte_string) == 0:
            raise EOFError()
        index, a_ac, t_ac, c_ac, g_ac = struct.unpack(cls.STRUCT_FORMAT, byte_string)
        return index, (a_ac, t_ac, c_ac, g_ac)
    
    @classmethod
    def text2vac(cls, text_filepath, vac_filepath):
        with open(vac_filepath, 'wb') as vac_file, \
                open(text_filepath, 'rt') as text_file:
            for line in text_file:
                index, ac_string = line.rstrip().split(cls.COL_SEP)
                ac_tuple = tuple(map(int, ac_string.split(',')))
                cls.write_record(vac_file, int(index), ac_tuple)
    
    @classmethod
    def vac2text(cls, vac_filepath, text_filepath):
        with open(vac_filepath, 'rb') as vac_file, \
                open(text_filepath, 'wt') as text_file:
            while True:
                try:
                    index, ac_tuple = cls.read_record(vac_file)
                    ac_string = ','.join(map(str, ac_tuple))
                    text_file.write('%d\t%s\n' % (index, ac_string))
                except EOFError:
                    break
