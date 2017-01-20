import gzip
import struct

import numpy as np

from .common import *


class Vac:
    """
    Class for handling VAC file. VAC stands for Variant Allele Count.
    VAC is binary file, where each record represent one SNV from VCF.
    SNV record is 12 bytes long and has following format:
    index 4B -> absolute position of SNV in genome
    A allele count 2B
    T allele count 2B
    C allele count 2B
    G allele count 2B
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
    Max value of 2B
    """
    
    def __init__(self, fai_list, verbose=False):
        """
        :param fai_list:
        :param verbose:
        """
        # create dict for fast access
        self.fai_dict = fai_list2dict(fai_list)
        
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
    
    def __snv2vac(self, chrom, pos, ref, alt_list, info_list):
        """
        :param chrom: reference sequence name
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
        
        index = pos2index(chrom, pos, self.fai_dict)
        ac_map = dict(zip(base_list, count_list))
        return self.__pack2bytes(index, ac_map)
    
    @staticmethod
    def __pack2bytes(index, ac_map):
        # index is an integer
        record = struct.pack('<I', index)
        for base in BASES:
            if base in ac_map:
                # reference or alternative base
                # base allele count is a short integer
                record += struct.pack('<H', ac_map[base])
            else:
                # base is not present
                record += bytes(2)
        
        return record
    
    def vcf2vac(self, vcf_filename, out_filename):
        """
        Converts VCF to binary VAC file
        :param vcf_filename: input VCF filename, file must be gzipped
        :param out_filename: output VAC filename
        """
        variant_counter = 0
        snp_counter = 0
        with gzip.open(vcf_filename, "rt") as vcf_file, \
                open(out_filename, "wb") as ac_file:
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
                    byte_string = self.__snv2vac(
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt_list=alt_list,
                        info_list=info_list
                    )
                    ac_file.write(byte_string)
                
                if self.verbose and variant_counter % 10000 == 0:
                    print("variant %d" % variant_counter)
        
        if self.verbose:
            print("total variants %d" % variant_counter)
            print("total SNVs %d" % snp_counter)
    
    @classmethod
    def text2vac(cls, text_filename, vac_filename):
        with open(text_filename, "r") as text_file, \
                open(vac_filename, "wb") as vac_file:
            for line in text_file:
                data = list(map(int, line.strip().split(cls.COL_SEP)))
                vac_file.write(struct.pack('<IHHHH', *data))
    
    @staticmethod
    def vac2list(vac_filename):
        data = []
        with open(vac_filename, "rb") as vac_file:
            while True:
                byte_string = vac_file.read(12)
                if len(byte_string) == 0:
                    break
                data.append(struct.unpack('<IHHHH', byte_string))
        
        return data
