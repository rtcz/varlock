import sys
import numpy as np
import json
from bitarray import bitarray
import struct
import gzip

VCF_GZ_FILE = "in/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
CHR_FAI_FILE = "in/hs37d5.fa.fai"
CHR_LIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY']
BASES = ["A", "T", "G", "C"]


# VCF reference ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707


def parse_fai(fai_filename):
    data = []
    ref_keys = ["name", "length", "offset", "linebases"]
    with open(fai_filename, "r") as fai_file:
        for line in fai_file:
            raw_record = line.rstrip().split()
            ref_name = raw_record[0]
            ref_data = map(int, raw_record[1:])
            # noinspection PyTypeChecker
            ref = [ref_name] + list(ref_data)
            record_map = dict(zip(ref_keys, ref))
            data.append(record_map)
    return data


def load_json(file_name):
    with open(file_name, "r") as json_file:
        return json.load(json_file)


def dump_json(data, file_name):
    with open(file_name, "w") as json_file:
        json.dump(data, json_file, indent=4, sort_keys=True)


def pos2int(chr_id, chr_pos, fai_list):
    int_pos = 0
    for chr in fai_list:
        if chr["name"] == chr_id:
            return int(chr_pos) + int_pos
        else:
            int_pos += chr["length"]
    raise ValueError("chr %s not found" % chr_id)


def int2pos(pos, fai_list):
    for chr in fai_list:
        new_pos = pos - chr["length"]
        if new_pos < 0:
            return chr["name"], pos
        else:
            pos = new_pos
    raise ValueError("pos %d not found" % pos)


def strip_prefix(value, prefix):
    if value.startswith(prefix):
        return value[len(prefix):]
    else:
        return value


def strip_chr(value):
    return strip_prefix(value, 'chr')


def index2pos(index, fai_list):
    """
    :param index:
    :param fai_list:
    :return: 0-based position
    """
    ref_pos = index
    for ref in fai_list:
        new_ref_pos = ref_pos - ref["length"]
        if new_ref_pos < 0:
            return strip_chr(ref["name"]), ref_pos
        else:
            ref_pos = new_ref_pos
    raise ValueError("reference position for index %d not found" % index)


def pos2index(ref_name, ref_pos, fai_list):
    sum_ref_pos = 0
    for ref in fai_list:
        if strip_chr(ref["name"]) == strip_chr(ref_name):
            return int(ref_pos) + sum_ref_pos
        else:
            sum_ref_pos += ref["length"]
    raise ValueError("sequence name %s not found" % ref_name)


def is_snp(ref, alt):
    """
    :param ref: reference string
    :param alt: array of alleles
    :return: True if variant is SNP
    """
    is_snp = False
    if len(ref) == 1 and all(len(value) == 1 for value in alt):
        is_snp = True
    return is_snp


MAX_BASE_COUNT = 2 ** 16 - 1


def vcf2vac(vcf_filename, fai_filename, out):
    fai_list = parse_fai(fai_filename)
    line_counter = 0
    snp_counter = 0
    with gzip.open(vcf_filename, "r") as vcf_file, \
            open(out, "wb") as ac_file:
        for vcf_line in vcf_file:
            line_counter += 1
            
            if vcf_line[0] == "#":
                # skip header line
                continue
            
            data = vcf_line.split()
            ref = data[3]
            alt_list = data[4].split(",")
            
            if is_snp(ref, alt_list):
                snp_counter += 1
                
                chr = data[0]
                pos = int(data[1]) - 1  # vcf has 1-based index, convert to 0-based index
                info = data[7].split(";")
                # AC=value1,value2
                info_ac_list = map(int, info[0][3:].split(","))
                # AN=value
                info_an = int(info[2][3:])
                
                ref_ac = info_an - sum(info_ac_list)
                base_name_list = [ref] + alt_list
                base_count_list = [ref_ac] + info_ac_list
                max_count = max(base_count_list)
                
                if max_count > MAX_BASE_COUNT:
                    # TODO ??? raise ValueException
                    ratio = float(max_count) / MAX_BASE_COUNT
                    base_count_list = np.round(np.array(base_count_list) * ratio).astype(int)
                    print "decreasing counts for chr %s pos d%" % (chr, pos)
                
                int_pos = pos2int(chr, pos, fai_list)
                base_map = dict(zip(base_name_list, base_count_list))
                
                ac_file.write(struct.pack('<I', int_pos))
                
                for base in BASES:
                    if base in base_map:
                        # reference or alternative base
                        ac_file.write(struct.pack('<H', base_map[base]))
                    else:
                        # base is not present
                        ac_file.write("\x00\x00")
                        # if snp_counter == 10:
                        #     break
            if line_counter % 100000 == 0:
                print "%d" % line_counter
    
    print "total lines %d" % line_counter
    print "total SNPs %d" % snp_counter


def vcf_insight():
    with gzip.open(VCF_GZ_FILE, "r") as vcf_file:
        counter = 0
        for line in vcf_file:
            print line.rstrip()
            counter += 1
            if counter == 300:
                break
    
    print "--------------------------"
    
    with gzip.open(VCF_GZ_FILE, "r") as vcf_file:
        counter = 0
        for line in vcf_file:
            if line[0] == "#":
                continue
            variant_data = line.rstrip().split("\t")
            chr = variant_data[0]
            pos = variant_data[1]
            ref = variant_data[3]
            alt = variant_data[4]
            
            info_data = variant_data[7].split(";")
            alt_data = alt.split(",")
            
            if is_snp(ref, alt.split(",")):  # and len(alt_data) > 1:
                print chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + "\t".join(
                    [info_data[i] for i in [0, 1, 2, 3, 8]])
                counter += 1
                if counter == 100:
                    break


def parse_ac(ac_filename):
    with open(ac_filename, "rb") as ac_file:
        fai_list = parse_fai(CHR_FAI_FILE)
        pos_bytes = ac_file.read(4)
        while pos_bytes != "":
            int_pos = struct.unpack('<I', pos_bytes)[0]
            ac_list = struct.unpack('<HHHH', ac_file.read(8))
            base_map = dict(zip(BASES, ac_list))
            
            pos_bytes = ac_file.read(4)
            
            print int_pos
            print base_map
            print int2pos(int_pos, fai_list)


# vcf_insight()
# print parse_fai(CHR_FAI_FILE)

# parse_ac()

vcf2vac(VCF_GZ_FILE, CHR_FAI_FILE, "in/chr22.ac")
# parse_ac("out/chr22.ac")


# with open(VCF_GZ_FILE, "r") as vcf_file:
#     counter = 0
#     for line in vcf_file:
#         if line[0] == '#':
#             counter += 1
#         else:
#             print counter
#             break

# chr22 length 51304566
# chr22 vcf total lines 1103800
# chr22 vcf data lines 1103800 - 253 = 1103547
# chr22 total SNPs 1059517

# HG putative SNP count 66084847
# cca 793MB
