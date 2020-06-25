import gzip

import pysam

from src.bam import open_bam
from src.vac import Vac
from src.fasta_index import FastaIndex
from io import BufferedReader


bam_filename = 'resources/chr20_10m.bam'
vcf_filename = 'resources/chr20_10m.vcf.gz'
out_filename = 'resources/chr20_10m.vac'

# bam_filename = 'resources/full_chr22.bam'
# vcf_filename = 'resources/full_chr22.vcf.gz'
# out_filenam = 'resources/full_chr22.vac'

# TODO test two variants of temporary files merging

with pysam.AlignmentFile(bam_filename, 'rb') as sam_file, \
        gzip.open(vcf_filename, 'rt') as vcf_file, \
        open(out_filename, 'wb') as out_vac_file:
    
    vac = Vac(FastaIndex(sam_file.header), verbose=True)
    vac.vcf2vac(vcf_file, out_vac_file)

# python3 -m cProfile -s tottime examples.vac
# full_chr22 cca 120s

# STATS full_chr22
# chr22 length 51304566
# chr22 vcf total lines 1103800
# chr22 vcf data lines 1103800 - 253 = 1103547
# chr22 total SNPs 1059517
# HG putative SNP count 66084847
# cca 793MB
