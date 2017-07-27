import gzip

import pysam

import varlock as vrl

bam_filename = 'examples/resources/sample.bam'
vcf_filename = 'examples/resources/sample.vcf.gz'
out_filenam = 'examples/resources/sample.vac'

# bam_filename = 'resources/full_chr22.bam'
# vcf_filename = 'resources/full_chr22.vcf.gz'
# out_filenam = 'resources/full_chr22.vac'

# TODO test two variants of temporary files merging

with pysam.AlignmentFile(bam_filename, 'rb') as sam_file, \
        gzip.open(vcf_filename, 'rt') as vcf_file, \
        open(out_filenam, 'wb') as out_vac_file:
    vac = vrl.Vac(vrl.FastaIndex(sam_file, keep_chr=False), verbose=True)
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
