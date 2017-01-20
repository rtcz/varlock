from varlock import *
import os

RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')
VCF_GZ_FILE = os.path.join(RESOURCES_DIR, "sample.vcf.gz")
CHR_FAI_FILE = os.path.join(RESOURCES_DIR, "hs37d5.fa.fai")
OUT_FILENAME = os.path.join(RESOURCES_DIR, "sample.vac")

vac = Vac(parse_fai(CHR_FAI_FILE), verbose=True)
vac.vcf2vac(vcf_filename=VCF_GZ_FILE, out_filename=OUT_FILENAME)

# CALLS
# python3 /data/projects/varlock/scripts/varlock/vac.py
# python3 -m cProfile -s tottime /data/projects/varlock/scripts/varlock/vac.py

# OUT
# chr22 length 51304566
# chr22 vcf total lines 1103800
# chr22 vcf data lines 1103800 - 253 = 1103547
# chr22 total SNPs 1059517
# HG putative SNP count 66084847
# cca 793MB

# INFO
# VCF reference:
# ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707
