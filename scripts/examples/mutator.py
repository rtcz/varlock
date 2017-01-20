import os
import random

from varlock import *

# python3 /data/projects/varlock/scripts/mutate.py
# python3 -m cProfile -s tottime /data/projects/varlock/scripts/varlock/vac.py

RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')
FAI_FILENAME = os.path.join(RESOURCES_DIR, "hg19.fa.fai")
IN_BAM_FILENAME = os.path.join(RESOURCES_DIR, "sample.bam")
IN_VAC_FILENAME = os.path.join(RESOURCES_DIR, "sample.vac")
OUT_BAM_FILENAME = os.path.join(RESOURCES_DIR, "sample.mut.bam")
OUT_DIFF_FILENAME = os.path.join(RESOURCES_DIR, "sample.diff")

# EOF BAM
# written alignments 517362
# unmapped alignments 0
# overlapping alignments 495465
# max snv alignments 667
# total snvs 1059517
# read snvs 1059417, 100 omitted
# mutations 2373839
# ns mutations 53775

# time 8.613

# TODO ??? add option for choosing GRCh37 vs hg19 BAM file

mutator = Mutator(parse_fai(FAI_FILENAME), rnd=random.Random(0), verbose=True)
mutator.mutate(
    in_vac_filename=IN_VAC_FILENAME,
    in_bam_filename=IN_BAM_FILENAME,
    out_bam_filename=OUT_BAM_FILENAME,
    out_diff_filename=OUT_DIFF_FILENAME
)
