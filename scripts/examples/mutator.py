import os
import random

from varlock import *

# python3 /data/projects/varlock/scripts/mutate.py
# python3 -m cProfile -s tottime /data/projects/varlock/scripts/varlock/vac.py

RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')
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

mutator = Mutator(rnd=random.Random(0), verbose=True)
with pysam.AlignmentFile(IN_BAM_FILENAME, "rb") as in_bam_file, \
        open(IN_VAC_FILENAME, "rb") as in_vac_file, \
        pysam.AlignmentFile(OUT_BAM_FILENAME, "wb", template=in_bam_file) as out_bam_file, \
        open(OUT_DIFF_FILENAME, "wb") as out_diff_file:
        mutator.mutate(
            in_vac_file=in_vac_file,
            in_bam_file=in_bam_file,
            out_bam_file=out_bam_file,
            out_diff_file=out_diff_file
        )
