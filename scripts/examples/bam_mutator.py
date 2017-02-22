import io
import random

import varlock as vrl

bam_filename = 'examples/resources/sample.bam'
mut_bam_filename = 'examples/resources/sample.mut.bam'
vac_filename = 'examples/resources/sample.vac'

# bam_filename = 'resources/full_chr22.bam'
# vac_filename = 'resources/full_chr22.vac'
# mut_bam_filename = 'resources/full_chr22.mut.bam'

mut = vrl.BamMutator(rnd=random.Random(0), verbose=True)
with io.BytesIO() as diff_file:
    mut.mutate(bam_filename, vac_filename, mut_bam_filename, diff_file)

print(mut.all_stats())

# python3 -m cProfile -s tottime examples.bam_mutator
