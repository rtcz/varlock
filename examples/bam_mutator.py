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

# python3 -m cProfile -s tottime examples/bam_mutator.py
"""
total of 10000 alignments processed
{'alignment_count': 10000, 'snv_count': 9760, 'mut_count': 540, 'diff_count': 34, 'max_coverage': 652, 'overlapping_count': 4461}
         1523549 function calls (1507313 primitive calls) in 6.369 seconds

   Ordered by: internal time

   ncalls  tottime  percall  cumtime  percall _filename:lineno(function)
    23800    0.491    0.000    0.809    0.000 common.py:84(ref_pos2seq_pos)
    10000    0.425    0.000    0.428    0.000 calignmentfile.pyx:1361(write)
       14    0.358    0.026    0.358    0.026 {method 'read' of '_io.TextIOWrapper' objects}
     9760    0.292    0.000    1.700    0.000 mutator.py:511(create_mut_map)
    23800    0.283    0.000    0.283    0.000 calignedsegment.pyx:1449(get_aligned_pairs)
    29280    0.245    0.000    0.438    0.000 common.py:10(multi_random)
    29318    0.241    0.000    0.241    0.000 {built-in method array}
        1    0.231    0.231    4.871    4.871 mutator.py:270(mutate)
    29280    0.222    0.000    0.669    0.000 fromnumeric.py:43(_wrapit)
    29280    0.198    0.000    0.867    0.000 fromnumeric.py:917(argmax)
    34961    0.127    0.000    0.138    0.000 calignmentfile.pyx:672(get_reference_name)
     9760    0.120    0.000    2.084    0.000 mutator.py:365(__mutate_pos)
"""

"""
total of 517362 alignments processed
{'alignment_count': 517362, 'snv_count': 1059416, 'max_coverage': 770, 'mut_count': 53775, 'diff_count': 8852, 'overlapping_count': 495465}
         133488494 function calls (132964896 primitive calls) in 343.185 seconds

   Ordered by: internal time

   ncalls  tottime  percall  cumtime  percall _filename:lineno(function)
  2688971   39.119    0.000   63.576    0.000 common.py:84(ref_pos2seq_pos)
  2688971   21.604    0.000   21.604    0.000 calignedsegment.pyx:1449(get_aligned_pairs)
  1059416   21.457    0.000  125.544    0.000 mutator.py:511(create_mut_map)
  3178286   18.125    0.000   18.125    0.000 {built-in method array}
  3178248   17.618    0.000   31.765    0.000 common.py:10(multi_random)
        1   16.889   16.889  341.432  341.432 mutator.py:270(mutate)
  3178248   16.639    0.000   50.057    0.000 fromnumeric.py:43(_wrapit)
   517362   14.372    0.000   14.477    0.000 calignmentfile.pyx:1361(write)
  3178248   14.329    0.000   64.386    0.000 fromnumeric.py:917(argmax)
  3828681   10.228    0.000   11.106    0.000 calignmentfile.pyx:672(get_reference_name)
  1059416    9.203    0.000  155.791    0.000 mutator.py:365(__mutate_pos)
  2273799    8.814    0.000   27.776    0.000 mutator.py:37(__is_before_index)
  1059416    8.271    0.000   25.467    0.000 vac_iterator.py:14(__next__)
  3178248    7.519    0.000    7.519    0.000 {method 'argmax' of 'numpy.ndarray' objects}
  1059416    7.396    0.000    7.396    0.000 fasta_index.py:87(index2pos)
  3178248    6.160    0.000   24.285    0.000 numeric.py:414(asarray)
  1554880    5.485    0.000   14.934    0.000 mutator.py:47(__is_after_index)
  3828681    5.022    0.000    5.022    0.000 fasta_index.py:102(pos2index)
  1059416    4.984    0.000    6.000    0.000 mutator.py:526(<listcomp>)
  4747682    4.966    0.000    7.523    0.000 common.py:116(get_base)
  1059416    4.451    0.000    9.383    0.000 common.py:105(get_base_pileup)
  3828681    4.416    0.000   19.347    0.000 calignedsegment.pyx:844(__get__)
  1059416    3.929    0.000    4.994    0.000 mutator.py:493(count_bases)
  3828681    3.825    0.000   14.931    0.000 calignmentfile.pyx:1652(getrname)
  1059416    3.791    0.000    5.746    0.000 vac.py:180(read_snv_record)
  1059415    3.490    0.000   56.993    0.000 mutator.py:478(__set_seq_positions)
"""