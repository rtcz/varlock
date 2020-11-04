#!/bin/bash

project_dir=/data/projects/varlock
mkdir -p ${project_dir}/log

# for pipenv shell
cd ${project_dir}/scripts/varlock

for bam_file in /data/projects/varlock/mapping/grch38_decoy_alt/masked_gnomad3nfe_onepanel/one_*.bam
do
  echo "$bam_file"
  bam_filename=$(basename "$bam_file" ".bam")
  # pipenv run python varlock.py mask \
  # $(pipenv --venv)/bin/python varlock.py mask \
  cmd=" \
    pipenv run python ${project_dir}/scripts/varlock/varlock.py unmask \
      --key ${project_dir}/masking/gnomad3nfe_onepanel/id_rsa \
      --ver_key ${project_dir}/masking/gnomad3nfe_onepanel/id_rsa.pub \
      --bam ${bam_file} \
      --diff ${project_dir}/masking/gnomad3nfe_onepanel/${bam_filename}.bdiff \
      --out_bam ${project_dir}/mapping/grch38_decoy_alt/unmasked_gnomad3nfe_onepanel/${bam_filename}.bam \
      --include_unmapped \
      --verbose \
  "
#  echo "$cmd"
#  break
  eval "$cmd" \
    1> "${project_dir}/log/unmask_${bam_filename}.out" \
    2> "${project_dir}/log/unmask_${bam_filename}.err" \
    &

#  echo "$cmd" | qsub \
#    -cwd \
#    -o "${project_dir}/log/mask_${bam_filename}.out" \
#    -e "${project_dir}/log/mask_${bam_filename}.err" \
#    -N "mask_${bam_filename}"
done
