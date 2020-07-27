#!/bin/bash

project_dir=/data/projects/varlock
mkdir -p ${project_dir}/log

# for pipenv shell
cd ${project_dir}/scripts/varlock

for bam_file in /data/projects/exome/mapping/grch38_decoy_alt/deduplicated/one_*.bam
do
  echo "$bam_file"
  bam_filename=$(basename "$bam_file" ".bam")
  # pipenv run python varlock.py mask \
  # $(pipenv --venv)/bin/python varlock.py mask \
  cmd=" \
    pipenv run python ${project_dir}/scripts/varlock/varlock.py mask \
      --key ${project_dir}/mask/id_rsa \
      --pub_key ${project_dir}/mask/id_rsa.pub \
      --bam ${bam_file} \
      --vac ${project_dir}/vac/gnomad3nfe_pass_sorted_multisnp_onepanel.vac \
      --out_bam ${project_dir}/mapping/grch38_decoy_alt/deduplicated/${bam_filename}.bam \
      --out_diff ${project_dir}/mask/${bam_filename}.bdiff \
      --seed 0 \
      --verbose \
  "
  eval "$cmd" \
    1> "${project_dir}/log/mask_${bam_filename}.out" \
    2> "${project_dir}/log/mask_${bam_filename}.err" \
    &

#  echo "$cmd" | qsub \
#    -cwd \
#    -o "${project_dir}/log/mask_${bam_filename}.out" \
#    -e "${project_dir}/log/mask_${bam_filename}.err" \
#    -N "mask_${bam_filename}"
done
