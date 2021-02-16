exome_dir='/data/projects/exome'
varlock_dir='/data/projects/varlock'
exome_vcf_dir="${exome_dir}/variant/grch38_decoy_alt-one/original"
varlock_vcf_dir="${varlock_dir}/variant/grch38_decoy_alt-one/masked_original"
pca_dir="${varlock_dir}/pca/gnomad3nfe_onepanel"

mkdir -p "${pca_dir}/tmp"

# 1; prepare samples for merge
while read -r sample;
do
    echo "sample: ${sample}"
#  echo "${exome_vcf_dir}/${sample}.vcf"
#  echo "${varlock_vcf_dir}/${sample}.vcf"

  bcftools view --apply-filters 'PASS' "${exome_vcf_dir}/${sample}.vcf" \
    --output-type 'b' \
    --output-file "${pca_dir}/tmp/${sample}.bcf"
  bcftools index -f "${pca_dir}/tmp/${sample}.bcf"

  bcftools view --apply-filters 'PASS' "${varlock_vcf_dir}/${sample}.vcf" \
    --output-type 'b' \
    --output-file "${pca_dir}/tmp/${sample}_masked.bcf"
  bcftools index -f "${pca_dir}/tmp/${sample}_masked.bcf"

  # add prefix to masked samples to make their name unique
  bcftools reheader --samples <(echo "${sample}_masked") "${pca_dir}/tmp/${sample}_masked.bcf" \
    --output "${pca_dir}/tmp/${sample}_masked.bcf.tmp"

  mv "${pca_dir}/tmp/${sample}_masked.bcf.tmp" "${pca_dir}/tmp/${sample}_masked.bcf"

  # reindex modified file
  bcftools index -f "${pca_dir}/tmp/${sample}_masked.bcf"

done < "$(dirname "$0")/sample_names.txt"

# 2; merge samples into single vcf
# shellcheck disable=SC2010
bcftools merge --filter-logic '+' --merge 'none' --output-type 'b' --missing-to-ref \
  --output "${pca_dir}/merged_samples.bcf" \
  --file-list <(ls ${pca_dir}/tmp/*.bcf | grep -v 'csi$')
bcftools index -f "${pca_dir}/merged_samples.bcf"


# 3; prepare vcf for plink, consider SNPs only
# see https://www.biostars.org/p/335605/
bcftools norm --multiallelics '-any' --check-ref 'w' \
  --fasta-ref "${varlock_dir}/reference/grch38_decoy_alt/grch38_decoy_alt.fa" \
  "${pca_dir}/merged_samples.bcf" | \
bcftools annotate --remove 'ID' --set-id +'%CHROM:%POS:%REF:%ALT' | \
bcftools norm --rm-dup 'both' | \
bcftools view -v 'snps' --output-type 'b' > "${pca_dir}/merged_samples_snp.bcf"

bcftools index "${pca_dir}/merged_samples_snp.bcf"



