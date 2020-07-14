varlock_dir='/data/projects/varlock'
pca_dir="${varlock_dir}/pca"

# https://www.biostars.org/p/335605/
plink --bcf "${pca_dir}/merged_samples_snp.bcf" \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x hg38 no-fail \
      --make-bed \
      --out "${pca_dir}/plink/merged_samples_snp";

# TODO add missing genotype