varlock_dir='/data/projects/varlock'
pca_dir="${varlock_dir}/pca/gnomad3nfe_onepanel"

mkdir -p "${pca_dir}/plink"

# https://www.biostars.org/p/335605/
plink --bcf "${pca_dir}/merged_samples_snp.bcf" \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x hg38 no-fail \
      --make-bed \
      --out "${pca_dir}/plink/merged_samples_snp";

# 8000 is tool's upper limit (maximum number of PCs is equal to lesser from number of samples or features)
plink --bfile "${pca_dir}/plink/merged_samples_snp" --pca 8000 "header" "tabs" "var-wts" \
  --out "${pca_dir}/merged_samples_snp"