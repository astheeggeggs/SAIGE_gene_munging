#!/usr/bin/env bash

# UKB genotype calls
ukb_bed="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr"
ukb_bim="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr"
ukb_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"

out="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_pruned_"

# Create the vcf fam
TRANCHE='200k'
pheno="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/${TRANCHE}/QC_phenotypes.tsv.gz"
vcf_samples="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/${TRANCHE}/UKBB_WES${TRANCHE}.txt"
zcat $pheno | awk '{print $1, $1}' | tail -n +2 > $vcf_samples

for i in `seq 1 22`;
do
      /well/lindgren/dpalmer/plink --bed ${ukb_bed}${i}_v2.bed --bim ${ukb_bim}${i}_v2.bim --fam ${ukb_fam} \
      --keep ${vcf_samples} --maf 0.05 --geno 0.02 \
      --indep 50 5 2 --out ${out}chr${i}
      /well/lindgren/dpalmer/plink --bed ${ukb_bed}${i}_v2.bed --bim ${ukb_bim}${i}_v2.bim --fam ${ukb_fam} \
      --extract ${out}chr${i}.prune.in --make-bed --out ${out}chr${i}
done 

echo ${out}chr2 > ~/tmp_plink.txt
for i in `seq 3 22`;
do
      echo ${out}chr${i} >> ~/tmp_plink.txt
done

/well/lindgren/dpalmer/plink --bfile ${out}chr1 --merge-list ~/tmp_plink.txt --make-bed --out ${out}autosomes

for i in `seq 2 22`;
do
      rm ${out}chr${i}.*
done

rm ~/tmp_plink.txt
