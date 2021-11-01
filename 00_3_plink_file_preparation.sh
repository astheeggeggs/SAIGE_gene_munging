#!/usr/bin/env bash
#$ -N plink_prepation
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe
#$ -pe shmem 10

module purge
out="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_1_22_including_rare_saige_input"

# merge together rare variant plink files
mergelist="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/merge_list.txt"
chr=1
echo "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chr${chr}_saige_input" > ${mergelist}
for chr in {2..22}; do
    echo "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chr${chr}_saige_input" >> ${mergelist}
done

common_bed="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211026_long_ukb_wes_200k_sparse_autosomes.bed"
common_bim="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211026_long_ukb_wes_200k_sparse_autosomes.bim"
common_fam="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_imp_eur_chr1_22_sparse_markers.fam"

# Merge in rare variants and remove SNPs with high missingness
./../plink --bim ${common_bim} --bed ${common_bed} --fam ${common_fam} --merge-list ${mergelist} --make-bed --out ${out}

# for chr in {1..22}; do
#     rm /well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chr${chr}_saige_input.*
# done
