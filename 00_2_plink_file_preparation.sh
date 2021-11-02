#!/usr/bin/env bash
#$ -N plink_prepation
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe
#$ -t 1-22

module purge
module load PLINK/2.00a2.3_x86_64

source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
chr=$(get_chr ${SGE_TASK_ID})

# Use the fam from the kinship plink file (which went through all of our QC pipeline), filter to just those
eur_tsv="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_eur_chr1_22_sparse_markers.tsv"

# UKB genotype calls
ukb_bed="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr${chr}_v2.bed"
ukb_bim="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr${chr}_v2.bim"
ukb_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"

out="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chr${chr}_saige_input"

# Filter to rare variants and remove SNPs with high missingness
plink2 --bed ${ukb_bed} --bim ${ukb_bim} --fam ${ukb_fam} \
    --keep ${eur_tsv} --maf 0.0000000000001 --geno 0.02 --max-maf 0.0099 \
    --make-bed --out ${out}
