#!/usr/bin/env bash
#$ -N tabix_preparation
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe@@short.hge
#$ -t 1-23

source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
CHR=$(get_chr ${SGE_TASK_ID})
TRANCHE='200k'

module purge
module load BCFtools

# vcfFile="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
vcfFile="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr${CHR}.vcf.bgz"
# /users/gms/whv244/htslib_1.11/bin/tabix -p vcf -C ${vcfFile}
tabix -p vcf -C ${vcfFile}
