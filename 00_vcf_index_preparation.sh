#!/usr/bin/env bash
#$ -N tabix_preparation
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe
#$ -t 1-22

source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
chr=$(get_chr ${SGE_TASK_ID})

vcfFile="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
/users/gms/whv244/htslib_1.11/bin/tabix -p vcf -C ${vcfFile}
