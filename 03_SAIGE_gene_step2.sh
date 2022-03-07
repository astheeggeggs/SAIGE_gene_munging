#!/usr/bin/env bash
#$ -N SAIGE_gene_step2
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe@@short.hge
#$ -t 23-23

module purge
module load Anaconda3/2020.07
source activate /well/lindgren/users/mmq446/conda/skylake/envs/saige

source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
chr=$(get_chr ${SGE_TASK_ID})

# Find and replace @ with the current chromosome
vcfFile=${vcfFile//[@]/${chr}}
vcfFileIndex=${vcfFileIndex//[@]/${chr}}
groupFile=${groupFile//[@]/${chr}}

vcf_chr='chr'$(get_chr ${SGE_TASK_ID})

Rscript SAIGE_gene_step2_wrapper.r \
    --vcfFile ${vcfFile} \
    --vcfFileIndex ${vcfFileIndex} \
	--phenotypeFolder ${phenotypeFolder} \
	--phenotype ${phenotype} \
    --groupFile ${groupFile} \
    --chr ${vcf_chr} \
    --minMAF ${minMAF} \
    --maxMAFforGroupTest ${maxMAFforGroupTest}
