#!/usr/bin/env bash
#$ -N SAIGE_gene_step1
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe@@short.hge

module purge
module load Anaconda3/2020.07
source activate /well/lindgren/users/mmq446/conda/skylake/envs/saige

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

Rscript SAIGE_gene_step1_wrapper.r --phenofile ${phenofile} --phenotype ${phenotype} \
	--covarColList ${covars} --sampleIDColinphenoFile ${ID_col} --traitType ${trait_type}
