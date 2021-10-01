#!/usr/bin/env bash
#$ -N hail_shell
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/hail.log
#$ -e /well/lindgren/dpalmer/logs/hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q long.qc
#$ -t 1-22

set -o errexit
set -o nounset

module purge
source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
source /well/lindgren/dpalmer/ukb_utils/bash/hail_utils.sh

spark_dir="/well/lindgren/dpalmer/data/tmp/spark"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer"
set_up_hail

chr=$(get_chr ${SGE_TASK_ID})

python 00_WES_data_preparation.py --chr ${chr}
print_update "Finished running Hail for chr${chr}" "${SECONDS}"