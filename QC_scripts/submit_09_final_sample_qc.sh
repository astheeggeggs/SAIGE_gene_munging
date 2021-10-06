#!/usr/bin/env bash
#$ -N hail_shell
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/hail.log
#$ -e /well/lindgren/dpalmer/logs/hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 24
#$ -q short.qe
#$ -t 1-24

set -o errexit
set -o nounset

module purge
source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
source /well/lindgren/dpalmer/ukb_utils/bash/hail_utils.sh

spark_dir="/well/lindgren/dpalmer/data/tmp/spark9"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer"
set_up_hail

chr=$(get_chr ${SGE_TASK_ID})

python 09_final_sample_qc.py --chr ${chr}
print_update "Finished running Hail for chr${chr}" "${SECONDS}"
