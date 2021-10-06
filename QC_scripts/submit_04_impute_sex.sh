#!/usr/bin/env bash
#$ -N hail_shell
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/hail.log
#$ -e /well/lindgren/dpalmer/logs/hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 20
#$ -q short.qe
#$ -t 1-1

set -o errexit
set -o nounset

module purge
source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
source /well/lindgren/dpalmer/ukb_utils/bash/hail_utils.sh

spark_dir="/well/lindgren/dpalmer/data/tmp/spark4"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer"
set_up_hail

python 04_impute_sex.py
print_update "Finished running Hail" "${SECONDS}"
