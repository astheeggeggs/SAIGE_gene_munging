#!/usr/bin/env bash
#$ -N hail_shell
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/hail.log
#$ -e /well/lindgren/dpalmer/logs/hail.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 10
#$ -q short.qe
#$ -t 1-24

set -o errexit
set -o nounset

module purge
source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh
source /well/lindgren/dpalmer/ukb_utils/bash/hail_utils.sh

spark_dir="/well/lindgren/dpalmer/data/tmp/spark5"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer"
set_up_hail

set_up_vep() {
  module load EnsEMBLCoreAPI/96.0-r20190601-foss-2019a-Perl-5.28.1 # required for LOFTEE
  module load VEP/95.0-foss-2018b-Perl-5.28.0 # required FOR VEP (NOTE: this steps throws some errors since the above module is already loaded. It works nonetheless.)
  module load samtools/1.8-gcc5.4.0 # required for LOFTEE 
  export PERL5LIB=$PERL5LIB:/well/lindgren/flassen/software/VEP/plugins_grch38/
}

set_up_vep

chr=$(get_chr ${SGE_TASK_ID})

python 05_annotate_variants_vep.py --chr ${chr}
print_update "Finished running Hail for chr${chr}" "${SECONDS}"
