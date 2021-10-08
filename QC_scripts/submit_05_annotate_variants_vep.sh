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

module load Anaconda3/2020.07
module load java/1.8.0_latest
source activate hail-new
_mem=$( get_hail_memory )
new_spark_dir=/well/lindgren/dpalmer/tmp/spark_test/
export PYSPARK_SUBMIT_ARGS="--conf spark.local.dir=${new_spark_dir} --conf spark.executor.heartbeatInterval=1000000 --conf spark.network.timeout=1000000  --driver-memory ${_mem}g --executor-memory ${_mem}g pyspark-shell"
export PYTHONPATH="${PYTHONPATH-}:/well/lindgren/dpalmer/ukb_utils/python:/well/lindgren/dpalmer:/well/lindgren/dpalmer/ukb_common/src"

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
