import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='20')
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche
CHR = str(args.chr)

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', 'GRCh38')

# Inputs:
QC_MT = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr' + CHR + '.mt'

# Outputs:
QC_HARDCALLS_MT = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr' + CHR + '.hardcalls.mt'

mt = hl.read_matrix_table(QC_MT)
mt.select_entries(mt.GT).repartition(512).write(QC_HARDCALLS_MT, overwrite=True)
