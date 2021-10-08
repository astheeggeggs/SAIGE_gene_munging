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

# Inputs
MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'

# Outputs
INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR +'.keep.variant.ht'

mt = hl.read_matrix_table(MT_HARDCALLS)

mt = hl.variant_qc(mt, name='variant_qc')
mt = mt.filter_rows((mt.variant_qc.AF[0] > 0.0) & (mt.variant_qc.AF[0] < 1.0))

ht_rows_filter = mt.rows()
ht_rows_filter.select().write(INITIAL_VARIANT_LIST, overwrite=True)
