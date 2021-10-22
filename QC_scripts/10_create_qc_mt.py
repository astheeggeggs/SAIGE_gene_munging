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
MT  = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '.mt'
FINAL_SAMPLE_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
FINAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

# Outputs:
QC_MT = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr' + CHR + '.mt'

ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0', delimiter=',', impute=True, types={'f0': hl.tstr})
ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

mt = hl.read_matrix_table(MT)
mt = mt.drop('qual', 'info', 'filters')

mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

mt.write(QC_MT, overwrite=True)
