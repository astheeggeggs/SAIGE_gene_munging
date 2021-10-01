import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='22')
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche
CHR = str(args.chr)

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

# Inputs
MT  = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '.mt'
INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR +'.keep.variant.ht'

PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'

# Outputs
INITIAL_SAMPLE_QC_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_chr' + CHR + '_initial_sample_qc.tsv.bgz'

variants_to_filter = hl.read_table(INITIAL_VARIANT_LIST)
sample_annotations = hl.read_table(PHENOTYPES_TABLE)

select_rows = [
	'age','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6',
	'ukbb.centre11001', 'ukbb.centre11002', 'ukbb.centre11003', 'ukbb.centre11004', 'ukbb.centre11005',
	'ukbb.centre11006', 'ukbb.centre11007', 'ukbb.centre11008', 'ukbb.centre11009', 'ukbb.centre11010',
	'ukbb.centre11011', 'ukbb.centre11012', 'ukbb.centre11013', 'ukbb.centre11014', 'ukbb.centre11016',
	'ukbb.centre11017', 'ukbb.centre11018', 'ukbb.centre11020', 'ukbb.centre11021', 'ukbb.centre11022',
	'ukbb.centre11023', 'sex2', 'sequencing.batch2'
	]

mt = hl.read_matrix_table(MT)
mt = mt.filter_rows(hl.is_defined(variants_to_filter[mt.row_key]))
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s].select(*select_rows))
mt = mt.annotate_cols(gq = hl.agg.stats(mt.GQ), dp = hl.agg.stats(mt.DP))
mt = hl.sample_qc(mt, name='sample_qc')

mt.count_rows()
mt.cols().select('phenotype', 'sample_qc', 'gq', 'dp').flatten().export(output=INITIAL_SAMPLE_QC_FILE)
