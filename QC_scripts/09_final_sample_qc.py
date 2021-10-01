import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, required=True)
args = parser.parse_args()

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'

PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'
IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/05_imputesex.ht'

INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
SEXCHECK_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/05_sexcheck.remove.sample_list'

INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR + '.keep.variant.ht'
FINAL_VARIANT_LIST <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list')


sample_annotations = hl.read_table(PHENOTYPES_TABLE)
impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE)
ht_initial_variants = hl.read_table(INITIAL_VARIANT_LIST)
ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_sexcheck_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')

SAMPLE_BEFORE_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.before.samples.tsv'
SAMPLE_AFTER_QC_FILE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.after.samples.tsv'

ht_final_variants = hl.import_table(FINAL_VARIANT_LIST,
	types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)


mt_before = hl.read_matrix_table(MT_HARDCALLS)
mt_before = mt_before.annotate_cols(phenotype = sample_annotations[mt_before.s])
mt_before = mt_before.annotate_cols(imputesex = impute_sex_annotations[mt_before.col_key])

# Filter down to the collection of Europeans and samples that were outliers in the initial sample QC.
mt_before = mt_before.filter_rows(hl.is_defined(ht_initial_variants[mt_before.row_key]))
# European list, or White-British list.
mt_before = mt_before.filter_cols(mt_before.genetic.eur  == 1)
mt_before = mt_before.filter_cols(hl.is_defined(ht_initial_samples[mt_before.col_key]))
mt_before = mt_before.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt_before.col_key]))
# Need to decide whether to remove samples with excess ultra-rare variants

mt_before = hl.variant_qc(mt_before, name = 'qc')

mt_before = mt_before.annotate_rows(
	qc = mt_before.qc.annotate(AC=mt_before.qc.AC[1],
	AF = mt_before.qc.AF[1],
	homozygote_count = mt_before.qc.homozygote_count[1]))

mt_before = mt_before.filter_rows((mt_before.qc.AF > 0) & (mt_before.qc.AF < 1))
mt_before = hl.sample_qc(mt_before)

n = mt_before.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt_after = mt_before.filter_rows(hl.is_defined(ht_final_variants[mt_before.row_key]))
mt_after = hl.sample_qc(mt_after)

n = mt_after.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

mt_before.cols().select("phenotype", "sample_qc").flatten().export(SAMPLE_BEFORE_QC_FILE)
mt_after.cols().select("phenotype", "sample_qc").flatten().export(SAMPLE_AFTER_QC_FILE)
