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
MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'
PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'
IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/05_imputesex.ht'

INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
SEXCHECK_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_sexcheck.remove.sample_list'
URV_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/07_URV.remove.sample_list'

INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR + '.keep.variant.ht'

# Outputs:
VARIANT_QC_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/08_final_qc.variants_chr' + CHR + '.tsv.bgz'

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE)
ht_initial_variants = hl.read_table(INITIAL_VARIANT_LIST)
ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_sexcheck_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')
ht_urv_samples = hl.import_table(URV_LIST, no_header=True, key='f0')

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])

# Filter down to the collection of Non-Finnish (genetic_eur_no_fin_oct2021) Europeans and samples that were outliers in the initial sample QC.
mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))
mt = mt.filter_cols(
	(hl.is_defined(mt.phenotype.genetic_eur_no_fin_oct2021)) & (mt.phenotype.genetic_eur_no_fin_oct2021)
	)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_urv_samples[mt.col_key]))

mt = mt.annotate_cols(imputesex = impute_sex_annotations[mt.col_key])
mt = hl.variant_qc(mt, name='qc')

mt = mt.annotate_rows(
	qc=mt.qc.annotate(AC=mt.qc.AC[1],
	AF=mt.qc.AF[1],
	homozygote_count=mt.qc.homozygote_count[1]))

mt = mt.annotate_rows(qc = mt.qc.annotate(p_value_hwe = hl.case()
	.when(mt.locus.in_autosome(), mt.qc.p_value_hwe)
	.default(hl.agg.filter(mt.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt.GT).p_value)))
)

mt = mt.annotate_rows(qc = mt.qc.annotate(p_value_hwe = hl.case()
	.when(mt.locus.in_autosome(), mt.qc.het_freq_hwe)
	.default(hl.agg.filter(mt.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt.GT).het_freq_hwe)))
)

mt.rows().select('qc').flatten().export(VARIANT_QC_FILE)
