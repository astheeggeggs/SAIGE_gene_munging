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
MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'
PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'
IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/05_imputesex.ht'
UKB_vep_output = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/' + TRANCHE + '/'
ANNOTATION_TABLE = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '_vep_qc.ht'
FINAL_SAMPLE_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list'
FINAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list'

# Outputs:
QC_MT = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr' + CHR + '.mt'
QC_HARDCALLS_MT = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr' + CHR + '.hardcalls.mt'
QC_HARDCALLS_VCF = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr' + CHR + '.hardcalls.vcf.bgz'

ht_final_samples = hl.import_table(FINAL_SAMPLE_LIST, no_header=True, key='f0')
ht_final_variants = hl.import_table(FINAL_VARIANT_LIST, types={'locus':hl.tlocus(reference_genome='GRCh38'), 'alleles':hl.tarray(hl.tstr)})
ht_final_variants = ht_final_variants.key_by(ht_final_variants.locus, ht_final_variants.alleles)

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
impute_sex_annotations = hl.read_table(IMPUTESEX_TABLE)
annotation_annotations = hl.read_table(ANNOTATION_TABLE)

mt = hl.read_matrix_table(MT)
mt = mt.drop('qual', 'info', 'filters')

mt = mt.filter_cols(hl.is_defined(ht_final_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_final_variants[mt.row_key]))

mt = mt.annotate_cols(phenotype = sample_annotations[mt.col_key])
mt = mt.annotate_cols(imputesex = impute_sex_annotations[mt.col_key])
mt = mt.annotate_rows(annotation = annotation_annotations[mt.row_key])

mt = hl.variant_qc(mt, name = 'qc')

mt = mt.annotate_rows(qc = mt.qc.annotate(p_value_hwe = hl.case()
	.when(mt.locus.in_autosome(), mt.qc.het_freq_hwe)
	.default(hl.agg.filter(mt.imputesex.impute_sex.is_female,
		hl.agg.hardy_weinberg_test(mt.GT).het_freq_hwe)))
)

mt = hl.sample_qc(mt)
mt.write(QC_MT, overwrite=True)
mt.select_entries(mt.GT).repartition(512).write(QC_HARDCALLS_MT, overwrite=True)

mt = hl.read_matrix_table(QC_HARDCALLS_MT)
mt = mt.select_rows()
mt = mt.select_cols()

# Drop all annotations except keys
export_vcf(mt, output=QC_HARDCALLS_VCF)
