import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', 'GRCh38')

# # Exome
# MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + args.chr + '.mt'

# IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex.ht'
# IMPUTESEX_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex.tsv'
# Y_NCALLED = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_ycalled.tsv'

# # Need to perform LD pruning to determine sex estimation.
# INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'

# PRUNED_CHRX_VARIANTS = 'gs://dalio_bipolar_w1_w2_hail_02/data/variants/04_chrX.prune.in'

# PHENOTYPES_TABLE = 'gs://dalio_bipolar_w1_w2_hail_02/data/samples/phenotypes.ht'

# ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
# ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
# sample_annotations = hl.read_table(PHENOTYPES_TABLE)

# ht_pruned_chrx_variants = ht_pruned_chrx_variants.annotate(**hl.parse_variant(ht_pruned_chrx_variants.f0, reference_genome='GRCh38'))
# ht_pruned_chrx_variants = ht_pruned_chrx_variants.key_by(ht_pruned_chrx_variants.locus, ht_pruned_chrx_variants.alleles)

# mt = hl.read_matrix_table(MT_HARDCALLS)
# mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
# mt = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))

# n = mt.count()

# print('n samples:')
# print(n[1])
# print('n variants:')
# print(n[0])

# imputed_sex = hl.impute_sex(mt.GT, female_threshold=0.6, male_threshold=0.6)
# mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
# mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

# mt.cols().select('impute_sex', 'phenotype').flatten().export(IMPUTESEX_FILE)
# # Want to change this to reflect the dataset that I have.
# mt.cols().write(IMPUTESEX_TABLE, overwrite=True)

# # Determine non-missing allele count on the y.
# mt = hl.read_matrix_table(MT_HARDCALLS)
# mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
# mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
# mt = hl.sample_qc(mt, name='qc')

# mt_cols = mt.cols()
# mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)

# Genotypes

# Inputs
# UKB genotype calls
ukb_bed_X="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chrX_v2.bed"
ukb_bim_X="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chrX_v2.bim"
ukb_bed_Y="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chrY_v2.bed"
ukb_bim_Y="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chrY_v2.bim"
ukb_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'
PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'
PRUNED_CHRX_VARIANTS = '/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chrX_pruned.prune.in'

# Outputs
IMPUTESEX_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex.ht'
IMPUTESEX_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_imputesex.tsv.bgz'
Y_NCALLED = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_ycalled.tsv.bgz'

ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_pruned_chrx_variants = hl.import_table(PRUNED_CHRX_VARIANTS, no_header=True)
ht_pruned_chrx_variants = ht_pruned_chrx_variants.transmute(rsid=ht_pruned_chrx_variants.f0).key_by('rsid')
sample_annotations = hl.read_table(PHENOTYPES_TABLE)

# Read in the plink file
mt = hl.import_plink(bed=ukb_bed_X, bim=ukb_bim_X, fam=ukb_fam, reference_genome='GRCh37')
mt = mt.key_rows_by(mt.rsid)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_pruned_chrx_variants[mt.row_key]))
mt = mt.key_rows_by(mt.locus, mt.alleles)
n = mt.count()

print('n samples:')
print(n[1])
print('n variants:')
print(n[0])

imputed_sex = hl.impute_sex(mt.GT, female_threshold=0.6, male_threshold=0.6)
mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_cols(impute_sex = imputed_sex[mt.s])

mt.cols().select('impute_sex', 'phenotype').flatten().export(IMPUTESEX_FILE)
mt.cols().write(IMPUTESEX_TABLE, overwrite=True)

# Now, look on the Y chromosome, and determine non-missing allele count on the y.
mt = hl.import_plink(bed=ukb_bed_Y, bim=ukb_bim_Y, fam=ukb_fam, reference_genome='GRCh37')
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_rows(mt.locus.in_y_nonpar() | mt.locus.in_y_par())
mt = hl.sample_qc(mt, name='qc')

mt_cols = mt.cols()
mt_cols.select(n_called=mt_cols.qc.n_called).export(Y_NCALLED)
