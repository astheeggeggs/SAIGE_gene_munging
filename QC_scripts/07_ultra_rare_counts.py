import hail as hl
from hail.plot import show
from pprint import pprint
hl.plot.output_notebook()

import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, required=True)
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche
CHR = str(args.chr)

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

MT_HARDCALLS = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_hardcalls_chr' + CHR + '.mt'

ANNOTATION_TABLE = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '_vep_qc.ht'
URV_FILE = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/07_URVs_chr' + CHR + '.tsv.bgz'
PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'

INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'

INITIAL_VARIANT_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/variants/02_prefilter_chr' + CHR +'.keep.variant.ht'
SEXCHECK_LIST = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/04_sexcheck.remove.sample_list'

ht_initial_variants = hl.read_table(INITIAL_VARIANT_LIST)
ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')
ht_sexcheck_samples = hl.import_table(SEXCHECK_LIST, no_header=True, key='f0')

# Annotate variants with counts from non-psychiatric version of Gnomad.
# Fill in variants not in Gnomad variant list.
# Annotate variants with LoF/damaging missense annotation.

mt = hl.read_matrix_table(MT_HARDCALLS)
mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
mt = mt.filter_cols(~hl.is_defined(ht_sexcheck_samples[mt.col_key]))
mt = mt.filter_rows(hl.is_defined(ht_initial_variants[mt.row_key]))

# Drop some fields that are not needed.
mt = mt.drop('rsid', 'info', 'filters')

sample_annotations = hl.read_table(PHENOTYPES_TABLE)
consequence_annotations = hl.read_table(ANNOTATION_TABLE)
output_vep_ht_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '_vep.ht'

mt = mt.annotate_cols(phenotype = sample_annotations[mt.s])
mt = mt.annotate_rows(consequence = consequence_annotations[mt.row_key])
mt = mt.annotate_rows(is_singleton = hl.agg.sum(mt.GT.n_alt_alleles()) == 1)
mt = mt.filter_rows((mt.is_singleton) & (~mt.consequence.inGnomAD))

mt = mt.annotate_cols(n_coding_URV_SNP = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_snp(mt.alleles[0], mt.alleles[1]) & (mt.consequence.annotation != "non-coding")),
                      n_coding_URV_indel = hl.agg.count_where(mt.GT.is_non_ref() & hl.is_indel(mt.alleles[0], mt.alleles[1]) & (mt.consequence.annotation != "non-coding")),
                      n_URV_pLoF = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence.annotation == "pLoF")),
                      n_URV_missense = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence.annotation == "missense")),
                      n_URV_synonymous = hl.agg.count_where(mt.GT.is_non_ref() & (mt.consequence.annotation == "synonymous")))

mt.cols().flatten().export(URV_FILE)
