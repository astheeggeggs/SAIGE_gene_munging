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

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

PHENOFILE_CTS = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/UKBB_WES200k_filtered_cts_phenotypes.tsv.gz'
PHENOFILE_BINARY = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/UKBB_WES200k_filtered_binary_phenotypes.tsv.gz'
PHENOTYPES_TABLE = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/' + TRANCHE + '/phenotypes.ht'

pheno_table_cts = hl.import_table(PHENOFILE_CTS, key='ID', types={'ID': hl.tstr}, impute=True, missing=['', 'NA'], force_bgz=True)
pheno_table_binary = hl.import_table(PHENOFILE_BINARY, key='ID', types={'ID': hl.tstr}, impute=True, missing=['', 'NA'], force_bgz=True)

drop_cols = list(set(list(pheno_table_cts.row_value)) & set(pheno_table_binary.row_value))
pheno_table_binary = pheno_table_binary.drop(*drop_cols)

pheno_table = pheno_table_cts.annotate(**pheno_table_binary[pheno_table_cts.key])
pheno_table.write(PHENOTYPES_TABLE, overwrite=True)

ht = hl.read_table(PHENOTYPES_TABLE)
n_samples = ht.count()

print('')
print('nSamples: ', '{:,}'.format(n_samples))
pprint(ht.describe())
