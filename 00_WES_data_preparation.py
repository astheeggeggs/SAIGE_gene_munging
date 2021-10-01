import hail as hl
import argparse
import pandas
import os

from ukb_utils import hail_init
from ukb_utils import genotypes

from ukb_common import *

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, required=True)
args = parser.parse_args()

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default=22)
args = parser.parse_args()

hail_init.hail_bmrc_init('logs/hail/hail_export.log', 'GRCh38')

TRANCHE = '200k'

print('chromosome ' + args.chr)
input_filtered_mt_path = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr + '.mt'
input_raw_vcf = '/well/ukbb-wes/pvcf/oqfe/ukbb-wes-oqfe-pvcf-chr' + args.chr + '.vcf.gz'
output_mt_path = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_and_corrected_chr' + args.chr + '.mt'
output_vcf_path = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_and_corrected_chr' + args.chr + '.vcf.bgz'

mt_raw = hl.import_vcf(input_raw_vcf, reference_genome ='GRCh38', force_bgz=True, array_elements_required=False)
bridge = hl.import_table('/well/lindgren/UKBIOBANK/DATA/WES/bridge_11867_12788.csv', delimiter=',').key_by('eid_12788')
mt_raw = mt_raw.annotate_cols(**bridge[mt_raw.col_key.s]).key_cols_by().drop('s')
mt_raw = mt_raw.transmute_cols(s = mt_raw.eid_11867).key_cols_by('s')
mt_raw = mt_raw.filter_rows(mt_raw.alleles.length() <= 6)
mt_raw = hl.split_multi_hts(mt_raw)

mt_filtered = hl.read_matrix_table(input_filtered_mt_path)
mt_raw_filtered = mt_raw.filter_cols(hl.is_defined(mt_filtered.cols()[mt_raw.col_key]))
mt_raw_filtered = mt_raw_filtered.filter_rows(hl.is_defined(mt_filtered.rows()[mt_raw_filtered.row_key]))

# Annotate entries with the dosage, which we can pass to SAIGE instead of hard calls.
# mt_raw_filtered = mt_raw_filtered.annotate_entries(DS = hl.pl_dosage(mt_raw_filtered.PL))
mt_raw_filtered = hl.variant_qc(mt_raw_filtered)

# Export the resultant mt and vcfs
mt_raw_filtered.write(output_mt_path)

# # Chat with Nik about the best way to run this last piece.
# mt_raw_filtered.select_entries(mt_raw_filtered.GT, mt_raw_filtered.DS).export_vcf(output=output_vcf_path, parallel="header_per_shard")

# mt_raw_filtered = hl.read_matrix_table(output_mt_path)
# print(mt_raw_filtered.count())
