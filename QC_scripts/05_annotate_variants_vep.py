import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

from ukb_common import *
from gnomad.utils.vep import process_consequences

# If running vep, ensure that set_up_vep is run,
# grab this by sourcing the hail_utils.sh file in /well/lindgren/dpalmer.

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='20')
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche
CHR = str(args.chr)

print('chromosome ' + args.chr)
print('tranche: ' + TRANCHE)

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', 'GRCh38')

# Inputs
vep_config = "/well/lindgren/dpalmer/wes_ko_ukbb/utils/configs/vep_env.json"
GNOMAD_SITES_38_HT = '/well/lindgren/flassen/ressources/gnomad/gnomad_v2_liftover/exomes/gnomad.exomes.r2.1.1.sites.' + CHR + '.liftover_grch38.vcf.bgz'
input_mt_path = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '.mt'

# Outputs
UKB_vep_output = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/' + TRANCHE + '/'
output_vep_ht_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + CHR + '_vep_qc.ht'

mt = hl.read_matrix_table(input_mt_path)
ht = mt.rows()
ht = hl.vep(ht, vep_config)
ht = process_consequences(ht)

# Annotate with gnomAD information for the relevant chromosome
gnomad_variants_ht = hl.import_vcf(GNOMAD_SITES_38_HT, reference_genome ='GRCh38', force_bgz=True, array_elements_required=False).rows()
ht = ht.annotate(inGnomAD = hl.is_defined(gnomad_variants_ht[ht.key]))
ht = ht.annotate(annotation=annotation_case_builder(ht.vep.worst_csq_for_variant_canonical))

ht.write(output_vep_ht_path, overwrite=True)