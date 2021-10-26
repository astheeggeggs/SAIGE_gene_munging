import hail as hl
import argparse
import pandas
import os

from ukb_utils import hail_init
from ukb_utils import genotypes

# from ukb_common import *
from SAIGE_gene_munging.utils.annotations import *
from gnomad.utils.vep import process_consequences

# If running vep, ensure that set_up_vep is run,
# grab this by sourcing the hail_utils.sh file in /well/lindgren/dpalmer.

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='20')
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', 'GRCh38')

def count_variants(vep_ht_path, vep_vcf_path):
    from gnomad.utils.vep import process_consequences
    ht = hl.read_table(vep_ht_path)
    ht = process_consequences(ht)
    ht = annotate_dbnsfp(ht, vep_vcf_path)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.annotate(
        variant_id=ht.locus.contig + ':' + hl.str(ht.locus.position) + '_' + ht.alleles[0] + '/' + ht.alleles[1],
        annotation=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical,
                                           ht.vep.worst_csq_for_variant_canonical,
                                           ht.dbnsfp)
        )
    ht = ht.filter(hl.literal({'pLoF', 'LC', 'damaging_missense', 'other_missense', 'synonymous'}).contains(ht.annotation))
    print(ht.count())

TRANCHE = args.tranche
UKB_vep_output = '/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/' + TRANCHE + '/'
vep_config = "/well/lindgren/dpalmer/wes_ko_ukbb/utils/configs/vep_env.json"
groups="pLoF,damaging_missense|LC,pLoF|damaging_missense|LC,pLoF|damaging_missense,damaging_missense,other_missense,synonymous"

print('chromosome ' + args.chr)
print('tranche: ' + TRANCHE)
# input_mt_path = '/well/lindgren/UKBIOBANK/nbaya/wes_' + TRANCHE + '/ukb_wes_qc/data/filtered/ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr + '.mt'
input_mt_path = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/final_mt/' + '10_european.strict_filtered_chr' + args.chr + '.mt'
output_vep_ht_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr + '_vep.ht'
output_genemap_ht_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr + '_gene_map.ht'
output_genemap_processed_ht_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr + '_gene_map_processed.ht'
output_saige_input_path = UKB_vep_output + 'SAIGE_gene_input/' + 'ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr

mt = hl.read_matrix_table(input_mt_path)
ht = mt.rows()
ht = hl.vep(ht, vep_config)

# VEP run on sites only vcf to generate these files.
# Here are the settings.

# readonly RAW_ROOT="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/"
# readonly RAW_FILE="/ukb_wes_200k_filtered_chr${chr}.vcf.bgz"
# readonly TMP_FILE="/ukb_wes_200k_filtered_tmp_chr${chr}.vcf.gz"

# readonly OUT_ROOT="data/vep/full" #output"
# readonly OUT_FILE1="/ukb_wes_200k_full_vep_chr${chr}.vcf" 

# # Check if outfile already exists
# if [ ! -f "${OUT_ROOT}${OUT_FILE1}" ]; then

#   # extract variants in format that can be read by VEP
#   if [ ! -f "${OUT_ROOT}${TMP_FILE}" ]; then
#     module load BCFtools/1.9-foss-2018b
#     bcftools query -f '%CHROM %POS %ID %REF %ALT %AC %QUAL\n' "${RAW_ROOT}${RAW_FILE}" -o "${OUT_ROOT}${TMP_FILE}"
#   fi

# # Load modules (NOTE: load order is important)
# module purge
# set_up_vep

# ## run VEP 95 on temporary file
# vep --input_file "${OUT_ROOT}${TMP_FILE}" \
# --dir_cache /well/lindgren/flassen/software/VEP/vep95/GRCh38 \
# --assembly GRCh38 \
# --cache \
# --species "homo_sapiens" \
# --biotype \
# --canonical \
# --vcf \
# --pick \
# --total_length \
# --sift b \
# --polyphen b \
# --dir_plugins /well/lindgren/flassen/software/VEP/plugins_grch38/ \
# --plugin LoF,loftee_path:/well/lindgren/flassen/software/VEP/plugins_grch38/,human_ancestor_fa:/well/lindgren/flassen/software/VEP/loftee_human_ancestor/GRCh38/human_ancestor.fa.gz,conservation_file:/well/lindgren/flassen/software/VEP/loftee_human_ancestor/GRCh38/loftee.sql,gerp_bigwig:/well/lindgren/flassen/software/VEP/loftee_human_ancestor/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
# --plugin dbNSFP,/well/lindgren/flassen/software/VEP/dbNSFP/GRCh38/dbNSFP4.1a_grch38.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationAssessor_score,MutationAssessor_pred,MutationTaster_score,MutationTaster_pred,PROVEAN_score,PROVEAN_pred,REVEL_score,CADD_raw,CADD_phred \
# --fields "Gene,Feature,Feature_type,Consequence,IMPACT,SYMBOL,SYMBOL_SOURCE,BIOTYPE,CANONICAL,CDS_position,Protein_position,EXON,INTRON,LoF_flags,LoF_filter,LoF,SIFT,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationAssessor_score,MutationAssessor_pred,MutationTaster_score,MutationTaster_pred,PROVEAN_score,PROVEAN_pred,REVEL_score,CADD_raw,CADD_phred" \
# --output_file "${OUT_ROOT}${OUT_FILE1}.tmp" \
# --force_overwrite \
# --no_stats \
# --verbose \
# --offline

ht = process_consequences(ht)
ht.write(output_vep_ht_path, overwrite=True)

ht = hl.read_table(output_vep_ht_path)

vep_vcf_path = '/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/vep/full/ukb_wes_200k_full_vep_chr' + args.chr + '.vcf'
ht = annotate_dbnsfp(ht, vep_vcf_path)

gene_map_ht = create_gene_map_ht(ht)
gene_map_ht.write(output_genemap_ht_path, overwrite=True)

gene_map_ht = hl.read_table(output_genemap_ht_path)
gene_map_ht = post_process_gene_map_ht(gene_map_ht)
gene_map_ht.write(output_genemap_processed_ht_path, overwrite=True)

gene_map_ht = hl.read_table(output_genemap_processed_ht_path)

# # Read in the protein coding genes
# ht_coding = hl.import_table('/well/lindgren/dpalmer/genesets/data/biomart/protein_coding_genes.tsv',
#     missing=['NA', ''], impute=True, delimiter=' ')
# ht_coding = ht_coding.transmute(gene_id = ht_coding.ensembl_gene_id)
# ht_coding = ht_coding.key_by(ht_coding.gene_id).select()
# ht_coding.count()

# # Filter to protein coding genes
# gene_map_ht = gene_map_ht.key_by('gene_id')
# gene_map_ht = gene_map_ht.filter(hl.is_defined(ht_coding[gene_map_ht.key]))

# gene_ht.select(
#     group=gene_ht.gene_id + '_' + gene_ht.gene_symbol + '_' + gene_ht.annotation,
#     variant=hl.delimit(gene_ht.variants, '\t')
# ).key_by().drop('start').export(output_saige_input_path, header=False)

# Create a distinct file for each annotation
for group in groups.split(','):
    gene_ht = gene_map_ht.filter(group==gene_map_ht.annotation)
    output_saige_input_path_tmp = output_saige_input_path + '_' + group.replace('|','_') + '_saige_gene.tsv.gz'
    gene_ht.select(
        group=gene_ht.gene_id + '_' + gene_ht.gene_symbol + '_' + gene_ht.annotation,
        variant=hl.delimit(gene_ht.variants, '\t')
    ).key_by().drop('start').export(output_saige_input_path_tmp, header=False)

output_vep_ht_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr + '_vep.ht'
output_summary_path = UKB_vep_output + 'ukb_wes_' + TRANCHE + '_filtered_chr' + args.chr + '_summary.ht'
count_variants(output_vep_ht_path, vep_vcf_path)

ht = hl.read_table(output_vep_ht_path)
ht = process_consequences(ht)
ht = annotate_dbnsfp(ht, vep_vcf_path)
ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
ht = ht.group_by(
    gene=ht.vep.worst_csq_by_gene_canonical.gene_symbol,
    consequence=annotation_case_builder(ht.vep.worst_csq_by_gene_canonical,
                                        ht.vep.worst_csq_for_variant_canonical,
                                        ht.dbnsfp)
    ).partition_hint(100).aggregate(n_variants=hl.agg.count())
ht.write(output_summary_path, overwrite=True)
