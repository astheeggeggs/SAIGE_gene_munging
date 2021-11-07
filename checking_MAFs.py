import hail as hl
import argparse
from ukb_utils import hail_init

parser = argparse.ArgumentParser()
parser.add_argument("--chr", type=str, default='20')
args = parser.parse_args()

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', 'GRCh38')

ht = hl.import_table(f"/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/variant_manifest/ukb_wes_200k_variant_manifest_chr{args.chr}.tsv.gz",
                    force=True, impute=True)
ht = ht.transmute(locus = hl.locus(ht.chrom, ht.position, reference_genome = 'GRCh38'), alleles = [ht.a1, ht.a2])
ht = ht.key_by(ht.locus, ht.alleles).drop('rsid')

# Annotate with all the information in the vcfs
## Read in, split multi and remove invariant sites
ht_split_vcf = hl.import_vcf(f"/well/lindgren/ukb_wes_tmp/initial_qc/split/tmp-ukb_wes_200k_split_chr{args.chr}.vcf.gz", force_bgz=True, array_elements_required=False).rows()
ht_split_vcf = ht_split_vcf.annotate(in_variant_list = True)

## Filter genotypes and run GATK to determine excessHet
ht_excessHet_vcf = hl.import_vcf(f"/well/lindgren/ukb_wes_tmp/initial_qc/gatk_annot/tmp-ukb_wes_200k_split_filtered_gatk_chr{args.chr}.vcf.gz", force_bgz=True, array_elements_required=False).rows()

## MAFs and collection of variants present in gene-bass
ht_genebass = hl.read_matrix_table('/well/lindgren/UKBIOBANK/nbaya/resources/ukbb-exome-public/variant_results.mt').rows()
ht_genebass = ht_genebass.annotate(in_variant_list = True)

## Variant annotations following Nik's QC
ht_initial_QC = hl.read_matrix_table(f'/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr{args.chr}.mt').rows()
ht_initial_QC = ht_initial_QC.annotate(in_variant_list = True)

## Variant annotations following Duncan's subsequent QC
ht_final_QC = hl.read_matrix_table(f'/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr{args.chr}.mt').rows()
ht_final_QC = ht_final_QC.annotate(in_variant_list = True)

ht = ht.annotate(
excessHet_vcf = ht_excessHet_vcf[ht.key],
split_multi = ht_split_vcf[ht.key],
genebass = ht_genebass[ht.key],
Nik_QC = ht_initial_QC[ht.key],
Duncan_QC = ht_final_QC[ht.key]
)

ht = ht.transmute(
    excessHet_vcf_AC = ht.excessHet_vcf.info.AC[0],
    excessHet_vcf_AQ = ht.excessHet_vcf.info.AQ[0],
    excessHet_vcf_AF = ht.excessHet_vcf.info.AF[0],
    excessHet_vcf_qual = ht.excessHet_vcf.qual,
    split_multi_AQ = ht.split_multi.info.AQ[0],
    split_multi_AF = ht.split_multi.info.AF[0],
    split_multi_filters = ht.split_multi.filters,
    splti_multi_qual = ht.split_multi.qual,
    Nik_QC_AQ = ht.Nik_QC.info.AQ[0],
    Nik_QC_AC = ht.Nik_QC.info.AC[0],
    Nik_QC_AF = ht.Nik_QC.info.AF[0],
    Duncan_QC_in_variant_list = ht.Duncan_QC.in_variant_list,
    genebass_AC = ht.genebass.AC, 
    genebass_AF = ht.genebass.AF,
    genebass_gene = ht.genebass.gene,
    genebass_annotation = ht.genebass.annotation,
    genebass_in_variant_list = ht.genebass.in_variant_list,
    ref = ht.alleles[0],
    alt = ht.alleles[1]
)

ht = ht.drop('filters')
ht.export(f'/well/lindgren/UKBIOBANK/dpalmer/debugging_chr{args.chr}.tsv.bgz')
