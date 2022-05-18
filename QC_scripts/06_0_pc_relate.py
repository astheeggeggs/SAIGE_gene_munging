import hail as hl
import argparse

from ukb_utils import hail_init
from ukb_utils import genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--tranche", type=str, default='200k')
args = parser.parse_args()

TRANCHE = args.tranche
PC_RELATE_KINSHIP_THRESHOLD = 0.08838835
SAMPLE_LIST_PC_RELATED = '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/06_pc_relate.related.sample_list'
PC_RELATE_OUTPUT = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/06_pc_relate.tsv.bgz'
INITIAL_SAMPLES = '/well/lindgren/UKBIOBANK/dpalmer/wes_' + TRANCHE + '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list'

hail_init.hail_bmrc_init_local('logs/hail/hail_export.log', 'GRCh38')

# # Inputs
# # UKB genotype calls
# ukb_bed="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_pruned_autosomes.bed"
# ukb_bim="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_pruned_autosomes.bim"
# ukb_fam="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_pruned_autosomes.fam"

# ht_initial_samples = hl.import_table(INITIAL_SAMPLES, no_header=True, key='f0')

# # Read in the plink file
# mt = hl.import_plink(bed=ukb_bed, bim=ukb_bim, fam=ukb_fam, reference_genome='GRCh37')
# mt = mt.key_rows_by(mt.rsid)
# mt = mt.filter_cols(hl.is_defined(ht_initial_samples[mt.col_key]))
# mt = mt.key_rows_by(mt.locus, mt.alleles).checkpoint("/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_pruned_autosomes.mt")

# # Read in the collection of PCs from the UK Biobank phenotype file.
# # I've filtered to the columns and extracted the first 10 PCs in this file:
# # /well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_pcs.tsv
# # dt <- fread("/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv")
# # dt <- dt %>% select(eid, `22009-0.1`, `22009-0.2`, `22009-0.3`, `22009-0.4`, `22009-0.5`,
# # 	`22009-0.6`, `22009-0.7`, `22009-0.8`, `22009-0.9`, `22009-0.10`)
# # names(dt) <- c("eid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
# # fwrite(dt, "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_pcs.tsv", sep="\t")

# # ht_pcs = hl.import_table("/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_pcs.tsv",
# # 	impute=True, missing=['NA', ''])
# # ht_pcs = ht_pcs.annotate(eid=hl.str(ht_pcs.eid)).key_by('eid')
# # ht_pcs = ht_pcs.annotate(PCs = [ht_pcs.PC1, ht_pcs.PC2, ht_pcs.PC3, ht_pcs.PC4]).select("PCs")

mt = mt.read_matrix_table("/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_pruned_autosomes.mt")
mt = mt.annotate_rows(a=(hl.rand_unif(0,1) < 0.01))
mt = mt.annotate_cols(c=(hl.rand_unif(0,1) < 0.01))
mt_tmp = mt.filter_cols(mt.c)
mt_tmp = mt_tmp.filter_rows(mt_tmp.b)

_, scores_table, _ = hl.hwe_normalized_pca(mt_tmp.GT, k=2, compute_loadings=False)

n = hl.count(scores_table)
print(n)
# # Not tested, if this is slow or doesn't work well - consider passing the projected PCs from step 05_estimate_superpopulation.r
# pc_rel = hl.pc_relate(mt.GT, 0.01, k=4, statistics='kin')
# pairs = pc_rel.filter(pc_rel['kin'] > PC_RELATE_KINSHIP_THRESHOLD)

# related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
# result = mt.cols().select()
# result.filter(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=True).export(SAMPLE_LIST_PC_RELATED)

# pc_relate_table = pc_relate_table.flatten()
# pc_relate_table.export(PC_RELATE_OUTPUT)
