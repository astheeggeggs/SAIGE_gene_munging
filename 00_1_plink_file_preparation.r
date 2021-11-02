library(data.table)
library(dplyr)

# Inputs
ukb_fam <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"

# Outputs
out_all <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam"

# Separately, generate the fam file for all people, so that that it can be read by plink.
fwrite(fread(ukb_fam) %>% mutate(V6=-9), file=out_all, sep=' ', col.names=FALSE)

# The above is generated before running the QC pipeline, the below is based on the final sample set following the QC pipeline.

eur_fam <- "/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes.fam"
out_sparse_GRM_samples <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_imp_eur_chr1_22_sparse_markers.tsv"
out_sparse_fam <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_imp_eur_chr1_22_sparse_markers.fam"

fwrite(fread(eur_fam) %>% mutate(V1=V2) %>% select(V1, V2), file = out_sparse_GRM_samples, sep=' ', col.names=FALSE)
fwrite(fread(eur_fam) %>% mutate(V1=V2, V6=-9), file = out_sparse_fam, sep='\t', col.names=FALSE)