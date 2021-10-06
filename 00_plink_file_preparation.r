library(data.table)
library(dplyr)

eur_fam <- "/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/ukb_imp_eur_chr1_22_sparse_markers.fam"
ukb_fam <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
out <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_imp_eur_chr1_22_sparse_markers.fam"
out_all <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam"
out_sparse_GRM_samples <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_imp_eur_chr1_22_sparse_markers.tsv"

# Separately, generate the fam file for all people, so that that it can be read by plink.
fwrite(fread(ukb_fam) %>% mutate(V6=-9), file=out_all, sep=' ', col.names=FALSE)
fwrite(fread(eur_fam) %>% mutate(V1=V2) %>% select(V1, V2), file = out_sparse_GRM_samples, sep=' ')

# DEV: Create a plink command instead.
system("/well/lindgren/dpalmer/plink ")


