library(data.table)
library(dplyr)

eur_fam <- "/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/ukb_imp_eur_chr1_22_sparse_markers.fam"
ukb_fam <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
out <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_imp_eur_chr1_22_sparse_markers.fam"
fwrite(merge(fread(ukb_fam, key="V2"), fread(eur_fam, key="V2", select="V2")) %>% mutate(V6=-9), file=out, sep='\t', col.names=FALSE)
