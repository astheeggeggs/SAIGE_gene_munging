library(data.table)
library(dplyr)

# source activate /well/lindgren/users/mmq446/conda/skylake/envs/RSAIGE
# module load BCFtools

library(SAIGE, lib.loc='/well/lindgren/flassen/software/tmp/') 

TRANCHE <- '200k'

# Tables to write
output_folder <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/", TRANCHE)
filtered_output <- paste0(output_folder, '/UKBB_WES', TRANCHE , "_filtered_dec2021_phenotypes.tsv")
filtered_output_imp <- paste0(output_folder, '/UKBB_WES', TRANCHE , "_filtered_dec2021_imputed_phenotypes.tsv")  

# Tables to merge
# Quantitative phenotypes
file <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/curated_phenotypes_binary.tsv"
file <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/curated_phenotypes_cts.tsv"

# Only thing extra to include is the sequencing batch.
dt_sb <- fread("/well/lindgren/UKBIOBANK/ferreira/UKBB_phenotypes/UKBB_WES200k_physicalmeasures_phenotypes.txt", select = c("ID", "sequencing.batch"))
setnames(dt_sb, c("ID"), c("eid"))
dt_sb[, eid:=as.character(eid)]

dt <- fread(file)
dt[, eid:=as.character(eid)]
setkey(dt, "eid")
setkey(dt_sb, "eid")

dt <- merge(dt, dt_sb, all.x=TRUE)
dt[, sex := as.factor(sex)]
dt[, ukbb.centre:=as.factor(ukbb.centre)]
dt[, sequencing.batch:=as.factor(sequencing.batch - 1)]

# Read and merge in our new definitions of Europeans and non-Finnish Europeans.
dt_EUR_new <- fread("/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink_EUR.tsv") %>% transmute(eid = as.character(V1)) %>% mutate(genetic.eur.oct2021 = TRUE)
dt_EUR_no_FIN_new <- fread("/well/lindgren/dpalmer/ukb_get_EUR/data/final_EUR_list.tsv") %>% transmute(eid = as.character(V1)) %>% mutate(genetic.eur.no.fin.oct2021 = TRUE)
# Read in the full fam file which was used for ancestry assignment:
ukb_fam <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam"
dt_fam <- fread(ukb_fam) %>% transmute(eid = as.character(V1))
dt_fam <- data.table(dt_fam)
setkey(dt_fam, "eid")
setkey(dt_EUR_new, "eid")
setkey(dt_EUR_no_FIN_new, "eid")
dt_EUR_new <- merge(merge(dt_EUR_new, dt_EUR_no_FIN_new, all.x=TRUE), dt_fam, all.y=TRUE)
dt_EUR_new[, genetic.eur.oct2021 := ifelse(is.na(genetic.eur.oct2021), FALSE, genetic.eur.oct2021)]
dt_EUR_new[, genetic.eur.no.fin.oct2021 := ifelse(is.na(genetic.eur.no.fin.oct2021), FALSE, genetic.eur.no.fin.oct2021)]

dt <- merge(dt, dt_EUR_new, all.x=TRUE)

# Parentheses do not play nice with SAIGE.
names(dt) <- gsub("[\\(\\)]", "", names(dt))
names(dt) <- gsub("\\-", "_", names(dt))

# We also need to ensure that all samples with phenotype information are present in the vcf.
# We can do this by making use of the savvy library in SAIGE.
vcfFile <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr21.vcf.bgz")
vcfFileIndex <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr21.vcf.bgz.csi")
vcfField <- "GT"
isVariant <- setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
sampleListinvcf <- data.table(eid = getSampleIDlist_vcfMatrix())
setkey(sampleListinvcf, "eid")
dt <- merge(sampleListinvcf, dt)

# Remove all the withdrawn samples
withdrawn <- fread("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210809.csv") %>% transmute(eid = V1)

dt <- dt[which(!dt$eid %in% withdrawn$eid), ]

fwrite(dt, file=filtered_output, sep='\t')
system(paste("bgzip", filtered_output))

# Create a final two files that in addition is filtered to the collection of samples that have IMPUTED genotype information available.
dt <- fread(cmd = paste("zcat", paste0(filtered_output, ".gz")), key="ID")

imputed_EUR_fam <- fread(
    paste0("/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211026_long_ukb_wes_",
        TRANCHE, "_sparse_autosomes.fam")
    ) %>% transmute(eid = as.character(V2))
imputed_EUR_fam <- data.table(imputed_EUR_fam)
setkey(imputed_EUR_fam, "eid")
dt <- merge(dt, imputed_EUR_fam)

fwrite(dt, file=filtered_output_imp, sep='\t')
system(paste("bgzip", filtered_output_imp))
