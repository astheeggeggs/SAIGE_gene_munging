library(data.table)

# First, extract the sample IDs and sequencing batch information
ukb_fam_lindgren <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
ukb_sample_qc <- "/well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt"
ukb_k_means_EUR <- "/well/lindgren/UKBIOBANK/laura/k_means_clustering_pcs/ukbb_genetically_european_k4_4PCs_self_rep_Nov2020.txt"
ukb_pheno <- "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv"
ukb_batch <- "/well/lindgren/UKBIOBANK/ferreira/UKBB_phenotypes/UKBB_WES200k_biomarkers_phenotypes.txt"
ukb_RF_EUR <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink_EUR.tsv"
ukb_RF_NFE <- "/well/lindgren/dpalmer/ukb_get_EUR/data/final_EUR_list.tsv"

aux_pheno_out <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k/QC_phenotypes.tsv"

dt <- cbind(fread(ukb_fam_lindgren), fread(ukb_sample_qc))
drop_cols <- c("V2", "V3", "V4", "V5", "V6")
dt[, (drop_cols):=NULL]
names(dt)[1] <- "eid"
setkey(dt, 'eid')

# Grab my EUR and Non-Finnish European subset of the genotype data.
dt_EUR_new <- fread(ukb_RF_EUR, col.names = c('eid', 'V2'), key='eid')
dt_EUR_new[, V2 := NULL]
dt_EUR_new[, genetic.eur.oct2021 := TRUE]

dt_EUR_no_FIN_new <- fread(ukb_RF_NFE, col.names = 'eid', key='eid')
dt_EUR_no_FIN_new[, genetic.eur.no.fin.oct2021 := TRUE]

dt <- merge(dt, merge(dt_EUR_new, dt_EUR_no_FIN_new, all=TRUE), all=TRUE)
dt[, genetic.eur.oct2021 := ifelse(is.na(genetic.eur.oct2021), FALSE, genetic.eur.oct2021)]
dt[, genetic.eur.no.fin.oct2021 := ifelse(is.na(genetic.eur.no.fin.oct2021), FALSE, genetic.eur.no.fin.oct2021)]

# Then, go to the raw phenotype file and find the UKBB centre and self reported ancestry.
dt_pheno <- fread(ukb_pheno, sep=",", select = c("eid", "54-0.0", "21000-0.0"), key='eid')
dt <- merge(dt, dt_pheno, all=TRUE)

# Grab Laura's European definition
dt_PC <- fread(ukb_k_means_EUR, key='eid')
dt_PC[, genetically_european := ifelse(genetically_european == 1, TRUE, FALSE)]
dt <- merge(dt, dt_PC, all=TRUE)

# Grab the sequencing batch from one of the cts phenotype files
dt_batch <- fread(ukb_batch, select = c("ID", "sequencing.batch"))
setnames(dt_batch, "ID", "eid")
setkey(dt_batch, "eid")

dt <- merge(dt, dt_batch, all.x=TRUE)
setnames(dt, c("eid", "54-0.0", "21000-0.0"), c("ID", "ukbb.centre", "self.report.ethnicity"))

# Write the result to disk
fwrite(dt, aux_pheno_out, sep='\t')
system(paste("bgzip", aux_pheno_out))
