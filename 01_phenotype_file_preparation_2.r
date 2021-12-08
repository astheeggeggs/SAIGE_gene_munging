library(data.table)
library(dplyr)

# source activate /well/lindgren/users/mmq446/conda/skylake/envs/RSAIGE
# module load BCFtools

library(SAIGE, lib.loc='/well/lindgren/flassen/software/tmp/') 

TRANCHE <- '200k'

# Tables to write
output_folder <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/", TRANCHE)
cts_output <- paste0(output_folder, '/UKBB_WES', TRANCHE , "_cts_phenotypes.tsv")  
binary_output <- paste0(output_folder, '/UKBB_WES', TRANCHE, "_binary_phenotypes.tsv")

cts_filtered_output <- paste0(output_folder, '/UKBB_WES', TRANCHE , "_filtered_cts_phenotypes.tsv")  
binary_filtered_output <- paste0(output_folder, '/UKBB_WES', TRANCHE, "_filtered_binary_phenotypes.tsv")

# Tables to merge
# Quantitative phenotypes
folder <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/"
files <- c("UKBB_WES_plos_genetics_phenotypes.txt", "UKBB_WES_biomarkers_phenotypes.txt")
files <- paste0(folder, files)

# Only thing extra to include is the sequencing batch.

# Cols to drop in dt_tmp
covariate_cols <- c("age", "ukbb.centre", paste0("PC", seq(1,40)), "sex", "sequencing.batch") 

dt <- fread(files[1])
dt[, eid:=as.character(eid)]
setkey(dt, "eid")

if (length(files) > 1)
{
    for (file in files[2:length(files)]) {
    	dt_tmp <- fread(file)
        dt_tmp[, eid:=as.character(eid)]
    	setkey(dt_tmp, "eid")
    	dt_tmp[, (covariate_cols):=NULL]
    	dt <- merge(dt, dt_tmp, all=TRUE)
    }
}

# Filter to the full set of non inverse rank normalised phenotypes, which we will pass to SAIGE

cts_phenotypes <- c(
    # Log of the biomarkers
    "Alanine_aminotransferase",
    "Albumin",
    "Alkaline_Phosphatase",
    "Apolipoprotein_A",
    "Apolipoprotein_B",
    "Aspartate_aminotransferase",
    "C_reactive_Protein",
    "Calcium",
    "Cholesterol",
    "Creatinine_Serum",
    "Creatinine_Urine",
    "Cystatin_C_Serum",
    "Direct_Bilirubin",
    "Gamma_glutamyltransferase",
    "Glucose",
    "HbA1c",
    "HDL_Cholesterol",
    "IGF_1",
    "Lipoprotein_A",
    "Microalbumin_Urine",
    "Oestradiol",
    "Phosphate",
    "Potassium_Urine",
    "Rheumatoid_factor",
    "SHBG",
    "Sodium_Urine",
    "Testosterone",
    "Total_Bilirubin",
    "Total_Protein",
    "Triglyceride",
    "Urate",
    "Urea",
    "Vitamin_D",
    # Residualised log biomarkers. Both sexes, followed by sex specific residualisation of log biomarkers.
    "Alanine_aminotransferase_residual",
    "Alanine_aminotransferase_F_residual",
    "Alanine_aminotransferase_M_residual",
    "Albumin_residual",
    "Albumin_F_residual",
    "Albumin_M_residual",
    "Alkaline_Phosphatase_residual",
    "Alkaline_Phosphatase_F_residual",
    "Alkaline_Phosphatase_M_residual",
    "Apolipoprotein_A_residual",
    "Apolipoprotein_A_F_residual",
    "Apolipoprotein_A_M_residual",
    "Apolipoprotein_B_residual",
    "Apolipoprotein_B_F_residual",
    "Apolipoprotein_B_M_residual",
    "Aspartate_aminotransferase_residual",
    "Aspartate_aminotransferase_F_residual",
    "Aspartate_aminotransferase_M_residual",
    "C_reactive_Protein_residual",
    "C_reactive_Protein_F_residual",
    "C_reactive_Protein_M_residual",
    "Calcium_residual",
    "Calcium_F_residual",
    "Calcium_M_residual",
    "Cholesterol_residual",
    "Cholesterol_F_residual",
    "Cholesterol_M_residual",
    "Creatinine_Serum_residual",
    "Creatinine_Serum_F_residual",
    "Creatinine_Serum_M_residual",
    "Creatinine_Urine_residual",
    "Creatinine_Urine_F_residual",
    "Creatinine_Urine_M_residual",
    "Cystatin_C_Serum_residual",
    "Cystatin_C_Serum_F_residual",
    "Cystatin_C_Serum_M_residual",
    "Direct_Bilirubin_residual",
    "Direct_Bilirubin_F_residual",
    "Direct_Bilirubin_M_residual",
    "Gamma_glutamyltransferase_residual",
    "Gamma_glutamyltransferase_F_residual",
    "Gamma_glutamyltransferase_M_residual",
    "Glucose_residual",
    "Glucose_F_residual",
    "Glucose_M_residual",
    "HbA1c_residual",
    "HbA1c_F_residual",
    "HbA1c_M_residual",
    "HDL_Cholesterol_residual",
    "HDL_Cholesterol_F_residual",
    "HDL_Cholesterol_M_residual",
    "IGF_1_residual",
    "IGF_1_F_residual",
    "IGF_1_M_residual",
    "Lipoprotein_A_residual",
    "Lipoprotein_A_F_residual",
    "Lipoprotein_A_M_residual",
    "Microalbumin_Urine_residual",
    "Microalbumin_Urine_F_residual",
    "Microalbumin_Urine_M_residual",
    "Oestradiol_residual",
    "Oestradiol_F_residual",
    "Oestradiol_M_residual",
    "Phosphate_residual",
    "Phosphate_F_residual",
    "Phosphate_M_residual",
    "Potassium_Urine_residual",
    "Potassium_Urine_F_residual",
    "Potassium_Urine_M_residual",
    "Rheumatoid_factor_residual",
    "Rheumatoid_factor_F_residual",
    "Rheumatoid_factor_M_residual",
    "SHBG_residual",
    "SHBG_F_residual",
    "SHBG_M_residual",
    "Sodium_Urine_residual",
    "Sodium_Urine_F_residual",
    "Sodium_Urine_M_residual",
    "Testosterone_residual",
    "Testosterone_F_residual",
    "Testosterone_M_residual",
    "Total_Bilirubin_residual",
    "Total_Bilirubin_F_residual",
    "Total_Bilirubin_M_residual",
    "Total_Protein_residual",
    "Total_Protein_F_residual",
    "Total_Protein_M_residual",
    "Triglyceride_residual",
    "Triglyceride_F_residual",
    "Triglyceride_M_residual",
    "Urate_residual",
    "Urate_F_residual",
    "Urate_M_residual",
    "Urea_residual",
    "Urea_F_residual",
    "Urea_M_residual",
    "Vitamin_D_residual",
    "Vitamin_D_F_residual",
    "Vitamin_D_M_residual"
    )

binary_phenotypes <- c(
    "BC_combined",
    "CAD_combined",
    "COPD_combined",
    "CLD_combined",
    "CC_combined",
    "DEM_combined",
    "INF_combined"       
    "LC_combined",
    "NAFLD_combined",
    "RF_combined",
    "RF_acute_combined",
    "RF_chronic_combined",
    "STR_combined"       
    "STR_hem_combined",
    "STR_isc_combined"   
    )

# Phenotypes we don't have from the biomarkers.
cts_phenotypes <- c(
    "Visceral_adipose_tissue_volume_(VAT)",
    "Total_adipose_tissue_volume",
    "Abdominal_fat_ratio",
    "Liver_proton_density_fat_fraction_(AMRA)",
    "Body_mass_index_(BMI)",
    "Hip_circumference",
    "Standing_height",
    "Waist_circumference",
    "Body_fat_percentage"
    )

binary_phenotypes <- c(
    "hypothalamic_amenorrhea",
    "POI",
    "Alzheimers_disease",
    "depression",
    "autism",
    "ADHD",
    "ischaemic_heart_disease", # What's the difference between this and CAD?
    "Crohns_disease",
    "IBD",
    "Cirrhosis",
    "NASH",
    "psoriasis",
    "hyperandrogenism",
    "hematuria",
    "proteinuria",
    "oligomenorrhea",
    "habitual_aborter",
    "ectopic_pregnancy",
    "Preeclampsia",
    "GDM",
    "intrahepatic_cholestasis_in_pregnancy",
    "polycystic_kidney_disease",
    "T2D", # Grab these from Jenny's definition
    "T1D", # Grab from Jenny's definition
    "GDM2",
    "kallmann_syndrome", # This one is carefully defined - grab from here.
    "E230", #? Hypopituitarism?
    "PCOS4" # This one is carefully definied - grab from here.
    )

retain_cols <- c("ID", drop_cols, cts_phenotypes, binary_phenotypes)
dt <- dt[, ..retain_cols]
dt[, sex := as.factor(sex)]
dt[, ukbb.centre:=as.factor(ukbb.centre)]
dt[, genotyping.array:=as.factor(genotyping.array)]
dt[, sequencing.batch:=as.factor(sequencing.batch)]

# Read and merge in our new definitions of Europeans and non-Finnish Europeans.
dt_EUR_new <- fread("/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink_EUR.tsv") %>% transmute(ID = as.character(V1)) %>% mutate(genetic.eur.oct2021 = TRUE)
dt_EUR_no_FIN_new <- fread("/well/lindgren/dpalmer/ukb_get_EUR/data/final_EUR_list.tsv") %>% transmute(ID = as.character(V1)) %>% mutate(genetic.eur.no.fin.oct2021 = TRUE)
# Read in the full fam file which was used for ancestry assignment:
ukb_fam <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam"
dt_fam <- fread(ukb_fam) %>% transmute(ID = as.character(V1))
dt_fam <- data.table(dt_fam)
setkey(dt_fam, "ID")
setkey(dt_EUR_new, "ID")
setkey(dt_EUR_no_FIN_new, "ID")
dt_EUR_new <- merge(merge(dt_EUR_new, dt_EUR_no_FIN_new, all.x=TRUE), dt_fam, all.y=TRUE)
dt_EUR_new[, genetic.eur.oct2021 := ifelse(is.na(genetic.eur.oct2021), FALSE, genetic.eur.oct2021)]
dt_EUR_new[, genetic.eur.no.fin.oct2021 := ifelse(is.na(genetic.eur.no.fin.oct2021), FALSE, genetic.eur.no.fin.oct2021)]

dt <- merge(dt, dt_EUR_new, all.x=TRUE)
dt[, genetic.eur := ifelse(genetic.eur == 1, TRUE, FALSE)]

filter_to <- c("ID", "age", paste0("PC", seq(1,40)), "ukbb.centre", "sex",
    "genotyping.array", "sequencing.batch",  "white.british", "genetic.eur",
    "genetic.eur.oct2021", "genetic.eur.no.fin.oct2021",
    cts_phenotypes, binary_phenotypes)

dt <- cbind(dt[, ..filter_to], model.matrix(~ ukbb.centre + sex + genotyping.array + sequencing.batch, dt))
dt[, "(Intercept)":= NULL]
filter_to <- names(dt)
filter_to_cts <- setdiff(filter_to, binary_phenotypes)
filter_to_binary <- setdiff(filter_to, cts_phenotypes)

dt_cts <- dt[, ..filter_to_cts]
dt_binary <- dt[, ..filter_to_binary]

# Parentheses do not play nice with SAIGE.
names(dt_cts) <- gsub("[\\(\\)]", "", names(dt_cts))
names(dt_cts) <- gsub("\\-", "_", names(dt_cts))
names(dt_binary) <- gsub("[\\(\\)]", "", names(dt_binary))
names(dt_binary) <- gsub("\\-", "_", names(dt_binary))

fwrite(dt_cts, file=cts_output, sep='\t')
fwrite(dt_binary, file=binary_output, sep='\t')

system(paste("bgzip", cts_output))
system(paste("bgzip", binary_output))

# We also need to ensure that all samples with phenotype information are present in the vcf.
# We can do this by making use of the savvy library in SAIGE.

vcfFile <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr21.vcf.bgz")
vcfFileIndex <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr21.vcf.bgz.csi")
vcfField <- "GT"
isVariant <- setvcfDosageMatrix(vcfFile, vcfFileIndex, vcfField)
sampleListinvcf <- data.table(ID = getSampleIDlist_vcfMatrix())
setkey(sampleListinvcf, "ID")

dt_cts <- merge(sampleListinvcf, dt_cts)
dt_binary <- merge(sampleListinvcf, dt_binary)

# Remove all the withdrawn samples
withdrawn <- fread("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210809.csv") %>% transmute(ID = V1)

dt_cts <- dt_cts[which(!dt_cts$ID %in% withdrawn$ID), ]
dt_binary <- dt_binary[which(!dt_binary$ID %in% withdrawn$ID), ]

fwrite(dt_cts, file=cts_filtered_output, sep='\t')
fwrite(dt_binary, file=binary_filtered_output, sep='\t')

system(paste("bgzip", cts_filtered_output))
system(paste("bgzip", binary_filtered_output))

# Create a final two files that in addition is filtered to the collection of samples that have IMPUTED genotype information available.
dt_cts <- fread(cmd = paste("zcat", paste0(cts_filtered_output, ".gz")), key="ID")
dt_binary <- fread(cmd = paste("zcat", paste0(binary_filtered_output, ".gz")), key="ID")

imputed_EUR_fam <- fread(
    paste0("/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211026_long_ukb_wes_",
        TRANCHE, "_sparse_autosomes.fam")
    ) %>% transmute(ID = V2)
imputed_EUR_fam <- data.table(imputed_EUR_fam)
setkey(imputed_EUR_fam, "ID")
dt_cts <- merge(dt_cts, imputed_EUR_fam)
dt_binary <- merge(dt_binary, imputed_EUR_fam)

fwrite(dt_cts, file=cts_filtered_output_imp, sep='\t')
fwrite(dt_binary, file=binary_filtered_output_imp, sep='\t')

system(paste("bgzip", cts_filtered_output_imp))
system(paste("bgzip", binary_filtered_output_imp))

