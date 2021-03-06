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

cts_filtered_output_imp <- paste0(output_folder, '/UKBB_WES', TRANCHE , "_filtered_cts_phenotypes_with_imputed_genos.tsv")  
binary_filtered_output_imp <- paste0(output_folder, '/UKBB_WES', TRANCHE, "_filtered_binary_phenotypes_with_imputed_genos.tsv")

# Tables to merge
# Quantitative phenotypes
folder <- "/well/lindgren/UKBIOBANK/ferreira/UKBB_phenotypes/"
files <- c(
	paste0("UKBB_WES", TRANCHE, "_abdominalcomposition_phenotypes.txt"),
	paste0("UKBB_WES", TRANCHE, "_biomarkers_phenotypes.txt"),
	paste0("UKBB_WES", TRANCHE, "_physicalmeasures_phenotypes.txt")
	)
files <- paste0(folder, files)

# Cols to drop in dt_tmp
drop_cols <- c("age", "age2", "age3", "Inferred.Gender", "white.british", "genetic.eur",
	"genotyping.array", "ukbb.centre", paste0("PC", seq(1,40)), "sex", "sequencing.batch") 

dt <- fread(files[1])
dt[, ID:=as.character(ID)]
setkey(dt, "ID")

for (file in files[2:length(files)]) {
	dt_tmp <- fread(file)
    dt_tmp[, ID:=as.character(ID)]
	setkey(dt_tmp, "ID")
	dt_tmp[, (drop_cols):=NULL]
	dt <- merge(dt, dt_tmp)
}

# Binary traits
dt_tmp <- fread(paste0(folder, "UKBB_IMPUTED500k_GenetEUR_BinaryTraits.txt"))
dt_tmp[, (drop_cols):= NULL]
dt_tmp[, "missing":= NULL]
dt_tmp[, ID:=as.character(ID)]
setkey(dt_tmp, "ID")

dt <- merge(dt, dt_tmp, all.x=TRUE)

# Filter to the full set of non inverse rank normalised phenotypes, which we will pass to SAIGE

cts_phenotypes <- c(
    "Visceral_adipose_tissue_volume_(VAT)",
    "Total_adipose_tissue_volume",
    "Abdominal_fat_ratio",
    "Liver_proton_density_fat_fraction_(AMRA)",
    "Alanine_aminotransferase",
    "Albumin",
    "Alkaline_phosphatase",
    "Apolipoprotein_A",
    "Apolipoprotein_B",
    "Aspartate_aminotransferase",
    "C-reactive_protein",
    "Calcium",
    "Cholesterol",
    "Creatinine",
    "Cystatin_C",
    "Direct_bilirubin",
    "Gamma_glutamyltransferase",
    "Glucose",
    "Glycated_haemoglobin_(HbA1c)",
    "HDL_cholesterol",
    "IGF-1",
    "LDL_direct",
    "Lipoprotein_A",
    "Oestradiol",
    "Phosphate",
    "Rheumatoid_factor",
    "SHBG",
    "Testosterone",
    "Total_bilirubin",
    "Total_protein",
    "Triglycerides",
    "Urate",
    "Urea",
    "Vitamin_D",
    "Body_mass_index_(BMI)",
    "Hip_circumference",
    "Standing_height",
    "Waist_circumference",
    "Body_fat_percentage"
    )

binary_phenotypes <- c(
    "colorectal_cancer",
    "Trachea_bronchus_lung_cancer",
    "breast_cancer",
    "hypothalamic_amenorrhea",
    "POI",
    "dementia",
    "Alzheimers_disease",
    "depression",
    "autism",
    "ADHD",
    "renal_failure",
    "coronary_artery_disease",
    "ischaemic_heart_disease",
    "stroke_hemorrhagic",
    "stroke",
    "ischaemic_stroke",
    "chronic_obstructive_pulmonary_disease",
    "Crohns_disease",
    "IBD",
    "Cirrhosis",
    "NASH",
    "NAFLD",
    "psoriasis",
    "hyperandrogenism",
    "hematuria",
    "proteinuria",
    "acute_renal_failure",
    "chronic_kidney_disease",
    "male_infertility",
    "oligomenorrhea",
    "habitual_aborter",
    "female_infertility",
    "ectopic_pregnancy",
    "Preeclampsia",
    "GDM",
    "intrahepatic_cholestasis_in_pregnancy",
    "polycystic_kidney_disease",
    "T2D",
    "T1D",
    "GDM2",
    "kallmann_syndrome",
    "E230", #? Hypopituitarism?
    "PCOS1",
    "PCOS2",
    "PCOS3",
    "PCOS4"
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

