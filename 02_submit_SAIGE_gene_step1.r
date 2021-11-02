# R wrapper to submit all SAIGE gene step1 jobs
# Need to loop over all of the required phenotype columns in each of the files that Teresa passed to us.

# Phenotypes
TRANCHE <- '200k'

# Tables to write
folder <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/", TRANCHE)
cts_phenotype_file <- paste0(folder, '/UKBB_WES', TRANCHE , "_filtered_cts_phenotypes.tsv.gz")  
binary_phenotype_file <- paste0(folder, '/UKBB_WES', TRANCHE, "_filtered_binary_phenotypes.tsv.gz")
# cts_phenotype_file <- paste0(folder, '/UKBB_WES', TRANCHE , "_filtered_cts_phenotypes_with_imputed_genos.tsv.gz")  
# binary_phenotype_file <- paste0(folder, '/UKBB_WES', TRANCHE, "_filtered_binary_phenotypes_with_imputed_genos.tsv.gz")

cts_phenotypes <- c(
    "Visceral_adipose_tissue_volume_VAT",
    "Total_adipose_tissue_volume",
    "Abdominal_fat_ratio",
    "Liver_proton_density_fat_fraction_AMRA",
    "Alanine_aminotransferase",
    "Albumin",
    "Alkaline_phosphatase",
    "Apolipoprotein_A",
    "Apolipoprotein_B",
    "Aspartate_aminotransferase",
    "C_reactive_protein",
    "Calcium",
    "Cholesterol",
    "Creatinine",
    "Cystatin_C",
    "Direct_bilirubin",
    "Gamma_glutamyltransferase",
    "Glucose",
    "Glycated_haemoglobin_HbA1c",
    "HDL_cholesterol",
    "IGF_1",
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
    "Body_mass_index_BMI",
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
    "E230",
    "PCOS1",
    "PCOS2",
    "PCOS3",
    "PCOS4"
    )

cts_phenotypes <- c(
    # "Body_mass_index_BMI",
    # "Hip_circumference",
    # "Standing_height",
    "Waist_circumference" #,
    # "Body_fat_percentage"
    )

binary_phenotypes <- c(
    "coronary_artery_disease" #,
    # "T2D"
    )

# Covariates
covariate_cols <- c(
	"age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
	"sex2","sequencing.batch2"
	)
covars <- paste(covariate_cols, collapse=",")
ID_col <- "ID"
submission_script <- "02_SAIGE_gene_step1.sh"

# Quantitative traits
trait_type <- "quantitative"

# Loop and submit the jobs separately using qsub
for (phenotype in cts_phenotypes) {
	job_options <- paste(
	            "-v",
	            paste0(
	                "phenofile=\"", cts_phenotype_file, "\",",
	                "phenotype=\"", phenotype, "\",",
	                "covars=\"", covars, "\",",
	                "ID_col=\"", ID_col, "\",",
	                "trait_type=\"", trait_type, "\""
	               )
                )
	job_submission <- paste("qsub", job_options, submission_script)
	system(job_submission)
}

# Binary traits
trait_type <- "binary"

# Loop and submit the jobs separately using qsub
for (phenotype in binary_phenotypes) {
	job_options <- paste(
	            "-v",
	            paste0(
	                "phenofile=\"", binary_phenotype_file, "\",",
	                "phenotype=\"", phenotype, "\",",
	                "covars=\"", covars, "\",",
	                "ID_col=\"", ID_col, "\",",
	                "trait_type=\"", trait_type, "\""
	               )
                )
	job_submission <- paste("qsub", job_options, submission_script)
	system(job_submission)
}
