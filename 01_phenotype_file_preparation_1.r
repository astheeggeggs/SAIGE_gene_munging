source("utils/extract_phenotypes.r") # functions to extract phenotypes
source("utils/phenotypes_plos_genetics.r") # creates a plos_genetics list
source("utils/phenotypes_biomarkers_rivas.r") # creates a biomarker_fields_dt data.table
source("utils/phenotypes_remaining_of_interest.r") # creates a remaining_phenotypes list

# Read in the data and extract the phenotypes with these encodings, and IDs.
phenotype_file <- "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv"

dt_plos_genetics <- curate_binary_phenotypes(plos_genetics, phenotype_file)
dt_biomarkers <- curate_biomarker_phenotypes()
dt_remaining <- curate_binary_phenotypes(remaining_phenotypes, phenotype_file)

# Obtain the T1D, T2D, and GDM phenotypes using utils/phenotypes_diabetes_ukb_cases_cole.r 
# script from Joanna Cole.
source("utils/phenotypes_diabetes_ukb_cases_cole.r")
cols_to_retain <- c("eid", "DM_T1D", "DM_T2D", "DM_GD", "DM_ctrl_excl")
dt_diabetes <- data.table(dt_diabetes)
setnames(dt_diabetes, c("f.eid", "probable_t1dm", "probable_t2dm", "possible_gdm", "dm_unlikely"), cols_to_retain)
dt_diabetes = dt_diabetes[, ..cols_to_retain]
# Set dm_unlikely to controls, and the rest to NA.
dt_diabetes[, DM_T1D := ifelse((DM_ctrl_excl == 0 & DM_T1D == 0), NA, DM_T1D)]
dt_diabetes[, DM_T2D := ifelse((DM_ctrl_excl == 0 & DM_T2D == 0), NA, DM_T2D)]
dt_diabetes[, DM_GD := ifelse((DM_ctrl_excl == 0 & DM_GD == 0), NA, DM_GD)]

dt_diabetes[, DM_T1D := ifelse(DM_T1D == 1, TRUE, FALSE)]
dt_diabetes[, DM_T2D := ifelse(DM_T2D == 1, TRUE, FALSE)]
dt_diabetes[, DM_GD := ifelse(DM_GD == 1, TRUE, FALSE)]

# Merge and save all of these for read in.
setkey(dt_plos_genetics, "eid")
setkey(dt_biomarkers, "eid")
setkey(dt_remaining, "eid")
setkey(dt_diabetes, "eid")

covariates <- setdiff(intersect(names(dt_plos_genetics), names(dt_remaining)), "eid")

to_retain_remaining <- setdiff(names(dt_remaining), covariates)
to_retain_diabetes <- c("eid", "DM_T1D", "DM_T2D", "DM_GD")

dt_remaining <- dt_remaining[, ..to_retain_remaining]
dt_diabetes <- dt_diabetes[, ..to_retain_diabetes]

to_retain_biomarkers <- c(biomarker_fields_dt$biomarker)
to_retain_biomarkers <- c(
	to_retain_biomarkers,
	paste0(to_retain_biomarkers, "_residual"),
	paste0(to_retain_biomarkers, "_M_residual"),
	paste0(to_retain_biomarkers, "_F_residual")
	)
to_retain_biomarkers <- c("eid", to_retain_biomarkers)

dt_biomarkers <- dt_biomarkers[, ..to_retain_biomarkers]

dt <- merge(merge(dt_plos_genetics, dt_remaining), merge(dt_diabetes, dt_biomarkers, all=TRUE))

# Next, we use Teresa's collection of curated ICD10 codes to extract remaining phenotypes of interest.
kallmann_syndrome - This one is carefully defined - grab from Teresa’s file.
PCOS4 - This one is carefully defined - grab from Teresa’s file.
# Does Teresa have a detailed def of this?
POI (Teresa has a detailed definition of this)

NI_non_cancer_list <- list(
	PCOS = c("1350") # 1350	polycystic ovaries/polycystic ovarian syndrome
)
