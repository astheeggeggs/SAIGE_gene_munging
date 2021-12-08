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
setnames(dt, c("f.eid", "probable_t1dm", "probable_t2dm", "possible_gdm", "dm_unlikely"), cols_to_retain)
dt[, ..(cols_to_retain)]
# Set dm_unlikely to controls, and the rest to NA.
dt[, DM_T1D := ]
dt[, DM_T2D := ]
dt[, DM_GD := ]

# Merge and save all of these for read in.

# Next, we use Teresa's collection of curated ICD10 codes to extract remaining phenotypes of interest.
kallmann_syndrome - This one is carefully defined - grab from Teresa’s file.
PCOS4 - This one is carefully defined - grab from Teresa’s file.
# Does Teresa have a detailed def of this?
POI (Teresa has a detailed definition of this)

NI_non_cancer_list <- list(
	PCOS = c("1350") # 1350	polycystic ovaries/polycystic ovarian syndrome
)
