source("utils/extract_phenotypes.r") # functions to extract phenotypes
source("utils/phenotypes_plos_genetics.r") # creates a plos_genetics list
source("utils/phenotypes_biomarkers_rivas.r") # creates a biomarker_fields_dt data.table
source("utils/phenotypes_remaining_of_interest.r") # creates a remaining_phenotypes list

# Read in the data and extract the phenotypes with these encodings, and IDs.
# Old version
# phenotype_file <- "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv"
phenotype_file <- "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv"

# Create a data.table to use in the methods
methods_list <- c(plos_genetics, remaining_phenotypes)
methods_table <- data.table(
	outcome_code=character(),
	ICD10=character(), ICD9=character(),
	OPCS4=character(),
	NI_non_cancer=character(), NI_cancer=character(), NI_operation=character(),
	self_report=character(),
	read_v2=character(), read_v3=character()
	)

for (name in names(methods_list)) {
	methods_table <- rbind(
		methods_table,
		data.table(
			outcome_code = name,
			ICD10 = ifelse(!is.null(methods_list[[name]][["codings"]][["ICD10"]]), paste(methods_list[[name]][["codings"]][["ICD10"]], collapse="; "), ""),
			ICD9 = ifelse(!is.null(methods_list[[name]][["codings"]][["ICD9"]]), paste(methods_list[[name]][["codings"]][["ICD9"]], collapse="; "), ""),
			OPCS4 = ifelse(!is.null(methods_list[[name]][["codings"]][["OPCS4"]]), paste(methods_list[[name]][["codings"]][["OPCS4"]], collapse="; "), ""),
			NI_non_cancer = ifelse(!is.null(methods_list[[name]][["codings"]][["NI_non_cancer"]]), paste(methods_list[[name]][["codings"]][["NI_non_cancer"]], collapse="; "), ""),
			NI_cancer = ifelse(!is.null(methods_list[[name]][["codings"]][["NI_cancer"]]), paste(methods_list[[name]][["codings"]][["NI_cancer"]], collapse="; "), ""),
			NI_operation = ifelse(!is.null(methods_list[[name]][["codings"]][["NI_operation"]]), paste(methods_list[[name]][["codings"]][["NI_operation"]], collapse="; "), ""),
			self_report = ifelse(any(grepl("self_reported", names(methods_list[[name]][["codings"]]))), paste(methods_list[[name]][["codings"]][grep("self_reported", names(methods_list[[name]][["codings"]]))], collapse="; "), ""),
			read_v2 = ifelse(!is.null(methods_list[[name]][["read_codings"]][["read_v2"]]), paste(methods_list[[name]][["read_codings"]][["read_v2"]], collapse="; "), ""),
			read_v3 = ifelse(!is.null(methods_list[[name]][["read_codings"]][["read_v3"]]), paste(methods_list[[name]][["read_codings"]][["read_v3"]], collapse="; "), "")
		)
	)
}

fwrite(methods_table, file="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/ukb_binary_phenotype_definitions.tsv", sep="\t")

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
dt_biomarkers[, eid := as.character(eid)]
dt_diabetes[, eid := as.character(eid)]

# Merge and save all of these for read in.
setkey(dt_plos_genetics, "eid")
setkey(dt_biomarkers, "eid")
setkey(dt_remaining, "eid")
setkey(dt_diabetes, "eid")

dt_bin <- merge(merge(dt_plos_genetics, dt_remaining), dt_diabetes, all=TRUE)
source("utils/phenotypes_cts_traits.r")
dt_cts[, eid := as.character(eid)]
dt_cts <- merge(dt_biomarkers, dt_cts, all=TRUE)

fwrite(dt_bin, file = "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/curated_phenotypes_binary.tsv", sep='\t')
fwrite(dt_cts, file = "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/curated_phenotypes_cts.tsv", sep='\t')

# bgzip the resultant .tsv file.

# Next, we use Teresa's collection of curated ICD10 codes to extract remaining phenotypes of interest.
# kallmann_syndrome - This one is carefully defined - grab from Teresa’s file.
# PCOS4 - This one is carefully defined - grab from Teresa’s file.
# Also, infertility traits aren't quite right - need to consider male and female specific traits separately.
# Need to account for statins in three of the biomarkers too - follow the Nature genetics UKB biomarkers paper.

