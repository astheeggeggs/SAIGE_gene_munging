library(data.table)
library(dplyr)
library(readxl)

reformat_clinical <- function(primary_care_read_codes = "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt") {
    dt <- fread(primary_care_read_codes)
    dt_v2 <- dt[read_2 != "",c('eid', 'read_2')]
    dt_v3 <- dt[read_3 != "",c('eid', 'read_3')]

    paste_unique <- function(v) { return(paste(unique(v), collapse=" ")) }

    dt_v2 <- dcast(dt_v2, eid ~ ., fun.aggregate=paste_unique, value.var="read_2")
    dt_v3 <- dcast(dt_v3, eid ~ ., fun.aggregate=paste_unique, value.var="read_3")

    # Name the columns
    setnames(dt_v2, ".", "read_v2_string")
    setnames(dt_v3, ".", "read_v3_string")

    # Pad each end with a space, so we can apply the same methods as below.
    dt_v2[, read_v2_string := gsub("^(.*)$", " \\1 ", read_v2_string)]
    dt_v3[, read_v3_string := gsub("^(.*)$", " \\1 ", read_v3_string)]

    # Set keys
    setkey(dt_v2, "eid")
    setkey(dt_v3, "eid")

    # Merge the result, for combination 
    dt <- merge(dt_v2, dt_v3, all=TRUE)

    return(dt)
}

curate_binary_phenotypes <- function(
    phenotype_list,
    # phenotype_file = "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv", # Old version
    phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
    primary_care_read_codes = "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt",
    create_primary_care = TRUE
    )
{

	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
		return(cols)
	}


	dt <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# Covariates that we want
	PCs <- c("22009")
	age <- c("21003")
	sex <- c("31")
	centre <- c("54")

	# Extract the relevant columns for each of the encodings
	ICD10s <- c("41202", "41204", "40006", "40001", "40002")
	ICD9s <- c("41203", "41205", "40013")
	
	OPCS4s <- c("41200", "41210")

	NI_non_cancer <- c("20002")
	NI_cancer <- c("20001")
	NI_operation <- c("20004")

	self_reported_diagnosis_by_doctor_CAD_STR = c("6150")
	self_reported_diagnosis_by_doctor_COPD = c("6152")

	cols <-  get_cols(
		c(
			PCs, age, sex, centre, # Covariates
			ICD10s,
			ICD9s,
			OPCS4s,
			NI_non_cancer,
			NI_cancer,
			NI_operation,
			self_reported_diagnosis_by_doctor_COPD,
			self_reported_diagnosis_by_doctor_CAD_STR
			),
		dt)

	select_cols <- rep("character", (length(cols) + 1))
	names(select_cols) <- c("eid", cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt <- fread(phenotype_file, na.strings=NULL, select=select_cols)

	ICD10_cols <-  get_cols(ICD10s, dt)
	ICD9_cols <-  get_cols(ICD9s, dt)
	OPCS4_cols <- get_cols(OPCS4s, dt)

	NI_non_cancer_cols <- get_cols(NI_non_cancer, dt)
	NI_cancer_cols <- get_cols(NI_cancer, dt)
	NI_operation_cols <- get_cols(NI_operation, dt)

	self_reported_diagnosis_by_doctor_CAD_STR_cols <- get_cols(self_reported_diagnosis_by_doctor_CAD_STR, dt)
	self_reported_diagnosis_by_doctor_COPD_cols <- get_cols(self_reported_diagnosis_by_doctor_COPD, dt)

	# Combine vector into a single string
	dt[, ICD10_string := do.call(paste,.SD), .SDcols=ICD10_cols]
	dt[, ICD9_string := do.call(paste,.SD), .SDcols=ICD9_cols]
	dt[, OPCS4_string := do.call(paste,.SD), .SDcols=OPCS4_cols]
	dt[, NI_non_cancer_string := do.call(paste,.SD), .SDcols=NI_non_cancer_cols]
	dt[, NI_cancer_string := do.call(paste,.SD), .SDcols=NI_cancer_cols]
	dt[, NI_operation_string := do.call(paste,.SD), .SDcols=NI_operation_cols]
	dt[, self_reported_diagnosis_by_doctor_CAD_STR_string := do.call(paste,.SD), .SDcols=self_reported_diagnosis_by_doctor_CAD_STR_cols]
	dt[, self_reported_diagnosis_by_doctor_COPD_string := do.call(paste,.SD), .SDcols=self_reported_diagnosis_by_doctor_COPD_cols]

	# Remove trailing whitespace
	dt[, ICD10_string := gsub("( )+", " ", ICD10_string)]
	dt[, ICD9_string := gsub("( )+", " ", ICD9_string)]
	dt[, OPCS4_string := gsub("( )+", " ", OPCS4_string)]
	dt[, NI_non_cancer_string := gsub("( )+", " ", NI_non_cancer_string)]
	dt[, NI_cancer_string := gsub("( )+", " ", NI_cancer_string)]
	dt[, NI_operation_string := gsub("( )+", " ", NI_operation_string)]
	dt[, self_reported_diagnosis_by_doctor_CAD_STR_string := gsub("( )+", " ", self_reported_diagnosis_by_doctor_CAD_STR_string)]
	dt[, self_reported_diagnosis_by_doctor_COPD_string := gsub("( )+", " ", self_reported_diagnosis_by_doctor_COPD_string)]

	# Pad start and end with space
	dt[, ICD10_string := gsub("^(.*)$", " \\1 ", ICD10_string)]
	dt[, ICD9_string := gsub("^(.*)$", " \\1 ", ICD9_string)]
	dt[, OPCS4_string := gsub("^(.*)$", " \\1 ", OPCS4_string)]
	dt[, NI_non_cancer_string := gsub("^(.*)$", " \\1 ", NI_non_cancer_string)]
	dt[, NI_cancer_string := gsub("^(.*)$", " \\1 ", NI_cancer_string)]
	dt[, NI_operation_string := gsub("^(.*)$", " \\1 ", NI_operation_string)]
	dt[, self_reported_diagnosis_by_doctor_CAD_STR_string := gsub("^(.*)$", " \\1 ", self_reported_diagnosis_by_doctor_CAD_STR_string)]
	dt[, self_reported_diagnosis_by_doctor_COPD_string := gsub("^(.*)$", " \\1 ", self_reported_diagnosis_by_doctor_COPD_string)]

    # If we want to include the primary care data, read in and format
    if (create_primary_care) {
        dt_clin <- reformat_clinical(primary_care_read_codes=primary_care_read_codes)
        dt_clin[, eid := as.character(eid)]
        setkey(dt_clin, "eid")
        setkey(dt, "eid")
        dt <- merge(dt, dt_clin, all=TRUE)
    }

	# grep out boolean phenotypes according to the rules
	for (phenotype in names(phenotype_list)) {
		cat(paste0(phenotype, "...\n"))
		combined <- paste0(phenotype, "_combined")
		dt[, (combined):=FALSE]
		for (i in 1:length(phenotype_list[[phenotype]][["codings"]])) 
		{	
			coding <- phenotype_list[[phenotype]][["codings"]][[i]]
			coding_name <- names(phenotype_list[[phenotype]][["codings"]])[i]
			cat(paste0(coding_name, "..."))
			dt_col <- paste0(coding_name, "_string")
			new_colname <- paste0(phenotype, "_", coding_name)
			grep_exp <- paste0(" ", paste(coding, collapse=" | "), " ")
			if (!(dt_col %in% names(dt))) {
				stop("Error: column not found in data.table")
			}
			dt[, (new_colname):=grepl(grep_exp, dt[[dt_col]])]
			cat(paste(sum(dt[[new_colname]]), "\n"))
			dt[, (combined):=ifelse((dt[[combined]] | dt[[new_colname]]), TRUE, FALSE)]
			cat(paste(sum(dt[[combined]]), "\n"))
		}

        if (create_primary_care) {
            cat("creating encoding including primary care information...\n")
            combined <- paste0(phenotype, "_combined")
            combined_with_read <- paste0(phenotype, "_combined_primary_care")
            dt[, (combined_with_read):=dt[[combined]]]
            if (length(phenotype_list[[phenotype]][["read_codings"]]) > 0) {
                for (i in 1:length(phenotype_list[[phenotype]][["read_codings"]])) {
                    coding <- phenotype_list[[phenotype]][["read_codings"]][[i]]
                    coding_name <- names(phenotype_list[[phenotype]][["read_codings"]])[i]
                    cat(paste0(coding_name, "..."))
                    dt_col <- paste0(coding_name, "_string")
                    new_colname <- paste0(phenotype, "_", coding_name)
                    coding <- gsub("\\.", "\\\\.", coding)
                    grep_exp <- paste0(" ", paste(coding, collapse=" | "), " ")
                    if (!(dt_col %in% names(dt))) {
                        stop("Error: column not found in data.table")
                    }
                    dt[, (new_colname):=grepl(grep_exp, dt[[dt_col]])]
                    cat(paste(sum(dt[[new_colname]]), "\n"))
                    dt[, (combined_with_read):=ifelse((dt[[combined_with_read]] | dt[[new_colname]]), TRUE, FALSE)]
                    cat(paste(sum(dt[[combined_with_read]]), "\n"))
                }
            }
        }
	}

	# Define the final phenotypes
	phenotypes_with_ctrl_excl <- grep("_ctrl_excl", names(phenotype_list), value=TRUE)

	for (phenotype in gsub("_ctrl_excl", "", phenotypes_with_ctrl_excl))
	{
		cat(paste0(phenotype, "...\n"))
		current_phenotypes <- grep(phenotype, names(phenotype_list), value=TRUE)
		if (length(current_phenotypes) == 1) {
			cat("skip this one, no exclusions...\n")
			next
		} else {
			case_phenotypes <- setdiff(current_phenotypes, paste0(phenotype, "_ctrl_excl"))
			for (case_phenotype in case_phenotypes) {
				cat(paste0(case_phenotype, "...\n"))
				new_colname <- paste0(case_phenotype, "_combined")
				dt[, (new_colname):=ifelse(
					(dt[[paste0(phenotype, "_ctrl_excl_combined")]] & !dt[[paste0(case_phenotype, "_combined")]]),
					NA, dt[[paste0(case_phenotype, "_combined")]])]
			}
		}
	}

    if (create_primary_care)
    {
        for (phenotype in gsub("_ctrl_excl", "", phenotypes_with_ctrl_excl))
        {
            cat(paste0(phenotype, "...\n"))
            current_phenotypes <- grep(phenotype, names(phenotype_list), value=TRUE)
            if (length(current_phenotypes) == 1) {
                cat("skip this one, no exclusions...\n")
                next
            } else {
                case_phenotypes <- setdiff(current_phenotypes, paste0(phenotype, "_ctrl_excl"))
                for (case_phenotype in case_phenotypes) {
                    cat(paste0(case_phenotype, "...\n"))
                    new_colname <- paste0(case_phenotype, "_combined_primary_care")
                    dt[, (new_colname):=ifelse(
                        (dt[[paste0(phenotype, "_ctrl_excl_combined_primary_care")]] & !dt[[paste0(case_phenotype, "_combined_primary_care")]]),
                        NA, dt[[paste0(case_phenotype, "_combined_primary_care")]])]
                }
            }
        }
    }

	PC_cols <-  get_cols(PCs, dt)
	sex_col <-  get_cols(sex, dt)[1]
	centre_col <- get_cols(centre, dt)[1]
	age_col <- get_cols(age, dt)[1]
	new_PC_cols <- gsub(".*\\.", "PC", PC_cols)
	setnames(dt, PC_cols, new_PC_cols)
	setnames(dt, c(sex_col, age_col, centre_col), c("sex", "age", "ukbb.centre"))
	cols_to_retain <- grep("combined", names(dt), value=TRUE)
	cols_to_retain <- c("eid", "sex", "age", "ukbb.centre", new_PC_cols,
		setdiff(cols_to_retain, grep("ctrl_excl", cols_to_retain, value=TRUE)))
	return(dt[, ..cols_to_retain])
}

curate_biomarker_phenotypes <- function(
    # phenotype_file = "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv", # Old version
    phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
    biomarker_file="/well/lindgren/UKBIOBANK/DATA/Biomarker_data/ukb27722.csv",
    biomarker_fields=biomarker_fields_dt,
    filter_before_fitting="/well/lindgren/UKBIOBANK/samvida/icp_phewas/eids_passed_QC_210526.txt") {

    get_cols <- function(codes, dt, na.filter=FALSE)
    {
    cols <- c()
    for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
    return(cols)
    }

    dt_pheno <- fread(phenotype_file, na.strings=NULL, nrow=1)
    dt_bio <- fread(biomarker_file, na.string=NULL, nrow=1)

    # Covariates that we want
    PCs <- c("22009")
    age <- c("21003")
    sex <- c("31")
    ancestry <- c("21000")
    centre <- c("54")
    batch <- c("22000")

    # Biomarker specific
    sampling_time_blood <- c("3166")
    sampling_time_urine <- c("20035")
    fasting_time_hrs <- c("74")
    assessment_date <- c("53")

    # Biomarkers
    bio_fields <- biomarker_fields[["field"]]
    bio_dates <- biomarker_fields[["date_field"]]
    bio_report <- biomarker_fields[["reportability_field"]]
    dilution <- "30897"

    # Standard covariates
    PC_cols <- get_cols(PCs, dt_pheno)
    sex_col <-  get_cols(sex, dt_pheno)[1]
    centre_col <- get_cols(centre, dt_pheno)[1]
    age_col <- get_cols(age, dt_pheno)[1]
    ancestry_col <- get_cols(ancestry, dt_pheno)[1]
    batch_col <- get_cols(batch, dt_pheno)[1]

    # Biomarker covariates
    sampling_blood_col <- get_cols(sampling_time_blood, dt_pheno)[1]
    sampling_urine_col <- get_cols(sampling_time_urine, dt_pheno)[1]
    fasting_time_col <- get_cols(fasting_time_hrs, dt_pheno)[1]
    assessment_date_col <- get_cols(assessment_date, dt_pheno)[1]

    # Biomarkers
    dilution_col <- get_cols(dilution, dt_bio)[1]
    biomarker_cols <-  c(get_cols(c(bio_fields, bio_dates, bio_report), dt_bio), dilution_col)
    # Filter down to the first instance
    biomarker_cols <- grep("*\\-0\\.0", biomarker_cols, value=TRUE)

    pheno_cols <- c("eid",
        c(PC_cols, sex_col, centre_col, age_col, ancestry_col, batch_col,
          sampling_blood_col, sampling_urine_col, fasting_time_col, assessment_date_col))
    biomarker_cols <- c("eid", biomarker_cols)

    # Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
    dt_pheno <- fread(phenotype_file, na.strings=NULL, select=pheno_cols, key='eid')
    dt_bio <- fread(biomarker_file, na.strings=NULL, select=biomarker_cols, key='eid')
    dt <- merge(dt_pheno, dt_bio, all=TRUE)

    # Rename all the columns
    PC_colnames <- gsub(".*\\.", "PC", PC_cols)
    setnames(dt, PC_cols, gsub(".*\\.", "PC", PC_colnames))
    # Standard covariates
    setnames(dt, c(sampling_blood_col, sampling_urine_col, fasting_time_col,
        assessment_date_col, dilution_col),
        c("sampling_time_blood", "sampling_time_urine",
          "fasting_time_hrs", "assessment_date", "sample_dilution_factor"))
    # Biomarker specific covariates
    setnames(dt, c(sex_col, centre_col, age_col, ancestry_col, batch_col),
        c("sex", "ukbb.centre", "age", "ancestry", "batch"))
    setnames(dt, paste0(biomarker_fields[["field"]], "-0.0"), biomarker_fields[["biomarker"]])
    to_remove <- which(is.na(biomarker_fields[["date_field"]]))
    setnames(dt, paste0(biomarker_fields[["date_field"]][-to_remove], "-0.0"),
        paste0(biomarker_fields[["biomarker"]][-to_remove], "_date"))
    setnames(dt, paste0(biomarker_fields[["reportability_field"]][-to_remove], "-0.0"),
        paste0(biomarker_fields[["biomarker"]][-to_remove], "_reportability"))

    # Define derived phenotypes
    age_bin_cuts <- seq(40, 71, by = 5)
    dt[, ukbb.centre := as.factor(ukbb.centre)]
    dt[, sex := as.factor(sex)]
    dt[, age_5yr := as.factor(cut(dt$age, age_bin_cuts, include.lowest = TRUE, right = FALSE))]
    dt[, age := as.factor(age)]
    dt[, ancestry := as.factor(ifelse(is.na(ancestry) | ancestry %in% c(0,-1,-3), NA, ancestry))]
    # Convert sampling time to nearest second in the day.
    dt[, sampling_time_blood := as.numeric(hms(gsub("[^ ]+ (.*)", "\\1", as.character(sampling_time_blood))))]
    dt[, sampling_time_urine := as.numeric(hms(gsub("[^ ]+ (.*)", "\\1", as.character(sampling_time_urine))))]
    # Quantize these times and convert to indicators
    dt[, sampling_time_blood_icosatile := as.factor(cut(sampling_time_blood,
            quantile(sampling_time_blood, probs = seq(0, 1, length.out = 21), type = 1, na.rm=TRUE),
            include.lowest = TRUE, right = FALSE, labels = paste0("q", 1:20)))]
    dt[, sampling_time_urine_icosatile := as.factor(cut(sampling_time_urine,
        quantile(sampling_time_urine, probs = seq(0, 1, length.out = 21), type = 1, na.rm=TRUE),
        include.lowest = TRUE, right = FALSE, labels = paste0("q", 1:20)))]
    # Quantize the sample dilution factor
    dt[, sample_dilution_factor_icosatile := as.factor(cut(sample_dilution_factor,
            quantile(sample_dilution_factor, probs = seq(0, 1, length.out = 21), na.rm=TRUE),
            include.lowest = TRUE, labels = paste0("q", 1:20)))]
    # Fix tails in fasting time
    dt[, fasting_time_hrs := as.factor(
        ifelse(fasting_time_hrs == 0, 1,
            ifelse(fasting_time_hrs > 18, 19, fasting_time_hrs)))]
    dt[, assessment_month := substr(assessment_date, 1, 7)]
    dt[, assessment_month := as.factor(
        ifelse(grepl("^2006-", assessment_month), "2006-01",
            ifelse(assessment_month %in%  c("2010-08", "2010-09", "2010-10"), "2010-08",
                assessment_month)))]
    dt[, genotyping.array := as.factor(ifelse(batch < 0, "BiLEVE", "Axiom"))]
    dt[, batch := as.factor(batch)]
    # Determine the day of the week
    # Send the non-reportable values to NA
    wday_factor <- function(x) {
        days <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
        return(as.factor(days[wday(x)]))
    }
    dt[, (paste0(biomarker_fields[["biomarker"]][-to_remove], "_day")) := lapply(.SD, "wday_factor"),
        .SDcols = paste0(biomarker_fields[["biomarker"]][-to_remove], "_date")]
    
    # Log the biomarker values
    dt[, (biomarker_fields[["biomarker"]]) := lapply(.SD, "log"), .SDcols = biomarker_fields[["biomarker"]]]
    
    # Filter to the collection of samples that Samvida used.
    dt_filter <- fread(filter_before_fitting, header=FALSE, col.names="eid", key="eid")
    dt_to_fit <- merge(dt, dt_filter)

    dt_EUR_no_FIN_new <- fread("/well/lindgren/dpalmer/ukb_get_EUR/data/final_EUR_list.tsv",
        header=FALSE, col.names="eid", key="eid")
    dt <- merge(dt, dt_EUR_no_FIN_new)
    # Run the linear models
    for (biomarker in biomarker_fields[["biomarker"]]) {
        cat(paste0("fitting biomarker: ", biomarker, "..."))
        formula <- paste0(biomarker,
            " ~ sex + age + sex*age_5yr + ancestry + sex*ancestry + fasting_time_hrs + ",
            "sample_dilution_factor_icosatile + ukbb.centre + batch + assessment_month + ", 
            paste(paste0("PC", seq(1,40)), collapse = " + ")
            )
        if (grepl("Urine", biomarker)) {
            sampling_vars <- "sampling_time_urine_icosatile"
        } else {
            sampling_vars <- paste0("sampling_time_blood_icosatile + ", biomarker, "_day")
        }

        formula <- paste(formula, "+", sampling_vars)
        fit <- lm(as.formula(formula), dt_to_fit, na.action=na.exclude)
        cat("fitted\n")

        # Use the model to fill in the rest of the data in our non-Finnish European set.
        # There's a problem with this - some of the new data to fit will be at factor values not
        # present in the original fitted data.
        # Have to ensure that any levels that are not present in the data used to fit are
        # not present in the data that we predict in.
        # Split the formula into all the columns that are used.

        cols <- unique(strsplit(formula, split=" \\+ |\\*")[[1]][-1])
        dt_tmp <- copy(dt)

        # We run another check to avoid factor levels being present in the data to predict but 
        # not in the fitted data. This is because for low sample dilution factor values, many of 
        # the biomarkers do not have values.

        dt_tmp <- dt_tmp[is.na(dt_tmp[[biomarker]]) == FALSE]
        dt_tmp[, sample_dilution_factor_icosatile := droplevels(sample_dilution_factor_icosatile)]

        for (col in cols) {
            if (class(dt_tmp[[col]]) == "factor") {
                if (any(table(dt_to_fit[[col]]) == 0)) {
                    # Collection of levels that aren't present in the fit and so should be set to NA.
                    send_to_NA <- as.integer(names(table(dt_to_fit[[col]])[which(table(dt_to_fit[[col]]) == 0)]))
                    dt_tmp[[col]][which(dt_tmp[[col]] %in% send_to_NA)] <- NA
                }
            }
        }

        rownames(dt_tmp) <- dt_tmp[["eid"]]
        predicted <- predict(fit, dt_tmp, na.action=na.exclude)
        dt_predicted <- data.table(eid=as.integer(names(predicted)))
        dt_predicted[[paste0(biomarker, "_fitted")]] <- predicted
        setkey(dt_predicted, "eid")

        dt <- merge(dt, dt_predicted, all=TRUE)
        dt[, (paste0(biomarker, "_residual")) := dt[[biomarker]] - dt[[paste0(biomarker, "_fitted")]]]

        # Run the sex specific versions
        formula <- paste0(biomarker,
            " ~ age + ancestry + fasting_time_hrs + sample_dilution_factor_icosatile + ",
            "ukbb.centre + batch + assessment_month + ", 
            paste(paste0("PC", seq(1,40)), collapse = " + ")
            )
        if (grepl("Urine", biomarker)) {
            sampling_vars <- "sampling_time_urine_icosatile"
        } else {
            sampling_vars <- paste0("sampling_time_blood_icosatile + ", biomarker, "_day")
        }

        formula <- paste(formula, "+", sampling_vars)
        sex_list <- list()
        sex_list$F <- 0
        sex_list$M <- 1

        sex_fit <- list()
        for (sex_char in c("F", "M")) {
            sex_fit[[sex_char]] <- lm(as.formula(formula), dt_to_fit[sex == sex_list[[sex_char]]], na.action=na.exclude)
            cat(paste0("fitted ", sex_char, "...\n"))
        }

        cols <- unique(strsplit(formula, split=" \\+ |\\*")[[1]][-1])
        dt_tmp <- copy(dt)

        # We run another check to avoid factor levels being present in the data to predict but 
        # not in the fitted data. This is because for low sample dilution factor values, many of 
        # the biomarkers do not have values.

        dt_tmp <- dt_tmp[is.na(dt_tmp[[biomarker]]) == FALSE]
        dt_tmp_sex <- list()
        
        for (sex_char in c("F", "M"))
        {
            dt_tmp_sex[[sex_char]] <- dt_tmp[sex==sex_list[[sex_char]]]
            dt_tmp_sex[[sex_char]][, sample_dilution_factor_icosatile := droplevels(sample_dilution_factor_icosatile)]

            for (col in cols) {
                if (class(dt_tmp[[col]]) == "factor")
                {
                    table_sex <- table(dt_to_fit[sex == sex_list[[sex_char]]][[col]])
                    if (any(table_sex == 0))
                    {
                        # Collection of levels that aren't present in the fit and so should be set to NA.
                        send_to_NA <- as.integer(names(table_sex[which(table_sex == 0)]))
                        dt_tmp_sex[[sex_char]][[col]][which(dt_tmp_sex[[sex_char]][[col]] %in% send_to_NA)] <- NA
                    }
                }
            }
        }

        for (sex_char in c("F", "M"))
        {
            rownames(dt_tmp_sex[[sex_char]]) <- dt_tmp_sex[[sex_char]][["eid"]]
            predicted <- predict(sex_fit[[sex_char]], dt_tmp_sex[[sex_char]], na.action=na.exclude)
            dt_predicted <- data.table(eid=as.integer(names(predicted)))
            dt_predicted[[paste0(biomarker, "_", sex_char, "_fitted")]] <- predicted
            setkey(dt_predicted, "eid")
            dt <- merge(dt, dt_predicted, all=TRUE)
            dt[, (paste0(biomarker, "_", sex_char, "_residual")) := dt[[biomarker]] - dt[[paste0(biomarker, "_", sex_char, "_fitted")]]]
        }
    }
    return(dt)
}

# Use the excel spreadsheet to add the various read codes and ICD9 codes
extract_terms <- function(
    strings_to_map,
    mapping_file = "/well/lindgren/UKBIOBANK/dpalmer/ukb_primary_care_mappings/all_lkps_maps_v3.xlsx",
    mapping_file_sheet = "icd9_icd10",
    from_col="ICD10",
    to_col="ICD9")
{
    # Read in the excel spreadsheet converting x -> y
    # This includes all possible x terms and what they're mapped to.
    dt <- data.table(read_excel(mapping_file, sheet = mapping_file_sheet))
    return(setdiff(unique(dt[get(from_col) %in% strings_to_map][,get(to_col)]), "UNDEF"))
}

# Further curation of binary phenotypes to extract the read codes etc and include further
# putative cases
