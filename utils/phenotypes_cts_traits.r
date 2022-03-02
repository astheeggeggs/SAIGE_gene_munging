library(data.table)
# Read in the data, and pull out the phenotypes that we care about

cts_phenotypes <- data.table(
	phenotype_id = c(
		"22001",
		"21001",
		"23104",
		"50",
		"22407",
		"22434",
		"22436",
		"23099",
		"48",
		"49"),
	phenotype = c(
		"sex",
		"BMI",
		"BMI_imp",
		"standing_height",
		"VAT",
		"AFR",
		"AMRA",
		"body_fat_percentage",
		"waist_circumference",
		"hip_circumference"),
	description = c(
		"Sex",
		"Body mass index (BMI)",
		"Body mass index (BMI), measured by impedance",
		"Standing height",
		"Visceral adipose tissue volume (VAT)",
		"Abdominal fat ratio	Abdominal composition",
		"10P Liver PDFF (proton density fat fraction) - Liver proton density fat fraction (PDFF) measured as the average PDFF in up to nine (at least three) regions of interest (ROI) in the liver placed while avoiding any inhomogeneities, major vessels and bile ducts. This field replaces 'Liver proton density fat fraction (AMRA)'",
		"Body fat percentage",
		"Waist circumference",
		"Hip circumference"),
		# file = c(
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/MRI_DEXA_phenotypes/ukb44287.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/MRI_DEXA_phenotypes/ukb44287.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/MRI_DEXA_phenotypes/ukb44287.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
		# "/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv"),
	file = c(
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb49917_AbdComp_updates.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb49917_AbdComp_updates.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb49917_AbdComp_updates.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
		"/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv")
	)

for (current_file in unique(cts_phenotypes$file)) {
	cts_phenotypes_tmp <- cts_phenotypes[file == current_file]
	dt <- fread(current_file, key="eid")
	cols_to_extract <- c('eid')
	col_names <- c('eid')
	for (i in 1:nrow(cts_phenotypes_tmp)) {
		pheno_tmp <- grep(paste0("^", cts_phenotypes_tmp$phenotype_id[i], "-"), names(dt), value=TRUE)
		print(cts_phenotypes_tmp$phenotype[i])
		if (length(pheno_tmp) > 0) {
			pheno_count <- rep(0, length(pheno_tmp))
			j <- 1
			for (pheno in pheno_tmp) {
				pheno_count[j] <- sum(!is.na(dt[, ..pheno]))
				j <- j+1
			}
			cols_to_extract <- c(cols_to_extract, pheno_tmp[which(pheno_count == max(pheno_count))[1]])
			col_names <- c(col_names, cts_phenotypes_tmp$phenotype[i])
		}
	}

	if (!exists("dt_cts")) {
		dt_cts <- dt[, ..cols_to_extract]
	} else {
		dt_cts <- merge(dt_cts, dt[, ..cols_to_extract])
	}
	setnames(dt_cts, cols_to_extract, col_names)
}

# Finally, determine the derived phenotypes
dt_cts[, WHR := (waist_circumference / hip_circumference)]

fit <- lm(as.formula("WHR ~ BMI + sex"), dt_cts, na.action=na.exclude)

rownames(dt_cts) <- dt_cts[["eid"]]
predicted <- predict(fit, dt_cts, na.action=na.exclude)
dt_predicted <- data.table(eid=as.integer(names(predicted)))
dt_predicted[["WHR_fitted"]] <- predicted
setkey(dt_predicted, "eid")

dt_cts <- merge(dt_cts, dt_predicted, all=TRUE)
dt_cts[, WHR_adj_BMI := WHR - WHR_fitted]
dt_cts[, WHR_fitted := NULL]

sex_list <- list()
sex_list$F <- 0
sex_list$M <- 1

sex_fit <- list()
for (sex_char in c("F", "M")) {
    sex_fit[[sex_char]] <- lm(as.formula("WHR ~ BMI"), dt_cts[sex == sex_list[[sex_char]]], na.action=na.exclude)
    cat(paste0("fitted ", sex_char, "...\n"))
    dt_cts_tmp <- copy(dt_cts[sex == sex_list[[sex_char]]])
    rownames(dt_cts_tmp) <- dt_cts_tmp[["eid"]]
    predicted <- predict(sex_fit[[sex_char]], dt_cts_tmp, na.action=na.exclude)
	dt_predicted <- data.table(eid=as.integer(names(predicted)))
	dt_predicted[["WHR_fitted"]] <- predicted
	setkey(dt_predicted, "eid")
	dt_cts <- merge(dt_cts, dt_predicted, all=TRUE)
	dt_cts[, paste0("WHR_adj_BMI_", sex_char) := WHR - WHR_fitted]
	dt_cts[, WHR_fitted := NULL]
}
