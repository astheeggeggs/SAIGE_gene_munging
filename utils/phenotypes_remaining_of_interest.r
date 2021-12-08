# Take these directly from Teresa's file, and map to ICD9.
# Remove the phenotypes that we already have careful definitions from the plos genetics paper.

remaining_phenotypes <- list()

dt <- fread("/well/lindgren/UKBIOBANK/ferreira/convert_disease_codes/icd10_codes_v3.txt")
for (phenotype in unique(dt$disease)) {
	remaining_phenotypes[[phenotype]] <- list(codings = list(
		ICD10 = (dt %>% filter(disease == phenotype))$icd10
		)
	)
}

# Now remove all but the phenotypes that we want to retain (the rest have cleaner definitions in either Jenny or Teresa's work).
remaining_phenotypes <- remaining_phenotypes[c("ADHD", "Alzheimers_disease", "depression", "hypothalamic_amenorrhea", "autism",
           "Cirrhosis", "Crohns_disease", "IBD", "NASH", "psoriasis", "hematuria", "proteinuria",
           "oligomenorrhea", "Preeclampsia", "dementia", "ectopic_pregnancy", "polycystic_kidney_disease",
           "habitual_aborter", "POI")]

for (phenotype in names(remaining_phenotypes)) {
    remaining_phenotypes[[phenotype]]$outcome <- gsub("_", " ", phenotype)
    substr(remaining_phenotypes[[phenotype]]$outcome, 1, 1) <- toupper(substr(remaining_phenotypes[[phenotype]]$outcome, 1, 1))
}

rename_list <- list(
    ADHD = "ADHD",
    Alzheimers_disease = "AD",
    depression = "DEP",
    hypothalamic_amenorrhea = "HAM",
    autism = "AUT",
    Cirrhosis = "CIRR",
    Crohns_disease = "CD",
    IBD = "IBD",
    NASH = "NASH",
    psoriasis = "PSOR",
    hematuria = "HEM",
    proteinuria = "PRO",
    oligomenorrhea = "OLI",
    Preeclampsia = "PRE",
    dementia = "DEM",
    ectopic_pregnancy = "EP",
    polycystic_kidney_disease = "PKD",
    habitual_aborter = "HAB",
    POI = "POI",
    E230 = "HYP"
)

for (phenotype in names(remaining_phenotypes)) {
    names(remaining_phenotypes)[which(names(remaining_phenotypes) == phenotype)] <- rename_list[[phenotype]]
}

# Add some further ICD10 codes
remaining_phenotypes$PRO$codings$ICD10 <- c(remaining_phenotypes$PRO$codings$ICD10, "R80")
remaining_phenotypes$PRE$codings$ICD10 <- c(remaining_phenotypes$PRE$codings$ICD10, "O13", "016")
remaining_phenotypes$HAB$codings$ICD10 <- c(remaining_phenotypes$HAB$codings$ICD10, "N96")

for (phenotype in names(remaining_phenotypes)) {
    pheno_coding_tmp <- extract_terms(remaining_phenotypes[[phenotype]]$codings$ICD10)
    if (length(pheno_coding_tmp) > 0) {
        remaining_phenotypes[[phenotype]]$codings$ICD9 <- pheno_coding_tmp
    }
    
    pheno_coding_tmp_read_v3 <- extract_terms(remaining_phenotypes[[phenotype]]$codings$ICD10,
        mapping_file_sheet = "read_ctv3_icd10", from_col = "icd10_code", to_col = "read_code")
    pheno_coding_tmp_read_v2 <- extract_terms(remaining_phenotypes[[phenotype]]$codings$ICD10,
        mapping_file_sheet = "read_v2_icd10", from_col = "icd10_code", to_col = "read_code")

    if (length(pheno_coding_tmp_read_v2) > 0 | length(pheno_coding_tmp_read_v3) > 0) {
        remaining_phenotypes[[phenotype]]$read_codings <- list()
        if (length(pheno_coding_tmp_read_v2) > 0) {
            remaining_phenotypes[[phenotype]]$read_codings$read_v2 <- pheno_coding_tmp_read_v2
        }
        if (length(pheno_coding_tmp_read_v3) > 0) {
            remaining_phenotypes[[phenotype]]$read_codings$read_v3 <- pheno_coding_tmp_read_v3
        }
    }
}

NI_non_cancer_list <- list(
    DEM = c("1263"), # 1263 = dementia/alzheimers/cognitive impairment
    AD = c("1263"), # 1263 = dementia/alzheimers/cognitive impairment
    DEP = c("1269", "1531"), # 1531 = post-natal depression, 1286 = depression
    PSOR = c("1453"), # 1453 = psoriasis
    IBD = c("1461", "1462", "1463"), # 1461 = inflammatory bowel disease, 1462 = crohns disease, 1463 = ulcerative colitis
    CD = c("1462"), # 1462 = crohns disease
    CIRR = c("1158", "1604", "1506"), # 1158 = liver failure/cirrhosis, 1604 = alcoholic liver disease/alcoholic cirrhosis, 1506 = primary biliary cirrhosis
    EP = c("1558"), # 1558 = ectopic pregnancy
    PRE = c("1073"), # 1073 = gestational hypertension/pre-eclampsia
    # HAB = c("1559"), # 1559 = miscarriage
    POI = c("1551"), # 1551 = ovarian problem
    PKD = c("1427"), # 1427 = polycystic kidney
    HAM = c("1430", "1237"), # 1430 = hypopituitarism, 1237 = disorder of pituitary gland
    HYP = c("1430") # 1430 = hypopituitarism
)

for (phenotype in names(NI_non_cancer_list)) {
    remaining_phenotypes[[phenotype]]$codings$NI_non_cancer <- NI_non_cancer_list[[phenotype]]
}