# Extract the ICD encodings according to pLoS genetics paper
# https://journals.plos.org/plosgenetics/article/authors?id=10.1371/journal.pgen.1008405 Supplementary materials, Table B.

plos_genetics <- list()

plos_genetics$BC <- list(
	outcome = "Breast cancer",
	codings = list(
		ICD10 = c("C50", "C500", "C501", "C502", "C503", "C504", "C505", "C506", "C508", "C509", "Z853"),
		ICD9 = c("174", "1740", "1741", "1742", "1743", "1744", "1745", "1746", "1748", "1749", "V103"),
		NI_cancer = "1002"),
	description = "Codes for breast cancer, including personal history codes. Run in women only"
)

plos_genetics$CAD <- list(
	outcome = "CAD",
	codings = list(
		ICD10 = c("I20", "I200", "I201", "I208", "I209", "I21", "I210", "I211", "I212", "I213", "I214", "I219", "I21X", "I22", "I220", "I221", "I228", "I229", "I23", "I230", "I231", "I232", "I233", "I234", "I235", "I236", "I238", "I24", "I240", "I241", "I248", "I249", "I251", "I252", "I255", "I256", "I258", "I259"),
		ICD9 = c("410", "4109", "411", "4119", "412", "4129", "413", "4139", "4140", "4148", "4149"),
		OPCS4 = c("K40", "K401", "K402", "K403", "K404", "K408", "K409", "K41", "K411", "K412", "K413", "K414", "K418", "K419", "K42", "K421", "K422", "K423", "K424", "K428", "K429", "K43", "K431", "K432", "K433", "K434", "K438", "K439", "K44", "K441", "K442", "K448", "K449", "K45", "K451", "K452", "K453", "K454", "K455", "K456", "K458", "K459", "K46", "K461", "K462", "K463", "K464", "K465", "K468", "K469", "K49", "K491", "K492", "K493", "K494", "K498", "K499", "K501", "K75", "K751", "K752", "K753", "K754", "K758", "K759"),
		NI_non_cancer = c("1074", "1075"),
		NI_operation = c("1070", "1095", "1523"),
		self_reported_diagnosis_by_doctor_CAD_STR = c("1", "2")
		),
	description = "Codes for myocardial infarction, percutaneous transluminal coronary angioplasty, coronary artery bypass grafting, chronic ischemic heart disease and angina. Outcome definition from Nelson et al (29)", "SOFT CAD definition including angina. "
)

plos_genetics$CAD_ctrl_excl <- list(
	outcome = "CAD - control group exclusions",
	codings = list(
		ICD10 = c("I250", "I253", "I254"),
		ICD9 = c("4141")
		),
	description = "Participants were excluded from the CAD-control group if they had these codes pertaining to heart aneurysm and atherosclerotic cardiovascular disease. Outcome definition from Nelson et al (29)", "SOFT CAD definition including angina. "
)

plos_genetics$COPD <- list(
	outcome = "COPD",
	codings = list(
		ICD10 = c("J41", "J410", "J411", "J418", "J42", "J43", "J431", "J432", "J438", "J439", "J44", "J440", "J441", "J448", "J449"),
		ICD9 = c("491", "4910", "4911", "4912", "4918", "4919", "492", "4929"),
		NI_non_cancer = c("1112", "1113", "1472"),
		self_reported_diagnosis_by_doctor_COPD = "6"
		),
	description = "Codes for chronic bronchitis, emphysema, COPD, or complications specified from COPD"
)

plos_genetics$CLD <- list(
	outcome = "Chronic liver disease",
	codings = list(
		ICD10 = c("K702", "K703", "K704", "K717", "K721", "K74", "K740", "K741", "K742", "K743", "K744", "K745", "K746"),
		ICD9 = c("27103", "4562", "571", "5712", "5715", "57150", "57151", "57158", "57159", "5716"),
		NI_non_cancer = c("1604", "1158")
		),
	description = "Codes for fibrosis, sclerosis, cirrhosis of liver and liver failure, including if caused by alcohol, toxic liver disease, biliary cirrhosis, and other causes as well as codes defining complications specified as caused by these"
)

plos_genetics$CC <- list(
	outcome = "Colorectal Cancer",
	codings = list(
		ICD10 = c("C18", "C180", "C181", "C182", "C183", "C184", "C185", "C186", "C187", "C188", "C189", "C19", "C20"),
		ICD9 = c("153", "1530", "1531", "1532", "1533", "1534", "1535", "1536", "1537", "1538", "1539", "154", "1540", "1541"),
		NI_cancer = c("1020", "1022", "1023")
		),
	description = "Codes for cancers from caecum to rectum, including appendix"
)

plos_genetics$DEM <- list(
	outcome = "Dementia",
	codings = list(
		ICD10 = c("F00", "F000", "F001", "F002", "F009", "F01", "F010", "F011", "F012", "F013", "F018", "F019", "F03", "G30", "G300", "G301", "G308", "G309"),
		ICD9 = c("2900", "2901", "2902", "2903", "2904", "3310"),
		NI_non_cancer = "1263"
		),
	description = "Codes for dementia in Alzheimer's, vascular dementia, unspecified dementia, senile and presenile dementia"
)

plos_genetics$INF <- list(
	outcome = "Infertility",
	codings = list(
		ICD10 = c("N97", "N970", "N971", "N972", "N973", "N978", "N979", "N46"),
		ICD9 = c("628", "6280", "6281", "6282", "6283", "6284", "6288", "6289", "606", "6069"),
		NI_non_cancer = c("1403", "1404")
		),
	description = "Codes for male or female infertility of different anatomical origins. Sex-specific codes applied in relevant sex only"
)

plos_genetics$LC <- list(
	outcome = "Lung cancer",
	codings = list(
		ICD10 = c("C33", "C34", "C340", "C341", "C342", "C343", "C348", "C349", "Z851"),
		ICD9 = c("162", "1620", "1622", "1623", "1624", "1625", "1628", "1629", "V101"),
		NI_cancer = c("1001", "1027", "1028", "1080")
	),
	description = "Codes for cancers in trachea, bronchi, and lungs, including personal history codes"
)

plos_genetics$NAFLD <- list(
	outcome = "NAFLD",
	codings = list(
		ICD10 = "K760"
	),
	description = "Code for fatty liver disease"
)

plos_genetics$RF <- list(
	outcome = "Renal Failure",
	codings = list(
		ICD10 = c("N17", "N170", "N171", "N172", "N178", "N179", "N18", "N180", "N181", "N182", "N183", "N184", "N185", "N188", "N189", "N19", "I120", "I131", "I132", "Z992"),
		ICD9 = c("584", "5845", "5846", "5847", "5848", "5849", "585", "5859", "586", "5869"),
		OPCS4 = c("L746", "X40", "X401", "X402", "X403", "X404", "X405", "X406", "X407", "X408", "X409", "X41", "X411", "X412", "X418", "X419", "X42", "X421", "X428", "X429"),
		NI_non_cancer = c("1192", "1193", "1194"),
		NI_operation = c("1476", "1580", "1581", "1582")
		),
	description = "Codes for both acute, chronic and unspecified renal failure and chronic kidney disease. Also includes renal failure from hypertensive disease and various dialysis procedures in ICD-10 and OPCS-4 codes"
)

plos_genetics$RF_acute <- list(
	outcome = "Renal Failure - acute",
	codings = list(
		ICD10 = c("N17", "N170", "N171", "N172", "N178", "N179"),
		ICD9 = c("584", "5845", "5846", "5847", "5848", "5849")
		),
	description = "Codes that are specifically for acute renal failure"
)

plos_genetics$RF_chronic <- list(
	outcome = "Renal Failure - chronic",
	codings = list(
		ICD10 = c("N18", "N180", "N181", "N182", "N183", "N184", "N185", "N188", "N189"),
		ICD9 = c("585", "5859")
		),
	description = "Codes that are specifically for chronic kidney disease"
)

plos_genetics$RF_ctrl_excl <- list(
	outcome = "Renal Failure - control group exclusions",
	codings = list(
		ICD10 = c("N17", "N170", "N171", "N172", "N178", "N179", "N18", "N180", "N181", "N182", "N183", "N184", "N185", "N188", "N189", "N19", "I120", "I131", "I132", "Z992"),
		ICD9 = c("584", "5845", "5846", "5847", "5848", "5849", "585", "5859", "586", "5869"),
		OPCS4 = c("L746", "X40", "X401", "X402", "X403", "X404", "X405", "X406", "X407", "X408", "X409", "X41", "X411", "X412", "X418", "X419", "X42", "X421", "X428", "X429"),
		NI_non_cancer = c("1192", "1193", "1194"),
		NI_operation = c("1476", "1580", "1581", "1582")
		),
	description = "Participants were excluded from the renal failure control groups, including acute renal failure and chronic kidney disease, if they had these codes pertaining to renal failure"
)

plos_genetics$STR <- list(
	outcome = "Stroke",
	codings = list(
		ICD10 = c("I60", "I600", "I601", "I602", "I603", "I604", "I605", "I606", "I607", "I608", "I609", "I61", "I610", "I611", "I612", "I613", "I614", "I615", "I616", "I618", "I619", "I63", "I630", "I631", "I632", "I633", "I634", "I635", "I636", "I638", "I639", "I64"),
		ICD9 = c("430", "4309", "431", "4319", "434", "4340", "4341", "4349", "436", "4369"),
		NI_non_cancer = c("1081", "1086", "1491", "1583"),
		self_reported_diagnosis_by_doctor_CAD_STR = "3"
		),
	description = "Codes for subarachnoid and intracerebral hemorrhages and cerebral infarctions including cerebral thrombosis and embolism, and unspecified stroke. Does not included transient cerebral ischaemia but includes acute but ill-defined cerebrovascular disease"
)

plos_genetics$STR_ctrl_excl <- list(
	outcome = "Stroke - control group exclusions",
	codings = list(
		ICD10 = c("I60", "I600", "I601", "I602", "I603", "I604", "I605", "I606", "I607", "I608", "I609", "I61", "I610", "I611", "I612", "I613", "I614", "I615", "I616", "I618", "I619", "I63", "I630", "I631", "I632", "I633", "I634", "I635", "I636", "I638", "I639", "I64", "G45", "G450", "G451", "G452", "G453", "G454", "G458", "G459"),
		ICD9 = c("430", "4309", "431", "4319", "434", "4340", "4341", "4349", "436", "4369", "435", "4359"),
		NI_non_cancer = c("1081", "1082", "1086", "1491", "1583"),
		self_reported_diagnosis_by_doctor_CAD_STR = "3"
		),
	description = "Participants were excluded from the stroke control groups, including hemorrhagic and ischemic stroke, if they had these codes pertaining to stroke and transient ischemic attacks"
)

plos_genetics$STR_hem <- list(
	outcome = "Stroke - hemorrhagic",
	codings = list(
		ICD10 = c("I60", "I600", "I601", "I602", "I603", "I604", "I605", "I606", "I607", "I608", "I609", "I61", "I610", "I611", "I612", "I613", "I614", "I615", "I616", "I618", "I619"),
		ICD9 = c("430", "4309", "431", "4319"),
		NI_non_cancer = c("1086", "1491")
		),
	description = "Codes that denote hemorrhagic stroke"
)

plos_genetics$STR_isc <- list(
	outcome = "Stroke - ischemic",
	codings = list(
		ICD10 = c("I63", "I630", "I631", "I632", "I633", "I634", "I635", "I636", "I638", "I639"),
		ICD9 = c("434", "4340", "4341", "4349"),
		NI_non_cancer = "1583"
		),
	description = "Codes that denote ischemic stroke"
)

# Extract some more read codes for these plos genetics defined phenotypes. Note that we are currently not using these, but we could.
for (phenotype in names(plos_genetics)) {

    pheno_coding_tmp_read_v3 <- extract_terms(plos_genetics[[phenotype]]$codings$ICD10,
        mapping_file_sheet = "read_ctv3_icd10", from_col = "icd10_code", to_col = "read_code")
    pheno_coding_tmp_read_v2 <- extract_terms(plos_genetics[[phenotype]]$codings$ICD10,
        mapping_file_sheet = "read_v2_icd10", from_col = "icd10_code", to_col = "read_code")

    if (length(pheno_coding_tmp_read_v2) > 0 | length(pheno_coding_tmp_read_v3) > 0) {
        plos_genetics[[phenotype]]$read_codings <- list()
        if (length(pheno_coding_tmp_read_v2) > 0) {
            plos_genetics[[phenotype]]$read_codings$read_v2 <- pheno_coding_tmp_read_v2
        }
        if (length(pheno_coding_tmp_read_v3) > 0) {
            plos_genetics[[phenotype]]$read_codings$read_v3 <- pheno_coding_tmp_read_v3
        }
    }
}
