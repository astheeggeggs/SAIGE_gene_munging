library(data.table)
library(ggplot2)

# Examine the curated phenotypes
cts_output <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k/UKBB_WES200k_cts_phenotypes.tsv"  
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

dt <- fread(cts_output)

pdf(file="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k/UKBB_WES200k_cts_phenotypes.pdf")
for (phenotype in cts_phenotypes) {
	
}

binary_output <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k/UKBB_WES200k_binary_phenotypes.pdf")
dt <- fread(binary_output)

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
