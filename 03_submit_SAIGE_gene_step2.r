library(data.table)
options(scipen=999)

# R wrapper to submit all SAIGE gene step2 jobs
# Need to loop over all of the required phenotype columns in each of the files that Teresa passed to us.

# Phenotypes

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
    # "breast_cancer",
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
    # "male_infertility",
    "oligomenorrhea",
    # "habitual_aborter",
    "female_infertility",
    # "ectopic_pregnancy",
    # "Preeclampsia",
    "GDM",
    # "intrahepatic_cholestasis_in_pregnancy",
    "polycystic_kidney_disease",
    "T2D",
    "T1D",
    "GDM2",
    "kallmann_syndrome",
    "E230"#,
    # "PCOS1",
    # "PCOS2",
    # "PCOS3",
    # "PCOS4"
    )


cts_phenotypes <- c(
    "Body_mass_index_BMI",
    "Hip_circumference",
    "Standing_height",
    "Waist_circumference",
    "Body_fat_percentage"
    )

binary_phenotypes <- c(
    "coronary_artery_disease",
    "T2D"
    )

# annotations <- c("pLoF", "missense|LC", "pLoF|missense|LC", "synonymous", "missense")
annotations <- c("pLoF", "synonymous")
submission_script <- "03_SAIGE_gene_step2.sh"
# Quantitative traits
# Loop and submit the jobs separately using qsub
# For each phenotype, wes want to run gene-based tests for distinct MAF cutoffs

minMAF_maxMAFforGroupTest <- data.table(
    minMAF = c(0, 0, 0, 0.01),
    maxMAFforGroupTest = c(0.0001, 0.001, 0.01, 0.5)
    )
vcfFile <- "/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr@.vcf.bgz"
vcfFileIndex <- "/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr@.vcf.bgz.csi"
phenotypeFolder <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k"

for (annotation in annotations) {

    if (annotation == "pLoF") {
        queue <- "short.qc"
    } else {
        queue <- "long.qc"
    }

    groupFile <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/ukb_wes_200k_filtered_chr@_", annotation, "_saige_gene.tsv.gz")
    for (phenotype in cts_phenotypes) {
        for (row in 1:nrow(minMAF_maxMAFforGroupTest)) {
            job_options <-paste(
                "-q", queue,
                "-v",
                paste0(
                    "vcfFile=\"", vcfFile, "\",",
                    "vcfFileIndex=\"", vcfFileIndex, "\",",
                    "phenotypeFolder=\"", phenotypeFolder, "\",",
                    "phenotype=\"", phenotype, "\",",
                    "groupFile=\"", groupFile, "\",",
                    "minMAF=\"", minMAF_maxMAFforGroupTest$minMAF[row], "\",",
                    "maxMAFforGroupTest=\"", minMAF_maxMAFforGroupTest$maxMAFforGroupTest[row], "\""
                    )
                )
            job_submission <- paste("qsub", job_options, submission_script)
            system(job_submission)
        }
    }
}

# Binary traits
# Loop and submit the jobs separately using qsub
# For each phenotype, we want to run gene-based tests for distinct MAF cutoffs

for (annotation in annotations) {
    
    if (annotation == "pLoF") {
        queue <- "short.qe"
    } else {
        queue <- "long.qf"
    }

    groupFile <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/ukb_wes_200k_filtered_chr@_", annotation, "_saige_gene.tsv.gz")
    for (phenotype in binary_phenotypes) {
    	for (row in 1:nrow(minMAF_maxMAFforGroupTest)) {
            job_options <- paste(
                "-q", queue,
                "-v",
                paste0(
                    "vcfFile=\"", vcfFile, "\",",
                    "vcfFileIndex=\"", vcfFileIndex, "\",",
                    "phenotypeFolder=\"", phenotypeFolder, "\",",
                    "phenotype=\"", phenotype, "\",",
                    "groupFile=\"", groupFile, "\",",
                    "minMAF=\"", minMAF_maxMAFforGroupTest$minMAF[row], "\",",
                    "maxMAFforGroupTest=\"", minMAF_maxMAFforGroupTest$maxMAFforGroupTest[row], "\""
                    )
                )
            job_submission <- paste("qsub", job_options, submission_script)
            system(job_submission)
        }
    }
}
