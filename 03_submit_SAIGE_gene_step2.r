library(data.table)
options(scipen=999)

# R wrapper to submit all SAIGE gene step2 jobs
# Need to loop over all of the required phenotype columns in each of the files that Teresa passed to us.

TRANCHE <- "200k"

# Phenotypes
cts_phenotypes <- c(
    "Alanine_aminotransferase_residual",
    "Albumin_residual",
    "Alkaline_Phosphatase_residual",
    "Apolipoprotein_A_residual",
    "Apolipoprotein_B_residual",
    "Aspartate_aminotransferase_residual",
    "C_reactive_Protein_residual",
    "Calcium_residual",
    "Cholesterol_residual",
    "Creatinine_Serum_residual",
    "Creatinine_Urine_residual",
    "Cystatin_C_Serum_residual",
    "Direct_Bilirubin_residual",
    "Gamma_glutamyltransferase_residual",
    "Glucose_residual",
    "HbA1c_residual",
    "HDL_Cholesterol_residual",
    "IGF_1_residual",
    "Lipoprotein_A_residual",
    "Microalbumin_Urine_residual",
    "Oestradiol_residual",
    "Phosphate_residual",
    "Potassium_Urine_residual",
    "Rheumatoid_factor_residual",
    "SHBG_residual",
    "Sodium_Urine_residual",
    "Testosterone_residual",
    "Total_Bilirubin_residual",
    "Total_Protein_residual",
    "Triglyceride_residual",
    "Urate_residual",
    "Urea_residual",
    "Vitamin_D_residual",
    # "Alanine_aminotransferase_M_residual",
    # "Albumin_M_residual",
    # "Alkaline_Phosphatase_M_residual",
    # "Apolipoprotein_A_M_residual",
    # "Apolipoprotein_B_M_residual",
    # "Aspartate_aminotransferase_M_residual",
    # "C_reactive_Protein_M_residual",
    # "Calcium_M_residual",
    # "Cholesterol_M_residual",
    # "Creatinine_Serum_M_residual",
    # "Creatinine_Urine_M_residual",
    # "Cystatin_C_Serum_M_residual",
    # "Direct_Bilirubin_M_residual",
    # "Gamma_glutamyltransferase_M_residual",
    # "Glucose_M_residual",
    # "HbA1c_M_residual",
    # "HDL_Cholesterol_M_residual",
    # "IGF_1_M_residual",
    # "Lipoprotein_A_M_residual",
    # "Microalbumin_Urine_M_residual",
    # "Oestradiol_M_residual",
    # "Phosphate_M_residual",
    # "Potassium_Urine_M_residual",
    # "Rheumatoid_factor_M_residual",
    # "SHBG_M_residual",
    # "Sodium_Urine_M_residual",
    # "Testosterone_M_residual",
    # "Total_Bilirubin_M_residual",
    # "Total_Protein_M_residual",
    # "Triglyceride_M_residual",
    # "Urate_M_residual",
    # "Urea_M_residual",
    # "Vitamin_D_M_residual",
    # "Alanine_aminotransferase_F_residual",
    # "Albumin_F_residual",
    # "Alkaline_Phosphatase_F_residual",
    # "Apolipoprotein_A_F_residual",
    # "Apolipoprotein_B_F_residual",
    # "Aspartate_aminotransferase_F_residual",
    # "C_reactive_Protein_F_residual",
    # "Calcium_F_residual",
    # "Cholesterol_F_residual",
    # "Creatinine_Serum_F_residual",
    # "Creatinine_Urine_F_residual",
    # "Cystatin_C_Serum_F_residual",
    # "Direct_Bilirubin_F_residual",
    # "Gamma_glutamyltransferase_F_residual",
    # "Glucose_F_residual",
    # "HbA1c_F_residual",
    # "HDL_Cholesterol_F_residual",
    # "IGF_1_F_residual",
    # "Lipoprotein_A_F_residual",
    # "Microalbumin_Urine_F_residual",
    # "Oestradiol_F_residual",
    # "Phosphate_F_residual",
    # "Potassium_Urine_F_residual",
    # "Rheumatoid_factor_F_residual",
    # "SHBG_F_residual",
    # "Sodium_Urine_F_residual",
    # "Testosterone_F_residual",
    # "Total_Bilirubin_F_residual",
    # "Total_Protein_F_residual",
    # "Triglyceride_F_residual",
    # "Urate_F_residual",
    # "Urea_F_residual",
    # "Vitamin_D_F_residual",
    "BMI",
    "BMI_imp",
    "standing_height",
    "VAT",
    "AFR",
    "AMRA",
    "body_fat_percentage",
    "waist_circumference",
    "hip_circumference",
    "WHR",
    "WHR_adj_BMI" #,
    # "WHR_adj_BMI_M"
    # "WHR_adj_BMI_F"
    )  

binary_phenotypes <- c(
    "BC_combined",
    "BC_combined_primary_care",
    "CAD_combined",
    "CAD_combined_primary_care",
    "COPD_combined",
    "COPD_combined_primary_care",
    "CLD_combined",
    "CLD_combined_primary_care",
    "CC_combined",
    "CC_combined_primary_care",
    "DEM_combined",
    "DEM_combined_primary_care",
    "INF_combined",
    "INF_combined_primary_care",
    "LC_combined",
    "LC_combined_primary_care",
    "NAFLD_combined",
    "NAFLD_combined_primary_care",
    "RF_combined",
    "RF_combined_primary_care",
    "RF_acute_combined",
    "RF_acute_combined_primary_care",
    "RF_chronic_combined",
    "RF_chronic_combined_primary_care",
    "STR_combined",
    "STR_combined_primary_care",
    "STR_hem_combined",
    "STR_hem_combined_primary_care",
    "STR_isc_combined",
    "STR_isc_combined_primary_care",
    "ADHD_combined",
    "ADHD_combined_primary_care",
    "AD_combined",
    "AD_combined_primary_care",
    "DEP_combined",
    "DEP_combined_primary_care",
    "HAM_combined",
    "HAM_combined_primary_care",
    "AUT_combined",
    "AUT_combined_primary_care",
    "CIRR_combined",
    "CIRR_combined_primary_care",
    "CD_combined",
    "CD_combined_primary_care",
    "IBD_combined",
    "IBD_combined_primary_care",
    "NASH_combined",
    "NASH_combined_primary_care",
    "PSOR_combined",
    "PSOR_combined_primary_care",
    "HEM_combined",
    "HEM_combined_primary_care",
    "PRO_combined",
    "PRO_combined_primary_care",
    "OLI_combined",
    "OLI_combined_primary_care",
    "PRE_combined",
    "PRE_combined_primary_care",
    "EP_combined",
    "EP_combined_primary_care",
    "PKD_combined",
    "PKD_combined_primary_care",
    "HAB_combined",
    "HAB_combined_primary_care",
    "POI_combined",
    "POI_combined_primary_care",
    "HYP_combined",
    "HYP_combined_primary_care",
    "DM_T1D",
    "DM_T2D",
    "DM_GD")

annotations <- c('pLoF', 'damaging_missense|LC', 'pLoF|damaging_missense|LC', 'pLoF|damaging_missense',  'damaging_missense', 'other_missense', 'synonymous')
annotations <- c('pLoF', 'damaging_missense', 'other_missense', 'synonymous')
submission_script <- "03_SAIGE_gene_step2.sh"

# Quantitative traits
# Loop and submit the jobs separately using qsub
# For each phenotype, wes want to run gene-based tests for distinct MAF cutoffs

minMAF_maxMAFforGroupTest <- data.table(
    minMAF = c(0)#, 0, 0.01),
    maxMAFforGroupTest = c(0.01)#, 0.5, 0.5)
    )

vcfFile <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr@.vcf.bgz")
vcfFileIndex <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr@.vcf.bgz.csi")
phenotypeFolder <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/", TRANCHE)

for (annotation in annotations) {

    if (annotation %in% c('pLoF', 'damaging_missense|LC', 'pLoF|damaging_missense|LC', 'pLoF|damaging_missense',  'damaging_missense', 'synonymous')) {
        queue <- "short.qe@@short.hge"
    } else {
        queue <- "long.qf@@long.hgf"
    }

    annotation_filename <- gsub("\\|", "_", annotation)
    groupFile <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/", TRANCHE, "/SAIGE_gene_input/ukb_wes_", TRANCHE, "_filtered_chr@_", annotation_filename, "_saige_gene.tsv.gz")
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
    
    if (annotation %in% c('pLoF', 'damaging_missense|LC', 'pLoF|damaging_missense|LC', 'pLoF|damaging_missense',  'damaging_missense', 'synonymous')) {
        queue <- "short.qe@@short.hge"
    } else {
        queue <- "long.qf@@long.hgf"
    }

    annotation_filename <- gsub("\\|", "_", annotation)
    groupFile <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/", TRANCHE, "/SAIGE_gene_input/ukb_wes_", TRANCHE, "_filtered_chr@_", annotation_filename, "_saige_gene.tsv.gz")
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
