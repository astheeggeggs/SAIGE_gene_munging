library(data.table)
library(lubridate)

biomarker_fields_dt <- data.table(
  biomarker = c(
    "Alanine_aminotransferase",
    "Albumin",
    "Alkaline_Phosphatase",
    "Apolipoprotein_A",
    "Apolipoprotein_B",
    "Aspartate_aminotransferase",
    "C_reactive_Protein",
    "Calcium",
    "Cholesterol",
    "Creatinine_Serum",
    "Creatinine_Urine",
    "Cystatin_C_Serum",
    "Direct_Bilirubin",
    "Gamma_glutamyltransferase",
    "Glucose",
    "HbA1c",
    "HDL_Cholesterol",
    "IGF_1",
    "Lipoprotein_A",
    "Microalbumin_Urine",
    "Oestradiol",
    "Phosphate",
    "Potassium_Urine",
    "Rheumatoid_factor",
    "SHBG",
    "Sodium_Urine",
    "Testosterone",
    "Total_Bilirubin",
    "Total_Protein",
    "Triglyceride",
    "Urate",
    "Urea",
    "Vitamin_D"),
  field = c(
    "30620","30600","30610","30630","30640","30650","30710","30680","30690",
    "30700","30510","30720","30660","30730","30740","30750","30760","30770",
    "30790","30500","30800","30810","30520","30820","30830","30530","30850",
    "30840","30860","30870","30880","30670","30890"), 
  date_field = c(
    "30621","30601","30611","30631","30641","30651","30711","30681","30691",
    "30701",NA,     "30721","30661","30731","30741","30751","30761","30771",
    "30791",NA,     "30801","30811",NA,     "30821","30831",NA,     "30851",
    "30841","30861","30871","30881","30671","30891"),
  reportability_field = c(
    "30626","30606","30616","30636","30646","30656","30716","30686","30696",
    "30706",NA,     "30726","30666","30736","30746","30756","30766","30776",
    "30796",NA,     "30806","30816",NA,     "30826","30836",NA,     "30856",
    "30846","30866","30876","30886","30676","30896")
  )

# # Write the results.

# # Fit each of the models to all of the Non-Finnish European data.
# # Save the results.

# # Sanity check that we get the same numbers as Samvida.
# # I found a bug!
# # Now, get the new phenotypes that are required for the altered biomarkers
# # Let's just do it all again, but filtering the collection of samples that Samvida uses before running the regressions.

# samvida_covar <- fread("/well/lindgren/UKBIOBANK/samvida/icp_phewas/covariates_210527.txt")
# setkey(samvida_covar, "eid")

# covs_samvida <- c("age","age_5yr_dummy","age_dummy","ancestry","fasting_time_hrs","fasting_time_dummy",
#     "sdf","sdf_icosatile","UKB_assessment_centre","geno_batch","sampling_time_blood",
#     "sampling_time_blood_icosatile","sampling_time_urine","sampling_time_urine_icosatile",
#     "assessment_date","assessment_month_dummy","geno_array", paste0("PC", seq(1,40)))
# setnames(samvida_covar, covs_samvida, paste0(covs_samvida, "_samvida"))

# covs <- c("eid", "sex", "age", "age_5yr", "ancestry", "fasting_time_hrs",
#     "sample_dilution_factor_icosatile", "ukbb.centre", "batch", "assessment_month", 
#     paste0("PC", seq(1,40)), "sampling_time_urine_icosatile", "sampling_time_blood_icosatile")
# dt_covs <- dt[, ..covs]

# dt_both_covs <- merge(samvida_covar, dt_covs)

# # PCs
# for (i in seq(1,40)) {
#     print(max(abs(dt_both_covs[[paste0("PC", i)]] - dt_both_covs[[paste0("PC", i, "_samvida")]])))
# }

# all(dt_both_covs$age == dt_both_covs$age_samvida)
# all(dt_both_covs$age == dt_both_covs$age_dummy_samvida)
# all(dt_both_covs$age_5yr == dt_both_covs$age_5yr_dummy_samvida, na.rm=TRUE)
# all(dt_both_covs$ancestry == dt_both_covs$ancestry_samvida)
# all(dt_both_covs$UKB_assessment_centre_samvida == dt_both_covs$ukbb.centre)
# # Slight discrepancy between batch in the seq file and that in the phenotype file.
# levels(dt_both_covs$batch) <- list(
#     UKBiLEVEAX_b11 = -11, UKBiLEVEAX_b10 = -10, UKBiLEVEAX_b9 = -9, UKBiLEVEAX_b8 = -8,
#     UKBiLEVEAX_b7 = -7, UKBiLEVEAX_b6 = -6, UKBiLEVEAX_b5 = -5, UKBiLEVEAX_b4 = -4,
#     UKBiLEVEAX_b3 = -3, UKBiLEVEAX_b2 = -2, UKBiLEVEAX_b1 = -1,
#     Batch_b001 = 1, Batch_b002 = 2, Batch_b003 = 3, Batch_b004 = 4, Batch_b005 = 5, Batch_b006 = 6,
#     Batch_b007 = 7, Batch_b008 = 8, Batch_b009 = 9, Batch_b010 = 10, Batch_b011 = 11, Batch_b012 = 12,
#     Batch_b013 = 13, Batch_b014 = 14, Batch_b015 = 15, Batch_b016 = 16, Batch_b017 = 17, Batch_b018 = 18,
#     Batch_b019 = 19, Batch_b020 = 20, Batch_b021 = 21, Batch_b022 = 22, Batch_b023 = 23, Batch_b024 = 24,
#     Batch_b025 = 25, Batch_b026 = 26, Batch_b027 = 27, Batch_b028 = 28, Batch_b029 = 29, Batch_b030 = 30,
#     Batch_b031 = 31, Batch_b032 = 32, Batch_b033 = 33, Batch_b034 = 34, Batch_b035 = 35, Batch_b036 = 36,
#     Batch_b037 = 37, Batch_b038 = 38, Batch_b039 = 39, Batch_b040 = 40, Batch_b041 = 41, Batch_b042 = 42,
#     Batch_b043 = 43, Batch_b044 = 44, Batch_b045 = 45, Batch_b046 = 46, Batch_b047 = 47, Batch_b048 = 48,
#     Batch_b049 = 49, Batch_b050 = 50, Batch_b051 = 51, Batch_b052 = 52, Batch_b053 = 53, Batch_b054 = 54,
#     Batch_b055 = 55, Batch_b056 = 56, Batch_b057 = 57, Batch_b058 = 58, Batch_b059 = 59, Batch_b060 = 60,
#     Batch_b061 = 61, Batch_b062 = 62, Batch_b063 = 63, Batch_b064 = 64, Batch_b065 = 65, Batch_b066 = 66,
#     Batch_b067 = 67, Batch_b068 = 68, Batch_b069 = 69, Batch_b070 = 70, Batch_b071 = 71, Batch_b072 = 72,
#     Batch_b073 = 73, Batch_b074 = 74, Batch_b075 = 75, Batch_b076 = 76, Batch_b077 = 77, Batch_b078 = 78,
#     Batch_b079 = 79, Batch_b080 = 80, Batch_b081 = 81, Batch_b082 = 82, Batch_b083 = 83, Batch_b084 = 84,
#     Batch_b085 = 85, Batch_b086 = 86, Batch_b087 = 87, Batch_b088 = 88, Batch_b089 = 89, Batch_b090 = 90,
#     Batch_b091 = 91, Batch_b092 = 92, Batch_b093 = 93, Batch_b094 = 94, Batch_b095 = 95
#     )
# # 41 differences.
# all(dt_both_covs$batch == dt_both_covs$geno_batch_samvida)
# all(dt_both_covs$fasting_time_hrs == dt_both_covs$fasting_time_dummy_samvida)
# all(dt_both_covs$assessment_month == dt_both_covs$assessment_month_dummy_samvida)

# # Something to do with the edges.
# all(dt_both_covs$sampling_time_blood_icosatile == dt_both_covs$sampling_time_blood_icosatile_samvida) # Different
# all(dt_both_covs$sampling_time_urine_icosatile == dt_both_covs$sampling_time_urine_icosatile_samvida) # Different
# all(dt_both_covs$sample_dilution_factor_icosatile == dt_both_covs$sdf_icosatile_samvida)

# # Grab the biomarker information
# bio_list <- readRDS("/well/lindgren/UKBIOBANK/samvida/icp_phewas/biomarkers_QCd_210527.rds")
# for (biomarker in biomarker_fields_dt[["biomarker"]]) {
#     dt_tmp <- data.table(bio_list[[biomarker]])
#     setkey(dt_tmp, "eid")
#     if (grepl("Urine", biomarker)) {
#         cols <- c("eid", paste0(biomarker))
#     } else {
#         cols <- c("eid", paste0(biomarker, c("", "_date", "_reportability", "_day")))
#     }
#     dt_tmp <- merge(dt_tmp, dt[, ..cols])
#     print(max(abs(log(dt_tmp$value) - dt_tmp[[biomarker]])))
#     print(all(dt_tmp$reportability == dt_tmp[[paste0(biomarker, "_reportability")]]))
#     print(all(dt_tmp$date == dt_tmp[[paste0(biomarker, "_date")]]))
#     print(all(dt_tmp$assay_day == dt_tmp[[paste0(biomarker, "_day")]]))
# }

# for (biomarker in biomarker_fields_dt[["biomarker"]]) {
#     dt_bio_lm <- readRDS(paste0("/well/lindgren/UKBIOBANK/samvida/icp_phewas/results/adjusted_", biomarker, ".rds"))
#     dt_tmp <- data.table(dt_bio_lm$sex_combined$df)
#     dt_tmp[, eid := as.integer(eid)]
#     setkey(dt_tmp, "eid")
#     dt_tmp <- merge(dt_tmp, dt)
#     print(dt_tmp[["adj_value"]] - dt_tmp[[paste0(biomarker, "_residual")]])
# }

# # Found the error in
# # https://github.com/lindgrengroup/samvida_general/blob/main/UKB_scripts/biomarker_scripts/3_transform_adjust_values.R
# # Line 46. Samvida is fixing it now.
