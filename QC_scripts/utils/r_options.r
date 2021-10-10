TRANCHE <- '200k'

# 03_initial_sample_qc_plot.r
PLOTS <- '/well/lindgren/dpalmer/ukb_exome_qc/plots/'
# Define some thresholds 
T_sample_callRate <- 0.95
T_dpMean <- 19.5
T_gqMean <- 47.8

# 05_impute_sex_plot.r
T_impute_sex <- 0.6

# 07_ultra_rare_counts_plot.r
T_nURVSNP <- 100
T_nURVIndel <- 10

# 08_final_variant_qc_plot.r
T_variant_call_rate  <- 0.97
T_pHWE <- 1e-6

# # 14_final_variant_qc_filter.r
# VARIANT_LIST <- '../../variants_BipEx/14_final_qc.keep.variant_list'

# # 15_final_sample_qc_plot.r
# SAMPLE_BEFORE_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.before.samples.tsv'
# SAMPLE_AFTER_QC_FILE <- 'gsutil cat gs://dalio_bipolar_w1_w2_hail_02/data/samples/15_final_qc.after.samples.tsv'

# # 15_final_sample_qc_filter.r
# FINAL_SAMPLE_LIST <- '../../samples_BipEx/15_final_qc.keep.sample_list'
# FINAL_SAMPLE_LIST_REMOVE_SINGLETON_OUTLIERS <- '../../samples_BipEx/15_final_qc_remove_singleton_outliers.keep.sample_list'
