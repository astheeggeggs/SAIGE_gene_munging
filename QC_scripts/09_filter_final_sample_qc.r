library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")
source("utils/helpers.r")

SAMPLE_AFTER_COMBINED_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc.after.samples.tsv')

dt_after <- fread(SAMPLE_AFTER_COMBINED_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE) %>%
  mutate(phase='After Variant QC')

dt_keep <- dt_after[, 's']
print(paste0("Started with: ", nrow(dt_keep), " samples"))

n_SDs <- 4

# r_ti_tv
dt_keep_ti_tv <- group_by(dt_after, sequencing_batch) %>%
  summarise(mean = mean(sample_qc.r_ti_tv), sd = sd(sample_qc.r_ti_tv)) %>%
  inner_join(dt_after, by='sequencing_batch') %>%
  filter(
    ((sample_qc.r_ti_tv >= mean - n_SDs*sd) & 
     (sample_qc.r_ti_tv <= mean + n_SDs*sd)) | 
    is.na(sd)
    )

dt_keep <- dt_keep_ti_tv %>% inner_join(dt_keep, by='s')
print(paste0("Remove Ti/Tv outliers: ", nrow(dt_keep), " samples remain"))

# r_het_hom_var
dt_keep_het_hom_var <- group_by(dt_after, sequencing_batch) %>%
  summarise(mean = mean(sample_qc.r_het_hom_var), sd = sd(sample_qc.r_het_hom_var)) %>%
  inner_join(dt_after, by='sequencing_batch') %>%
  filter(
    ((sample_qc.r_het_hom_var >= mean - n_SDs*sd) & 
     (sample_qc.r_het_hom_var <= mean + n_SDs*sd)) |
    is.na(sd)
    )

dt_keep <- dt_keep_het_hom_var %>% inner_join(dt_keep, by='s')
print(paste0("Remove Het/HomVar outliers: ", nrow(dt_keep), " samples remain"))

# r_insertion_deletion
dt_keep_insertion_deletion <- group_by(dt_after, sequencing_batch) %>%
  summarise(mean = mean(sample_qc.r_insertion_deletion), sd = sd(sample_qc.r_insertion_deletion)) %>%
  inner_join(dt_after, by='sequencing_batch') %>%
  filter(
    ((sample_qc.r_insertion_deletion >= mean - n_SDs*sd) &
     (sample_qc.r_insertion_deletion <= mean + n_SDs*sd)) | 
    is.na(sd)
    )

dt_keep <- dt_keep_insertion_deletion %>% inner_join(dt_keep, by='s')
print(paste0("Remove Ins/Del outliers: ", nrow(dt_keep), " samples remain"))

# n_singletons
T_n_singletons <- 175
dt_keep_n_singletons <- dt_after %>% filter((sample_qc.n_singleton <= T_n_singletons) | is.na(sample_qc.n_singleton))
dt_keep <- dt_keep_n_singletons %>% inner_join(dt_keep, by='s')

print(paste0("Remove n_singletons outliers: ", nrow(dt_keep), " samples remain"))

dt_final_sample_summary <- data.table(
  Filter = c(
    "Samples after population filters",
    paste0("Within batch Ti/Tv ratio outside ", n_SDs, " standard deviations"),
    paste0("Within batch Het/HomVar ratio outside ", n_SDs, " standard deviations"),
    paste0("Within batch Insertion/Deletion ratio outside ", n_SDs, " standard deviations"),
    paste0("n singletons >", T_n_singletons),
    "Samples after final sample filters"),
  Samples = c(
    nrow(dt_after),
    nrow(dt_after) - nrow(dt_keep_ti_tv),
    nrow(dt_after) - nrow(dt_keep_het_hom_var),
    nrow(dt_after) - nrow(dt_keep_insertion_deletion),
    nrow(dt_after) - nrow(dt_keep_n_singletons),
    nrow(dt_keep))
  )

FINAL_SAMPLE_SUMMARY <- paste0(
  '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
  '/ukb_wes_qc/data/samples/09_final_sample.summary.tsv')

FINAL_SAMPLE_LIST <- paste0(
  '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
  '/ukb_wes_qc/data/samples/09_final_qc.keep.sample_list')

fwrite(dt_final_sample_summary, file=FINAL_SAMPLE_SUMMARY, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

# write out
fwrite(dt_keep, file=FINAL_SAMPLE_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)

