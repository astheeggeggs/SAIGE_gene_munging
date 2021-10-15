library(data.table)
library(dplyr)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--tranche", default='200k', help = "Which exome sequencing tranche?")
args <- parser$parse_args()

TRANCHE <- args$tranche

source("08_final_variant_qc_plot.r")

FINAL_VARIANT_LIST <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/variants/08_final_qc.keep.variant_list')
VARIANT_SUMMARY <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/variants/08_variant_count.tsv')

dt <- fread(COMBINED_VARIANT_QC_FILE, header=TRUE, sep='\t')
print(paste0("Initial number of variants: ", nrow(dt)))
current_variants <- nrow(dt)

dt_out <- dt %>% filter((qc.AF > 0) & (qc.AF < 1))
print(paste0("After removing invariant sites in the cleaned dataset: ", nrow(dt_out)))
print(nrow(dt_out) - current_variants)
current_variants <- nrow(dt_out)

dt_out <- dt_out %>% filter(qc.call_rate >= T_variant_call_rate)
print(paste0("After ensuring that the overall call rate is good: ", nrow(dt_out)))
print(nrow(dt_out) - current_variants)
current_variants <- nrow(dt_out)

dt_out <- dt_out %>% filter(qc.p_value_hwe > T_pHWE)
print(paste0("After ensuring that sites pass HWE p-value threshold: ", nrow(dt_out)))
print(nrow(dt_out) - current_variants)
current_variants <- nrow(dt_out)

dt <- fread(cmd = paste0('zcat ', COMBINED_VARIANT_QC_FILE, '.bgz'), header=TRUE, sep='\t')
print("After removing invariant sites in the cleaned dataset: ")
print(nrow(dt) - nrow(dt %>% filter((qc.AF > 0) & (qc.AF < 1))))

print("After ensuring that the overall call rate is good: ")
print(nrow(dt) - nrow(dt %>% filter(qc.call_rate >= T_variant_call_rate)))

print("After ensuring that sites pass HWE p-value threshold: " )
print(nrow(dt) - nrow(dt %>% filter(qc.p_value_hwe > T_pHWE)))

dt_final_variant_summary <- data.table(Filter = c("Variants after initial filter",
												  "Invariant sites after sample filters",
												  paste0("Overall variant call rate < ", T_variant_call_rate),
												  "Variants failing HWE filter",
												  "Variants after filters"),
									   Variants = c(nrow(dt),
									   			    nrow(dt %>% filter((qc.AF <= 0) | (qc.AF >= 1))),
									   			    nrow(dt %>% filter(qc.call_rate < T_variant_call_rate)),
									   			    nrow(dt %>% filter(qc.p_value_hwe <= T_pHWE)),
									   			    nrow(dt_out)))


fwrite(dt_final_variant_summary, file=VARIANT_SUMMARY, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
fwrite(dt_out %>% select("locus", "alleles"), file=FINAL_VARIANT_LIST, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
