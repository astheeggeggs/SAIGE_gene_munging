library(dplyr)
library(data.table)

# module purge
# module load R/3.6.2-foss-2019b
# Rscript 03_initial_sample_qc_filter.r --tranche 200k
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--tranche", default='200k', help = "Which exome sequencing tranche?")
args <- parser$parse_args()

TRANCHE <- args$tranche

# Output
SAMPLE_LIST_INITIAL_QC = paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list')
SAMPLE_LIST_INITIAL_SAMPLE_COUNT = paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/03_sample_count.tsv')

# Run the plotting again to ensure that the thresholds are as in the plots.
source("03_initial_sample_qc_plot.r")
names(dt) <- gsub("sample_qc\\.", "", names(dt))

dt_out <- filter(dt, call_rate > T_sample_callRate) %>%
    filter(dp_stats.mean > T_dpMean) %>%
    filter(gq_stats.mean > T_gqMean)
dt_out <- dt_out %>% select(s)
print(nrow(dt_out))

fwrite(dt_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)

# Create the table too
dt_summary_count <- data.table(
    "Filter" = c("Initial samples in vcf",
                 paste0("Sample call rate < ", T_sample_callRate),
                 paste0("Mean DP < ", T_dpMean),
                 paste0("Mean GQ < ", T_gqMean),
                 "Samples after sample QC filters"),
    "Samples" = c(nrow(dt),
                nrow(filter(dt, call_rate <= T_sample_callRate)),
                nrow(filter(dt, dp_stats.mean <= T_dpMean)),
                nrow(filter(dt, gq_stats.mean <= T_gqMean)),
                nrow(dt_out)),
    "Batch 1" = c(nrow(dt %>% filter(sequencing_batch == "Batch 1")),
			    nrow(filter(dt, call_rate <= T_sample_callRate) %>% filter(sequencing_batch == "Batch 1")),
				nrow(filter(dt, dp_stats.mean <= T_dpMean) %>% filter(sequencing_batch == "Batch 1")),
				nrow(filter(dt, gq_stats.mean <= T_gqMean) %>% filter(sequencing_batch == "Batch 1")),
				sum(dt_out$s %in% (dt %>% filter(sequencing_batch == "Batch 1"))$s)),
	"Batch 2" = c(nrow(dt %>% filter(sequencing_batch == "Batch 2")),
                nrow(filter(dt, call_rate <= T_sample_callRate) %>% filter(sequencing_batch == "Batch 2")),
                nrow(filter(dt, dp_stats.mean <= T_dpMean) %>% filter(sequencing_batch == "Batch 2")),
                nrow(filter(dt, gq_stats.mean <= T_gqMean) %>% filter(sequencing_batch == "Batch 2")),
                sum(dt_out$s %in% (dt %>% filter(sequencing_batch == "Batch 2"))$s))
    )

fwrite(dt_summary_count, file=SAMPLE_LIST_INITIAL_SAMPLE_COUNT , quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
