rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--tranche", default='200k', help = "Which exome sequencing tranche?")
args <- parser$parse_args()

TRANCHE <- args$tranche

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")
source("utils/helpers.r")

save_figures <- TRUE

CHR <- 1

# Input files
SAMPLE_BEFORE_QC_FILE <-  paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc_chr', CHR, '.before.samples.tsv.bgz')
SAMPLE_AFTER_QC_FILE <-  paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc_chr', CHR, '.after.samples.tsv.bgz')

# Output files
SAMPLE_BEFORE_COMBINED_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc.before.samples.tsv')
SAMPLE_AFTER_COMBINED_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc.after.samples.tsv')

read_and_key <- function(gz_filepath, starts='sample_qc', tmp=FALSE) {
    dt <- fread(cmd = paste('zcat', gz_filepath),
        stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(s, starts_with(starts))
    if (tmp) {
        names(dt) <- c('s', paste0('tmp_', names(dt)[-1]))
    }
    setkey(dt, 's')
    return(dt)
}

update_entries <- function(dt) {
    dt <- dt %>% mutate(
        # Counts
        sample_qc.n_called = sample_qc.n_called + tmp_sample_qc.n_called,
        sample_qc.n_not_called = sample_qc.n_not_called + tmp_sample_qc.n_not_called,
        sample_qc.n_filtered = sample_qc.n_filtered + tmp_sample_qc.n_filtered,
        sample_qc.n_hom_ref = sample_qc.n_hom_ref + tmp_sample_qc.n_hom_ref,
        sample_qc.n_het = sample_qc.n_het + tmp_sample_qc.n_het,
        sample_qc.n_hom_var = sample_qc.n_hom_var + tmp_sample_qc.n_hom_var,
        sample_qc.n_non_ref = sample_qc.n_non_ref + tmp_sample_qc.n_non_ref,
        sample_qc.n_singleton = sample_qc.n_singleton + tmp_sample_qc.n_singleton,
        sample_qc.n_snp = sample_qc.n_snp + tmp_sample_qc.n_snp,
        sample_qc.n_insertion = sample_qc.n_insertion + tmp_sample_qc.n_insertion,
        sample_qc.n_deletion = sample_qc.n_deletion + tmp_sample_qc.n_deletion,
        sample_qc.n_transition = sample_qc.n_transition + tmp_sample_qc.n_transition,
        sample_qc.n_transversion = sample_qc.n_transversion + tmp_sample_qc.n_transversion,
        sample_qc.n_star = sample_qc.n_star + tmp_sample_qc.n_star
        ) %>% select(c(s, starts_with('sample_qc')))
    setkey(dt, 's')
    return(dt)
}

add_ratios <- function(dt) {
    dt <- dt %>% mutate(
        # Ratios
        sample_qc.call_rate = sample_qc.n_called / (sample_qc.n_called + sample_qc.n_not_called + sample_qc.n_filtered),
        sample_qc.r_ti_tv = ifelse(sample_qc.n_transversion == 0, NA, sample_qc.n_transition / sample_qc.n_transversion),
        sample_qc.r_het_hom_var = ifelse(sample_qc.n_hom_var == 0, NA, sample_qc.n_het / sample_qc.n_hom_var),
        sample_qc.r_insertion_deletion = ifelse(sample_qc.n_deletion == 0, NA, sample_qc.n_insertion / sample_qc.n_deletion),
    ) %>% select(c('s', starts_with('sample_qc')))
    setkey(dt, 's')
    return(dt)
}

dt_before <- read_and_key(SAMPLE_BEFORE_QC_FILE)
dt_after <- read_and_key(SAMPLE_AFTER_QC_FILE)

for (CHR in c(seq(2,22), "X")) {
    # Input files
    cat(paste0("chromosome ", CHR, "\n"))
    SAMPLE_BEFORE_QC_FILE <-  paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
        '/ukb_wes_qc/data/samples/09_final_qc_chr', CHR, '.before.samples.tsv.bgz')
    SAMPLE_AFTER_QC_FILE <-  paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
        '/ukb_wes_qc/data/samples/09_final_qc_chr', CHR, '.after.samples.tsv.bgz')

    dt_before_tmp <- read_and_key(SAMPLE_BEFORE_QC_FILE, tmp=TRUE)
    dt_after_tmp <- read_and_key(SAMPLE_AFTER_QC_FILE, tmp=TRUE)

    dt_before  <- merge(dt_before, dt_before_tmp)
    dt_after <- merge(dt_after, dt_after_tmp)

    dt_before <- update_entries(dt_before)
    dt_after <- update_entries(dt_after)
}

dt_before <- add_ratios(dt_before)
dt_after <- add_ratios(dt_after)

# Create plots
save_figures <- TRUE
dt_pheno <- create_pheno_dt(TRANCHE)

dt_before <- merge(dt_before, dt_pheno)
dt_after <- merge(dt_after, dt_pheno)

fwrite(dt_before, file=SAMPLE_BEFORE_COMBINED_QC_FILE, sep='\t')
fwrite(dt_after, file=SAMPLE_AFTER_COMBINED_QC_FILE, sep='\t')

dt_before <- fread(SAMPLE_BEFORE_COMBINED_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE) %>% 
    mutate(phase='Before Variant QC')

dt_after <- fread(SAMPLE_AFTER_COMBINED_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE) %>%
    mutate(phase='After Variant QC')

dt <- bind_rows(dt_before, dt_after) %>%
    mutate(phase=factor(phase, levels=c('Before Variant QC', 'After Variant QC')))

## Split by Location.
# Colour by sequencing batch
# Number of singletons.

y_labels <- ''
alpha <- 0.8

# By UKBB Centre
# Number of singletons
create_pretty_boxplots(dt, aes(y=sample_qc.n_singleton, x=factor(ukbb_centre)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figures, file=paste0(PLOTS, '09_nSingletons_by_centre'), y_label=y_labels,
    alpha=alpha, height=140)

# rHetHomVar
create_pretty_boxplots(dt, aes(y=sample_qc.r_het_hom_var, x=factor(ukbb_centre)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figures, file=paste0(PLOTS, '09_rHetHomVar_by_centre'), y_label=y_labels,
    alpha=alpha, height=140)

# rInsertionDeletion
create_pretty_boxplots(dt, aes(y=sample_qc.r_insertion_deletion, x=factor(ukbb_centre)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figures, file=paste0(PLOTS, '09_rInsertionDeletion_by_centre'), y_label=y_labels,
    alpha=alpha, height=140)

# rTiTv
create_pretty_boxplots(dt, aes(y=sample_qc.r_ti_tv, x=factor(ukbb_centre)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figures, file=paste0(PLOTS, '09_rTiTv_by_centre'), n_ticks=5, y_label=y_labels,
    alpha=alpha, height=140)

# By Batch
# Number of singletons
create_pretty_boxplots(dt, aes(y=sample_qc.n_singleton, x=factor(sequencing_batch)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figures, file=paste0(PLOTS, '09_nSingletons_by_batch'), y_label=y_labels,
    alpha=alpha, height=80)

# rHetHomVar
create_pretty_boxplots(dt, aes(y=sample_qc.r_het_hom_var, x=factor(sequencing_batch)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figures, file=paste0(PLOTS, '09_rHetHomVar_by_batch'), y_label=y_labels,
    alpha=alpha, height=80)

# rInsertionDeletion
create_pretty_boxplots(dt, aes(y=sample_qc.r_insertion_deletion, x=factor(sequencing_batch)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figures, file=paste0(PLOTS, '09_rInsertionDeletion_by_batch'), y_label=y_labels,
    alpha=alpha, height=80)

# rTiTv
create_pretty_boxplots(dt, aes(y=sample_qc.r_ti_tv, x=factor(sequencing_batch)),
    aes(color=factor(sequencing_batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figures, file=paste0(PLOTS, '09_rTiTv_by_batch'), n_ticks=5, y_label=y_labels,
    alpha=alpha, height=80)

q <- function(x, n_sd=4) {
    y <- subset(x, x < (mean(x) - n_sd * sd(x))| x > (mean(x) + n_sd * sd(x)))
    print(y)
    if (length(y) == 0) y <- NA
    return(y)
}

