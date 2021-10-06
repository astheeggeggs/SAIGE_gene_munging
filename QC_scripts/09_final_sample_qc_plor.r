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

save_figures <- TRUE

CHR <- 1

# Input files
SAMPLE_BEFORE_QC_FILE <-  paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc_chr', CHR, '.before.samples.tsv')
SAMPLE_AFTER_QC_FILE <-  paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc_chr', CHR, '.after.samples.tsv')

PHENOFILE_CTS <- paste0('/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/',
    TRANCHE, '/UKBB_WES200k_filtered_cts_phenotypes.tsv.gz')

# Output files
SAMPLE_BEFORE_COMBINED_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc.before.samples.tsv.bgz')
SAMPLE_AFTER_COMBINED_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/09_final_qc.after.samples.tsv.bgz')

read_and_key <- function(gz_filepath, starts_with='sample_qc', tmp=FALSE) {
    dt <- fread(cmd = paste('zcat', gz_filepath),
        stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(c(s, starts_with(starts))))
    if (tmp) {
        names(dt) <- c('s', paste0('tmp_', names(dt)[-1]))
    }
    setkey(dt, 's')
    return(dt)
}

update_entries <- function(dt) {
    dt <- dt %>% mutate(
        # Counts
        sample_qc.n_singleton = sample_qc.n_singleton + tmp_sample_qc.n_singleton

        # DEV: need to include the ratios here!
        ) %>% select(c(s, starts_with('sample_qc')))
    setkey(dt, 's')
    return(dt)
}

recode_phenos_and_write <- function(dt, filepath) {
    # Rename to obtain the non-Europeans
    dt <- dt %>% mutate(
        white.british = ifelse(white.british == 1, "White-British", "Non white-British"),
        genetic.eur = ifelse(genetic.eur == 1, "European", "Non-European"),
        sequencing.batch = ifelse(sequencing.batch == 1, "Batch 1", "Batch 2")
        )

    dt <- dt %>% mutate(
        white.british= ifelse(is.na(white.british), "Non white-British", white.british),
        genetic.eur = ifelse(is.na(genetic.eur), "Non-European", genetic.eur)
        )

    fwrite(dt, file=filepath, sep='\t')
    system(paste("bgzip", filepath))
    return(dt)
}

dt_before <- read_and_key(SAMPLE_BEFORE_QC_FILE)
dt_after <- read_and_key(SAMPLE_AFTER_QC_FILE)

for (CHR in seq(2,22)) {
    # Input files
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

save_figures <- TRUE

dt_pheno <- fread(cmd=paste('zcat', PHENOFILE_CTS)) %>% 
    mutate(s=ID) %>% 
    select(s, ukbb.centre, sex, genotyping.array, sequencing.batch, white.british, genetic.eur)
setkey(dt_pheno, 's')
dt_before <- merge(dt_before, dt_pheno)
dt_after <- merge(dt_after, dt_pheno)

dt_before <- recode_phenos_and_write(dt_before, SAMPLE_BEFORE_COMBINED_QC_FILE)
dt_after <- recode_phenos_and_write(dt_after, SAMPLE_AFTER_COMBINED_QC_FILE)


PHENOFILE_CTS <- paste0('/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/',
    TRANCHE, '/UKBB_WES200k_filtered_cts_phenotypes.tsv.gz')


dt_before <- fread(SAMPLE_BEFORE_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE) %>% 
    mutate(phase='Before Variant QC')

dt_after <- fread(SAMPLE_AFTER_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE) %>%
    mutate(phase='After Variant QC')

dt <- bind_rows(dt_before, dt_after) %>%
    mutate(phase=factor(phase, levels=c('Before Variant QC', 'After Variant QC')))

## Split by Location.
# Colour by sequencing batch
# Number of singletons.

y_labels <- 'Location'
y_labels <- ''
alpha <- 0.8

# Number of singletons
create_pretty_boxplots(dt, aes(y=sample_qc.n_singleton, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    save_figure=save_figure, file=paste0(PLOTS, '09_nSingletons_by_centre'), y_label=y_labels,
    alpha=alpha)

# rHetHomVar
create_pretty_boxplots(dt, aes(y=sample_qc.r_het_hom_var, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    save_figure=save_figure, file=paste0(PLOTS, '09_rHetHomVar_by_centre'), y_label=y_labels,
    alpha=alpha)

# rInsertionDeletion
create_pretty_boxplots(dt, aes(y=sample_qc.r_insertion_deletion, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    save_figure=save_figure, file=paste0(PLOTS, '09_rInsertionDeletion_by_centre'), y_label=y_labels,
    alpha=alpha)

# rTiTv
create_pretty_boxplots(dt, aes(y=sample_qc.r_ti_tv, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    save_figure=save_figure, file=paste0(PLOTS, '09_rTiTv_by_centre'), n_ticks=5, y_label=y_labels,
    alpha=alpha)

q <- function(x) {
  print(subset(x, x < (mean(x) - 3*sd(x))| x > (mean(x) + 3*sd(x))))
  y <- subset(x, x < (mean(x) - 3*sd(x))| x > (mean(x) + 3*sd(x)))
  if(length(y)==0) y <- NA
  return(y)
}

# Number of singletons
p <- create_pretty_boxplots(dt, aes(y=sample_qc.n_singleton, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
ggsave(paste0(PLOTS, '09_nSingletons_by_centre_outliers_highlighted.pdf'), p, width=160, height=90, units='mm')

# rHetHomVar
p <- create_pretty_boxplots(dt, aes(y=sample_qc.r_het_hom_var, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
ggsave(paste0(PLOTS, '09_rHetHomVar_by_centre_outliers_highlighted.pdf'), p, width=160, height=90, units='mm')

# rInsertionDeletion
p <- create_pretty_boxplots(dt, aes(y=sample_qc.r_insertion_deletion, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
ggsave(paste0(PLOTS, '09_rInsertionDeletion_by_centre_outliers_highlighted.pdf'), p, width=160, height=90, units='mm')

# rTiTv
p <- create_pretty_boxplots(dt, aes(y=sample_qc.r_ti_tv, x=factor(ukbb.centre)),
    aes(color=factor(sequencing.batch)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
    legend=TRUE) + stat_summary(fun.y=q, geom="point", position=position_dodge(1), color='grey')
ggsave(paste0(PLOTS, '09_rTiTv_by_centre_outliers_highlighted.pdf'), p, width=160, height=90, units='mm')
