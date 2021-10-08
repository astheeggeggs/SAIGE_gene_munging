library(data.table)
library(dplyr)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")
source("utils/helpers.r")

# Need to combine and plot - running sample QC metrics separately across chromosomes,
# so need to amalgamate afterwards.

CHR <- 1

# Input files
INITIAL_SAMPLE_QC_FILE <- paste0(
    '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/03_chr', CHR, '_initial_sample_qc.tsv.bgz')

# Output files
INITIAL_COMBINED_SAMPLE_QC_FILE <- paste0(
    '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/03_initial_sample_qc.tsv')

dt <- fread(cmd = paste('zcat', INITIAL_SAMPLE_QC_FILE),
    stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(c(s, starts_with('sample_qc'), starts_with('gq'), starts_with('dp')))
dt <- dt %>% mutate(
    dp.stdev_sum = (dp.stdev)^2 * dp.n,
    gq.stdev_sum = (gq.stdev)^2 * gq.n
    )

setkey(dt, 's')

for (CHR in seq(2,22)) {
    # Input files
    cat(paste0("chromosome ", CHR, "\n"))
    INITIAL_SAMPLE_QC_FILE <- paste0(
        '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
        '/ukb_wes_qc/data/samples/03_chr', CHR, '_initial_sample_qc.tsv.bgz')

    dt_tmp <- fread(
        cmd = paste('zcat', INITIAL_SAMPLE_QC_FILE),
        stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(c(s, starts_with('sample_qc'), starts_with('gq'), starts_with('dp')))
    names(dt_tmp) <- c('s', paste0('tmp_',names(dt_tmp)[-1]))
    setkey(dt_tmp, 's')

    dt  <- merge(dt, dt_tmp)

    dt <- dt %>% mutate(
        # Mins
        sample_qc.dp_stats.min = pmin(sample_qc.dp_stats.min, tmp_sample_qc.dp_stats.min),
        sample_qc.gq_stats.min = pmin(sample_qc.gq_stats.min, tmp_sample_qc.gq_stats.min),
        
        # Maxs
        sample_qc.dp_stats.max = pmax(sample_qc.dp_stats.max, tmp_sample_qc.dp_stats.max),
        sample_qc.gq_stats.max = pmax(sample_qc.gq_stats.max, tmp_sample_qc.gq_stats.max),
        
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
        sample_qc.n_star = sample_qc.n_star + tmp_sample_qc.n_star,

        # Means
        dp.sum = dp.sum + tmp_dp.sum,
        gq.sum = gq.sum + tmp_gq.sum,

        # Stdevs
        dp.stdev_sum = dp.stdev_sum + (tmp_dp.stdev)^2 * tmp_dp.n,
        gq.stdev_sum = dp.stdev_sum + (tmp_gq.stdev)^2 * tmp_gq.n,

        # N
        dp.n = dp.n + tmp_dp.n,
        gq.n = gq.n + tmp_gq.n
    ) %>% select(c(s, starts_with('sample_qc'), starts_with('dp'), starts_with('gq')))
    setkey(dt, 's')
}

dt <- dt %>% mutate(
    # Ratios
    sample_qc.call_rate = sample_qc.n_called / (sample_qc.n_called + sample_qc.n_not_called + sample_qc.n_filtered),
    sample_qc.r_ti_tv = ifelse(sample_qc.n_transversion == 0, NA, sample_qc.n_transition / sample_qc.n_transversion),
    sample_qc.r_het_hom_var = ifelse(sample_qc.n_hom_var == 0, NA, sample_qc.n_het / sample_qc.n_hom_var),
    sample_qc.r_insertion_deletion = ifelse(sample_qc.n_deletion == 0, NA, sample_qc.n_insertion / sample_qc.n_deletion),
    
    # Means
    sample_qc.dp_stats.mean = dp.sum / dp.n,
    sample_qc.gq_stats.mean = gq.sum / gq.n,

    # Stdevs
    sample_qc.dp_stats.stdev = sqrt(dp.stdev_sum / dp.n),
    sample_qc.gq_stats.stdev = sqrt(gq.stdev_sum / gq.n)
    ) %>% select(c('s', starts_with('sample_qc')))

setkey(dt, 's')

# Create plots
save_figures <- TRUE
dt_pheno <- create_pheno_dt(TRANCHE)
dt <- merge(dt, dt_pheno)

fwrite(dt, file=INITIAL_COMBINED_SAMPLE_QC_FILE, sep='\t')
system(paste("bgzip", INITIAL_COMBINED_SAMPLE_QC_FILE))

# Want to combine all the sample metrics (because we're limited to parallelisation within chr).

# PDFs, no splitting.
create_pretty_hist(dt, aes(x=sample_qc.call_rate), 'Call Rate', T_sample_callRate,
    title='Call Rate', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_03_callRate_hist'))
create_pretty_hist(dt, aes(x=sample_qc.dp_stats.mean), 'Mean Depth', T_dpMean,
    binwidth=0.25, xlim=c(NA, NA), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_03_dpMean_hist'))
create_pretty_hist(dt, aes(x=sample_qc.gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    binwidth=0.01, xlim=c(NA, NA), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_03_gqMean_hist'))

# CDFs, no splitting.
create_pretty_cumulative(dt, aes(sample_qc.call_rate), 'Call Rate', T_sample_callRate,
    xlim=c(NA, NA), title='Call Rate', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_03_callRate_cdf'))
create_pretty_cumulative(dt, aes(sample_qc.dp_stats.mean), 'Mean Depth', T_dpMean,
    xlim=c(NA,NA), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_03_dpMean_cdf'))
create_pretty_cumulative(dt, aes(sample_qc.gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    xlim=c(NA,NA), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_03_gqMean_cdf'))

legend_batch <- TRUE
save_figures <- TRUE
y_label_batch <- 'UKBB centre'
titles <- c('Call Rate',
    'Mean Depth (DP)',
    'Mean Genotype Quality (GQ)')
titles <- c('', '', '')
alpha <- 0.8
jitter_size <- 0.1

# Split by UKB Centre
create_pretty_boxplots(dt, aes(x=factor(ukbb_centre), y=sample_qc.call_rate), aes(color=factor(sequencing_batch)),
    T_sample_callRate, x_label='Call Rate', y_label=y_label_batch, key_label='Sequencing batch',
    xlim=quantile(dt$call_rate, c(0.01, 0.99)), legend=legend_batch, title=titles[1], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_callRate_by_centre'), n_ticks=5, alpha=alpha, jitter_size=jitter_size)
create_pretty_boxplots(dt, aes(x=factor(ukbb_centre), y=sample_qc.dp_stats.mean), aes(color=factor(sequencing_batch)),
    T_dpMean, x_label='Mean Depth', y_label=y_label_batch, key_label='Sequencing batch',
    xlim=quantile(dt$dp_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[2], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_dpMean_by_centre'), alpha=alpha, jitter_size=jitter_size)
create_pretty_boxplots(dt, aes(x=factor(ukbb_centre), y=sample_qc.gq_stats.mean), aes(color=factor(sequencing_batch)),
    T_gqMean, x_label='Mean Genotype Quality', y_label=y_label_batch, key_label='Sequencing batch',
    xlim=quantile(dt$gq_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[3], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_gqMean_by_centre'), alpha=alpha, jitter_size=jitter_size)

y_label_batch <- ''

# Split by LOCATION - NFE vs non-NFE
create_pretty_boxplots(dt, aes(x=factor(genetic_eur_no_fin_oct2021), y=sample_qc.call_rate), aes(color=factor(sequencing_batch)),
    T_sample_callRate, x_label='Call Rate', y_label=y_label_batch, key_label='Sequencing batch',
    xlim=quantile(dt$call_rate, c(0.01, 0.99)), legend=legend_batch, title=titles[1], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_callRate_by_NFE'), n_ticks=5, alpha=alpha, jitter_size=jitter_size)
create_pretty_boxplots(dt, aes(x=factor(genetic_eur_no_fin_oct2021), y=sample_qc.dp_stats.mean), aes(color=factor(sequencing_batch)),
    T_dpMean, x_label='Mean Depth', y_label=y_label_batch, key_label='Sequencing batch',
    xlim=quantile(dt$dp_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[2], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_dpMean_by_NFE'), alpha=alpha, jitter_size=jitter_size)
create_pretty_boxplots(dt, aes(x=factor(genetic_eur_no_fin_oct2021), y=sample_qc.gq_stats.mean), aes(color=factor(sequencing_batch)),
    T_gqMean, x_label='Mean Genotype Quality', y_label=y_label_batch, key_label='Sequencing batch',
    xlim=quantile(dt$gq_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[3], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_gqMean_by_NFE'), alpha=alpha, jitter_size=jitter_size)

# Split by self reported ancestry
create_pretty_boxplots(dt, aes(x=factor(self_report_ethnicity), y=sample_qc.call_rate), aes(color=factor(genetic_eur_no_fin_oct2021)),
    T_sample_callRate, x_label='Call Rate', y_label=y_label_batch, key_label='',
    xlim=quantile(dt$call_rate, c(0.01, 0.99)), legend=legend_batch, title=titles[1], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_callRate_by_anc'), n_ticks=5, alpha=alpha, jitter_size=jitter_size)
create_pretty_boxplots(dt, aes(x=factor(self_report_ethnicity), y=sample_qc.dp_stats.mean), aes(color=factor(genetic_eur_no_fin_oct2021)),
    T_dpMean, x_label='Mean Depth', y_label=y_label_batch, key_label='',
    xlim=quantile(dt$dp_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[2], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_dpMean_by_anc'), alpha=alpha, jitter_size=jitter_size)
create_pretty_boxplots(dt, aes(x=factor(self_report_ethnicity), y=sample_qc.gq_stats.mean), aes(color=factor(genetic_eur_no_fin_oct2021)),
    T_gqMean, x_label='Mean Genotype Quality', y_label=y_label_batch, key_label='',
    xlim=quantile(dt$gq_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[3], save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_03_gqMean_by_anc'), alpha=alpha, jitter_size=jitter_size)

