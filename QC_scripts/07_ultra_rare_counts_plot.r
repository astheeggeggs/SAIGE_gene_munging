library(hexbin)
library(gridExtra)
library(ggplot2)
library(ggExtra)
library(data.table)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--tranche", default='200k', help = "Which exome sequencing tranche?")
args <- parser$parse_args()

TRANCHE <- args$tranche
CHR <- 1

# Input files
URV_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
	'/ukb_wes_qc/data/samples/08_URVs_chr',  CHR, '.tsv.bgz')

# Output files
URV_COMBINED_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
	'/ukb_wes_qc/data/samples/08_URVs.tsv')

dt <- fread(cmd = paste('zcat', URV_FILE),
    stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(c(s, starts_with('n_')))
setkey(dt, 's')

for (CHR in seq(2,22)) {
    # Input files
    URV_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
	'/ukb_wes_qc/data/samples/07_URVs_chr',  CHR, '.tsv.bgz')

    dt_tmp <- fread(cmd = paste('zcat', URV_FILE),
        stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(c(s, starts_with('n_')))
    names(dt_tmp) <- c('s', paste0('tmp_',names(dt_tmp)[-1]))
    setkey(dt_tmp, 's')

    dt  <- merge(dt, dt_tmp)

    dt <- dt %>% mutate(
        # Counts
        n_coding_URV_SNP = n_coding_URV_SNP + tmp_n_coding_URV_SNP,
        n_coding_URV_indel = n_coding_URV_indel + tmp_n_coding_URV_indel,
        n_URV_pLoF = n_URV_pLoF + tmp_n_URV_pLoF,
        n_URV_missense = n_URV_missense + tmp_n_URV_missense,
        n_URV_synonymous = n_URV_synonymous + tmp_n_URV_synonymous
    ) %>% select(c(s, starts_with('n_')))
    setkey(dt, 's')
}

save_figures <- TRUE
dt_pheno <- create_pheno_dt(TRANCHE)
dt <- merge(dt, dt_pheno)

fwrite(dt, file=URV_COMBINED_FILE, sep='\t')
system(paste("bgzip", URV_COMBINED_FILE))

dt <- dt[sample(nrow(dt), replace=FALSE),]

# Scatters of URV-SNPs against URV-Indels.
p <- ggplot(dt, aes(x=n_coding_URV_SNP, y=n_coding_URV_indel, colour=factor(genetic.eur.no.fin.oct2021))) + 
geom_point(size=0.5, alpha=0.5) + 
scale_color_d3('category10') + theme_minimal() + labs(x='n URV SNPs', y='n URV indels', color='')
p <- ggExtra::ggMarginal(p, type = "density",
  xparams = list(adjust=1), yparams=list(adjust=1.5))
ggsave(paste0(PLOTS, '07_URVs_SNPs_vs_indels_by_NFE.pdf'), p, width=160, height=90, units='mm')

p <- ggplot(dt, aes(x=n_coding_URV_SNP, y=n_coding_URV_indel, colour=factor(white.british))) + 
geom_point(size=0.5, alpha=0.5) + 
scale_color_d3('category10') + theme_minimal() + labs(x='n URV SNPs', y='n URV indels', color='')
p <- ggExtra::ggMarginal(p, type = "density",
  xparams = list(adjust=1), yparams=list(adjust=1.5))
ggsave(paste0(PLOTS, '07_URVs_SNPs_vs_indels_by_british.pdf'), p, width=160, height=90, units='mm')

y_labels <- c('UKBB centre', 'Sequencing batch')
y_label_batch <- c('', '')
titles <- c(
    'Number of singletons, split by UKBB Centre',
    'Number of singletons, split by self reported ethnicity',
	'Number of singletons, split by NFE status'
    )

create_pretty_boxplots(dt, aes(x=factor(ukbb.centre), y=(n_coding_URV_SNP + n_coding_URV_indel)), aes(colour=factor(genetic.eur.no.fin.oct2021)), 
	y_label=y_label_batch[1], x_label='Number of Singletons', key_label='',
	title=titles[1], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'07_URVs_by_centre'), xlim=c(0,10),
	threshold=T_nURVSNP)
create_pretty_boxplots(dt, aes(x=factor(self.report.ethnicity), y=(n_coding_URV_SNP + n_coding_URV_indel)), aes(colour=factor(genetic.eur.no.fin.oct2021)), 
    y_label=y_label_batch[1], x_label='Number of Singletons', key_label='',
    title=titles[2], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'07_URVs_by_ethnicity'), xlim=c(0,10),
    threshold=T_nURVSNP)
create_pretty_boxplots(dt, aes(x=factor(genetic.eur.no.fin.oct2021), y=(n_coding_URV_SNP + n_coding_URV_indel)), aes(colour=factor(ukbb.centre)),
	y_label=y_label_batch[2], x_label='Number of Singletons', key_label='UKBB centre',
	title=titles[3], legend=TRUE, save_figure=save_figures,  file=paste0(PLOTS,'07_URVs_by_batch'),  xlim=c(0,10),
	threshold=T_nURVSNP)
