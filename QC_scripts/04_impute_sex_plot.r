rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
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

IMPUTESEX_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/04_imputesex.tsv.bgz')
Y_NCALLED_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/04_ycalled.tsv.bgz')
SEXCHECK_LIST <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/04_sexcheck.remove.sample_list')

dt <- fread(cmd = paste('zcat', IMPUTESEX_FILE), sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(imputed_sex=as.factor(ifelse(impute_sex.is_female == TRUE, 'Female', 'Male')))
dt_y <- fread(cmd = paste('zcat', Y_NCALLED_FILE), sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)
dt <- merge(dt, dt_y, by='s')

colors <- pal_d3('category20')(20)[c(1,2)]
fills <- pal_d3('category20')(20)[c(11,12)]

create_pretty_cumulative(dt, aes(impute_sex.f_stat), 'F-statistic', T_impute_sex,
    xlim=c(NA, NA), title='Cumulative Distribution of F-statistic', save_figure=TRUE, file=paste0(PLOTS,'04_F_stat_cdf'))

p <- ggplot(dt, aes(x=impute_sex.f_stat, fill=imputed_sex)) +
  geom_histogram(binwidth=0.01, alpha=0.8, color='#7f7f7f') +
  scale_fill_manual(values=fills, limits=c('Male', 'Female')) +
  labs(x='X chromosome F-statistic',
       y='Count',
       title='',
       fill='Imputed Sex') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
  theme_minimal() +
  theme(axis.title.x=element_text(margin=ggplot2::margin(t=10)),
        axis.title.y=element_text(margin=ggplot2::margin(r=10)),
        plot.title=element_text(hjust=0.5)) +
  geom_vline(xintercept=T_impute_sex, linetype='dashed')

ggsave(paste0(PLOTS, '04_imputesex_histogram', '.pdf'), p, width=160, height=90, units='mm')

dt <- dt %>% mutate(phenotype.Submitted_Gender = ifelse(phenotype.Submitted_Gender == "F", "Female", "Male"))
if(any(is.na(dt$phenotype.Submitted_Gender))) {
  dt <- dt %>% mutate(phenotype.Submitted_Gender = ifelse(is.na(phenotype.Submitted_Gender), "Unknown", phenotype.Submitted_Gender))
}

if(any(is.na(dt$phenotype.ukbb_centre))) {
  dt <- dt %>% mutate(phenotype.ukbb_centre = ifelse(is.na(phenotype.ukbb_centre), "Unknown", phenotype.ukbb_centre))
}

p <- ggplot(dt, aes(x=impute_sex.f_stat, y=factor(phenotype.ukbb_centre), colour=phenotype.Submitted_Gender)) +
  geom_jitter_rast(width=0, height=0.2, size=1, alpha=0.2, stroke=0.05) + 
  theme_minimal() +
  geom_vline(xintercept=T_impute_sex, linetype='dashed') +
  labs(x='X chromosome F-statistic',
       y='Location',
       color='Reported Sex') 

ggsave(paste0(PLOTS, '04_imputesex_scatter_box', '.pdf'), p, width=160, height=90, units='mm')

dt_false <- dt %>% filter(
  (impute_sex.f_stat > T_impute_sex & phenotype.Submitted_Gender == 'Female') |
  (impute_sex.f_stat < T_impute_sex & phenotype.Submitted_Gender == 'Male') |
  (impute_sex.f_stat > T_impute_sex & n_called < 100))

# Plots of sex estimates using Giulios plotting method.
p <- ggplot(dt, aes(x=impute_sex.f_stat, y=n_called, colour=factor(phenotype.Submitted_Gender))) +
geom_point_rast(size=0.5) + 
labs(x='X chromosome F-statistic', y='Number of calls in Y', color='Reported Sex') +
scale_color_d3('category10') +
scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
geom_point_rast(data=dt_false, aes(x=impute_sex.f_stat, y=n_called), size=0.5) + 
theme_minimal()
ggsave(paste0(PLOTS, '04_imputesex_scatter', '.pdf'), p, width=160, height=90, units='mm')

dt_out <- dt_false %>% select(s)
write.table(dt_out, file=SEXCHECK_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)
