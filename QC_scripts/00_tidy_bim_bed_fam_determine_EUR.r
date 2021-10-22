library(bigsnpr)
library(data.table)
library(dplyr)
library(ggplot2)
library(randomForest)
library(ggrastr)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")
source("utils/helpers.r")

# Function to sanity check the random matrix theory proejctions
project_onto_ref_PCs <- function(bed.ref, obj.bed, n_PCs = 20, n_PC_plot = 8,
  outdir = "/well/lindgren/dpalmer/ukb_get_EUR/", filename = "output", strand_flip=FALSE) {

  test <- bed_projectPCA(bed.ref, obj.bed, k = n_PCs, strand_flip = FALSE)#, ncores = nb_cores())

  plot(test$obj.svd.ref)
  ggsave(paste0(outdir, filename, "_scree.pdf"), width = 13, height = 7)
  plot(test$obj.svd, type = "loadings", loadings = 1:n_PCs, coeff = 0.4)
  ggsave(paste0(outdir, filename, "_loadings.pdf"), width = 13, height = 7)

  PC.ref <- predict(test$obj.svd.ref)
  proj1 <- test$simple_proj
  proj2 <- test$OADP_proj
  save("test", "PC.ref", "proj1", "proj2", file=paste0(outdir, filename, "_data.Rdata"))

  # shrinkage coefficients
  shrinkage <- unname(sapply(1:n_PCs, function(k) {
    MASS::rlm(proj2[, k] ~ proj1[, k] + 0, maxit=100)$coef
  }))
  print(round(shrinkage, 2))

  ind <- seq(1, min(20e3, nrow(proj1)))

  plot_grid(plotlist = lapply(1:floor(n_PC_plot/2), function(k) {
    k1 <- 2 * k - 1
    k2 <- 2 * k
    plot_grid(
      qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
        geom_point(aes(proj1[ind, k1], proj1[ind, k2]), color = "red") +
        theme_bigstatsr(0.5) +
        labs(x = paste0("PC", k1), y = paste0("PC", k2)),
      qplot(PC.ref[, k1], PC.ref[, k2], size = I(2)) +
        geom_point(aes(proj2[ind, k1], proj2[ind, k2]), color = "blue") +
        theme_bigstatsr(0.5) +
        labs(x = paste0("PC", k1), y = paste0("PC", k2)),
      scale = 0.95
    )
  }), nrow = floor(sqrt(n_PCs)))
  ggsave(file=paste0(outdir, filename, "_shrinkage.pdf"), width=15, height=15)
  return(round(shrinkage, 2))

}

# Generate the fam file for all people, so that that it can be read by plink.
ukb_fam_lindgren <- "/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"
ukb_fam <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam"
fwrite(fread(ukb_fam) %>% mutate(V6=-9), file=ukb_fam, sep=' ', col.names=FALSE)

# First, generate cleaned up genotype data
write(sapply(1:22, function(chr) {
  c(paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr", chr, "_v2.bed"),
    paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr", chr, "_v2.bim"),
    ukb_fam)
}), '/well/lindgren/dpalmer/ukb_get_EUR/plink_ukb_filepaths.txt', ncolumns = 3)

# If this dies, it's likely due to memory. Requires at least 10 slots on qe.
snp_plinkQC(
  plink.path = '/well/lindgren/dpalmer/plink',
  prefix.in = '/well/lindgren/dpalmer/ukb_get_EUR/plink_ukb_filepaths.txt',
  file.type = "--merge-list",
  prefix.out = "/well/lindgren/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined",
  geno = 0.01,
  autosome.only = TRUE
)

obj.bed <- bed('/well/lindgren/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined.bed')
bed.ref <- bed(download_1000G("/well/lindgren/dpalmer/ukb_get_EUR/data"))

outdir <- '/well/lindgren/dpalmer/ukb_get_EUR/'
filename <- 'ukb_projected_to_1kg_PCs'
n_PCs <- 10
save_figures <- TRUE
perform_plotting <- TRUE
creating_new_EUR_def <- TRUE

project_onto_ref_PCs(bed.ref, obj.bed, filename=filename, outdir=outdir, n_PCs=n_PCs)

# Grab 1000G information containing populations
# download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped",
#               destfile = "/well/lindgren/UKBIOBANK/dpalmer/1kg_for_EUR_assign/20130606_g1k.ped")
dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/1kg_for_EUR_assign/20130606_g1k.ped") %>% 
  mutate(sample.ID = `Individual ID`, population=Population) %>% 
  select(sample.ID, population) %>%
  # Determine the superpopulations
  mutate(super.population = case_when(
    population %in% c("CHB", "JPT", "CHS", "CDX", "KHV") ~ "EAS",
    population %in% c("CEU", "TSI", "FIN", "GBR", "IBS") ~ "EUR",
    population %in% c("YRI", "LWK", "MAG", "MSL", "ESN", "ASW", "ACB", "GWD") ~ "AFR",
    population %in% c("AMR", "PUR", "CLM", "PEL", "MXL") ~ "AMR",
    population %in% c("GIH", "PJL","BEB", "STU", "ITU") ~ "SAS",
    TRUE ~ "other"
  )
)

load(paste0(outdir, filename, "_data.Rdata"))
PC.ref <- data.table(PC.ref)
names(PC.ref) <- paste0("PC", seq(1,n_PCs))

dt <- cbind(data.table(bed.ref$fam), PC.ref)
setkey(dt, "sample.ID")
setkey(dt_pop, "sample.ID")
dt_1kg <- merge(dt, dt_pop)

proj2 <- data.table(proj2)
names(proj2) <- paste0("PC", seq(1,n_PCs))
dt_ukb <- cbind(data.table(obj.bed$fam), proj2)

# First determine the Europeans, then we can remove Finns by training on the Europeans.
# Need to remove the Finns.
# Merge in and obtain data.table for the PC locations, IDs, populations and superpopulations (dt_train).
# Use this to then assign individuals within the test data.

# Now, need to grab the population labels to train the random forest.
dt_train = dt_1kg %>%
  select(c(super.population, population, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

dt_predict = dt_ukb %>%
  select(c(sample.ID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

PCs_to_use <- paste0('PC', seq(1,4))

# Determine a classifier.
set.seed(160487)
rf <- randomForest(x=dt_train %>% select(PCs_to_use), y=as.factor(as.character(dt_train$super.population)), ntree=10000)
rf_probs <- predict(rf, dt_predict %>% select(PCs_to_use), type='prob')

# Take a look at the histogram of European assignment
# PDF of EUR cluster assignment
create_pretty_hist(data.table(rf_probs), aes(x=EUR), 'European cluster assignment certainty', T_European_RF,
    xlim=c(0.95, 1),
    title='Assignment of Europeans', save_figure=save_figures, file=paste0(PLOTS, '00_EUR_cluster_hist'))
# CDFs, no splitting.
create_pretty_cumulative(data.table(rf_probs), aes(x=EUR), 'European cluster assignment certainty', T_European_RF,
    xlim=c(0.95, 1), title='Assignment of Europeans', save_figure=save_figures, file=paste0(PLOTS, '00_EUR_cluster_cdf'))

check_thres <- function(row, threshold) {
  return(!any(row > threshold))
}

unsure <- apply(rf_probs, 1, check_thres, T_European_RF)
classification <- as.character(predict(rf, dt_predict %>% select(PCs_to_use)))
dt_predict$classification_loose <- as.factor(classification)
classification[unsure] <- 'unsure'
dt_predict$classification_strict <- as.factor(classification)

dt_predict <- dt_predict %>% mutate(sample.ID = as.character(sample.ID)) %>% select(sample.ID, classification_strict, classification_loose, starts_with("PC"))
dt_predict <- data.table(dt_predict)
setkeyv(dt_predict, c("sample.ID", paste0("PC", seq(1,10))))
setkeyv(dt_1kg, c("sample.ID", paste0("PC", seq(1,10))))
dt_classify <- merge(dt_predict, dt_1kg, all=TRUE)

PCs <- c(1,3,5)

if (perform_plotting)
{
    dt_plot <- dt_classify %>% filter(!is.na(classification_loose))
    # levels(dt_plot$classification_loose) <- c(levels(dt_plot$classification_loose) ,"1000 genomes")
    # dt_plot$classification_loose[is.na(dt_plot$classification_loose)] <- "1000 genomes"

    for (i in PCs)
    {
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_loose')
        p <- create_pretty_scatter(dt_plot, aes,
          save_figure=save_figures, file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_EUR_loose'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
        if (save_figures) {
            p <- p + geom_point(data=dt_1kg %>% filter(super.population == "EUR"),
              mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1)),
              inherit.aes=FALSE, shape=4, show.legend=FALSE)
            ggsave(file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_EUR_loose_1kg_labelled.pdf'), width=160, height=90, units='mm')
            ggsave(file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_EUR_loose_1kg_labelled.jpg'), width=160, height=90, units='mm', dpi=500)
        }

        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_strict')
        p <- create_pretty_scatter(dt_plot, aes,
          save_figure=save_figures, file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_EUR_strict'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
    }
}


if (creating_new_EUR_def) { 
    dt_out_classify_strict <- dt_classify %>% filter(classification_strict == 'EUR') %>% select(sample.ID)
    dt_out_classify_loose <- dt_classify %>% filter(classification_loose == 'EUR') %>% select(sample.ID)

    # Print the number of remaining samples.
    print(paste0("Number of European samples using strict classifier: ", dim(dt_out_classify_strict)[1]))
    print(paste0("Number of European samples using loose classifier: ", dim(dt_out_classify_loose)[1]))
}

# Next, need to run the procedure again, labelling and removing Finns.
# First, filter to the collection of retained Europeans
dt_ukb_EUR <- dt_predict %>% filter(classification_strict == "EUR") %>% transmute(V1 = sample.ID)
setkey(dt_ukb_EUR, "V1")
new_fam <- fread("/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam", colClasses=c(rep("character", 2), rep("integer", 4)), key='V1')
new_fam <- merge(dt_ukb_EUR, new_fam)
fwrite(new_fam %>% select(V1, V2), sep="\t", row.names=FALSE, col.names=FALSE, file="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink_EUR.tsv")

# Filter the 1000G to the Europeans and use plink to restrict.

dt_1kg_EUR <- dt_1kg %>% filter(super.population == "EUR") %>% transmute(V2 = sample.ID)
setkey(dt_1kg_EUR, "V2")
new_fam_1kg <- fread("/well/lindgren/dpalmer/ukb_get_EUR/data/1000G_phase3_common_norel.fam", colClasses=c("integer", "character", rep("integer", 4)), key='V2')
new_fam_1kg <- merge(dt_1kg_EUR, new_fam_1kg)
setcolorder(new_fam_1kg, paste0("V", seq(1,6)))
fwrite(new_fam_1kg %>% select(V1, V2), sep="\t", row.names=FALSE, col.names=FALSE, file="/well/lindgren/dpalmer/ukb_get_EUR/data/1000G_phase3_common_norel_EUR.tsv")

# Calls to plink

# 1000G Europeans
system(paste("/well/lindgren/dpalmer/plink",
  "--bfile /well/lindgren/dpalmer/ukb_get_EUR/data/1000G_phase3_common_norel",
  "--keep /well/lindgren/dpalmer/ukb_get_EUR/data/1000G_phase3_common_norel_EUR.tsv",
  "--make-bed --out /well/lindgren/dpalmer/ukb_get_EUR/data/1000G_phase3_common_norel_EUR_output"))

# UKBB Europeans
system(paste("/well/lindgren/dpalmer/plink",
  "--bfile /well/lindgren/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined",
  "--keep /well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink_EUR.tsv",
  "--make-bed --out /well/lindgren/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined_EUR"))

obj.bed <- bed('/well/lindgren/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined_EUR.bed')
bed.ref <- bed('/well/lindgren/dpalmer/ukb_get_EUR/data/1000G_phase3_common_norel_EUR_output.bed')

filename <- 'ukb_projected_to_1kg_EUR_PCs'
project_onto_ref_PCs(bed.ref, obj.bed, filename=filename, outdir=outdir, n_PCs=n_PCs)

load(paste0(outdir, filename, "_data.Rdata"))
PC.ref <- data.table(PC.ref)
names(PC.ref) <- paste0("PC", seq(1,n_PCs))

dt <- cbind(data.table(bed.ref$fam), PC.ref)
setkey(dt, "sample.ID")
setkey(dt_pop, "sample.ID")
dt_1kg <- merge(dt, dt_pop)

proj2 <- data.table(proj2)
names(proj2) <- paste0("PC", seq(1,n_PCs))
dt_ukb <- cbind(data.table(obj.bed$fam), proj2)

# Determine the Finns
# Merge in and obtain data.table for the PC locations, IDs, populations and superpopulations (dt_train).
# Use this to then assign individuals within the test data.

# Now, need to grab the population labels to train the random forest.
dt_train = dt_1kg %>%
  select(c(super.population, population, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)) %>% mutate(FIN = ifelse(dt_1kg$population == "FIN", "FIN", "non-FIN"))

dt_predict = dt_ukb %>%
  select(c(sample.ID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

PCs_to_use <- paste0('PC', seq(1))

# Determine a classifier.
set.seed(160487)
rf <- randomForest(x=dt_train %>% select(PCs_to_use), y=as.factor(as.character(dt_train$FIN)), ntree=10000)
rf_probs <- predict(rf, dt_predict %>% select(PCs_to_use), type='prob')

# Take a look at the histogram of Finnish assignment
# PDF of EUR cluster assignment
create_pretty_hist(data.table(rf_probs), aes(x=`non-FIN`), 'Non-Finnish cluster assignment certainty', T_Finnish_RF,
    xlim=c(0.95, 1),
    title='Assignment of Non-Finns', save_figure=save_figures, file=paste0(PLOTS, '00_FIN_cluster_hist'))
# CDFs, no splitting.
create_pretty_cumulative(data.table(rf_probs), aes(x=`non-FIN`), 'Non-Finnish cluster assignment certainty', T_Finnish_RF,
    xlim=c(0.95, 1), title='Assignment of Non-Finns', save_figure=save_figures, file=paste0(PLOTS, '00_FIN_cluster_cdf'))

unsure <- apply(rf_probs, 1, check_thres, T_Finnish_RF)
classification <- as.character(predict(rf, dt_predict %>% select(PCs_to_use)))
dt_predict$classification_loose <- as.factor(classification)
classification[unsure] <- 'unsure'
dt_predict$classification_strict <- as.factor(classification)

dt_predict <- dt_predict %>% mutate(sample.ID = as.character(sample.ID)) %>% select(sample.ID, classification_strict, classification_loose, starts_with("PC"))
dt_predict <- data.table(dt_predict)
setkeyv(dt_predict, c("sample.ID", paste0("PC", seq(1,10))))
setkeyv(dt_1kg, c("sample.ID", paste0("PC", seq(1,10))))
dt_classify <- merge(dt_predict, dt_1kg, all=TRUE)

PCs <- c(1)

if (perform_plotting)
{
    for (i in PCs)
    {
        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_loose')
        p <- create_pretty_scatter(dt_classify, aes,
          save_figure=save_figures, file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_FIN_loose'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
        if (save_figures) {
            p <- p + geom_point(data=dt_1kg %>% filter(population == "FIN"),
              mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1)),
              inherit.aes=FALSE, shape=4, show.legend=FALSE)
            ggsave(file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_FIN_loose_1kg_labelled.pdf'), width=160, height=90, units='mm')
            ggsave(file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_FIN_loose_1kg_labelled.jpg'), width=160, height=90, units='mm', dpi=500)
        }

        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_strict')
        p <- create_pretty_scatter(dt_classify, aes,
          save_figure=save_figures, file=paste0(PLOTS,'00_PC',i,'_PC',i+1,'_classify_FIN_strict'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
    }
}


if (creating_new_EUR_def) { 
    dt_out_classify_strict <- dt_classify %>% filter(classification_strict == 'non-FIN') %>% select(sample.ID)
    dt_out_classify_loose <- dt_classify %>% filter(classification_loose == 'non-FIN') %>% select(sample.ID)

    # Print the number of remaining samples.
    print(paste0("Number of Non-Finnish European samples using strict classifier: ", dim(dt_out_classify_strict)[1]))
    print(paste0("Number of Non-Finnish European samples using loose classifier: ", dim(dt_out_classify_loose)[1]))
}

# Write the result to disk.
fwrite(dt_out_classify_strict, file = "/well/lindgren/dpalmer/ukb_get_EUR/data/final_EUR_list.tsv", col.names=FALSE)
