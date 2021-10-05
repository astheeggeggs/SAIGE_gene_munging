library(bigsnpr)
library(data.table)
library(dplyr)
library(ggplot2)
library(randomForest)

# Function to sanity check the random matrix theory proejctions
project_onto_ref_PCs <- function(bed.ref, obj.bed, n_PCs = 20, n_PC_plot = 8,
  outdir = "/well/lindgren/dpalmer/ukb_get_EUR/", filename = "output", strand_flip=FALSE) {

  test <- bed_projectPCA(bed.ref, obj.bed, k = n_PCs, ncores = nb_cores(), strand_flip = FALSE)

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

# First, generate cleaned up genotype data
write(sapply(1:22, function(chr) {
  c(paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr", chr, "_v2.bed"),
    paste0("/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr", chr, "_v2.bim"),
    "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb11867_cal_chr1_v2_s488363_for_plink.fam")
}), '/well/lindgren/dpalmer/ukb_get_EUR/plink_ukb_filepaths.txt', ncolumns = 3)

snp_plinkQC(
  plink.path = '/well/lindgren/dpalmer/plink',
  prefix.in = '/well/lindgren/dpalmer/ukb_get_EUR/plink_ukb_filepaths.txt',
  file.type = "--merge-list",
  prefix.out = "/well/lindgren/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined",
  geno = 0.01,
  autosome.only = TRUE,
)

obj.bed <- bed('/well/lindgren/dpalmer/ukb_get_EUR/data/ukb_autosomes_combined.bed')
bed.ref <- bed(download_1000G("/well/lindgren/dpalmer/ukb_get_EUR/data"))

project_onto_ref_PCs(bed.ref, obj.bed, filename='ukb_projected_to_1kg_PCs')

dt <- data.table(bed.ref$fam)
setkey(dt, "sample.ID")

dt_pop <- fread("/well/lindgren/UKBIOBANK/dpalmer/1kg_for_EUR_assign/20130606_g1k.ped")
# Need to remove the Finns.
# Merge in and obtain data.table for the PC locations, IDs, populations and superpopulations (dt_train).
# Use this to then assign individuals within the test data.

# Now, need to grab the population labels to train the random forest.
dt_train = filter(dt_1kg, PHENOTYPE_COARSE=='1KG') %>%
  select(c(SUPER_POPULATION, POPULATION, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

dt_predict = filter(dt_1kg, PHENOTYPE_COARSE!='1KG') %>%
  select(c(s, POPULATION, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10))

PCs_to_use <- paste0('PC', seq(1,4))

# Determine a classifier.
set.seed(160487)
rf <- randomForest(x=dt_train[PCs_to_use], y=as.factor(as.character(dt_train$SUPER_POPULATION)), ntree=10000)
rf_probs <- predict(rf, dt_predict[PCs_to_use], type='prob')

check_thres <- function(row, threshold) {
  return(!any(row > threshold))
}

unsure <- apply(rf_probs, 1, check_thres, T_European_RF)
classification <- as.character(predict(rf, df_predict[PCs_to_use]))
df_predict$classification_loose <- as.factor(classification)
classification[unsure] <- 'unsure'
df_predict$classification_strict <- as.factor(classification)

df_classify <- df_predict %>% select(s, classification_strict, classification_loose) %>% inner_join(df_1kg, by='s')

# Include option to use an existing classification

if (!creating_new_EUR_def) {
    # Read in the previously strictly defined Europeans
    european_samples <- fread(EUROPEAN_SAMPLES_STRICT, sep='\t', stringsAsFactors=FALSE, header=FALSE, data.table=FALSE)
    df_classify$classification_strict <- ifelse(df_classify$s %in% european_samples$V1, 'European', 'Non-European')
}

if (perform_plotting) {
    for (i in PCs) {

        if (creating_new_EUR_def)
        {
            aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_loose')
            p <- create_pretty_scatter(df_classify, aes,
              save_figure=save_figures, file=paste0(PLOTS,'10_PC',i,'_PC',i+1,'_classify_EUR_loose'), n_x_ticks=5,
              x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
            p <- p + geom_point(data=df_1kg %>% filter(SUPER_POPULATION == "EUR"),
              mapping=aes_string(x=paste0('PC',i), y=paste0('PC',i+1)),
              inherit.aes=FALSE, shape=4, show.legend=FALSE)
            print(p)
        }

        aes <- aes_string(x=paste0('PC',i), y=paste0('PC',i+1), color='classification_strict')
        p <- create_pretty_scatter(df_classify, aes,
          save_figure=save_figures, file=paste0(PLOTS,'10_PC',i,'_PC',i+1,'_classify_EUR_strict'), n_x_ticks=5,
          x_label=paste0('Principal Component ',i), y_label=paste0('Principal Component ', i+1))
    }
}


