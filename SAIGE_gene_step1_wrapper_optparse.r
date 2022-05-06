#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(SAIGE)
require(optparse)

print(sessionInfo())

option_list <- list(
    make_option("--phenofile", type="character", default="",
        help = "Which phenotype file to pass")
    make_option("--phenotype", type="character", default="",
        help = "Which phenotype column to use in the phenotype file that was passed")
    make_option("--covarColList", type="character",
        default = paste0(
            "age,age2,age3,Inferred.Gender,genotyping.array,",
            "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21"),
        help = "Collection of covariates in the phenotype file to use")
    make_option("--sampleIDColinphenoFile", type="character", default = "eid",
        help = "Sample ID column in the phenotype file")
    make_option("--traitType", type="character", default = "quantitative",
        help = "What type of trait is the phenotype? Binary or quantitative?")
    make_option("--outdir", type="character",
        default = "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k")
)

# List of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

# Fred has run and created sparse GRMs using step 0 of SAIGE gene. 
# The resultant files are here:
# /well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input
# We use the European subset

sparse_GRM_folder <- "/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input"
sparse_GRM <- file.path(sparse_GRM_folder, "211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
sparse_GRM_SampleID <- paste0(sparse_GRM, ".sampleIDs.txt")
# plink file used to determine the sparse GRM - not actually used in the function call below.
plinkFile <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_1_22_including_rare_saige_input"

# Step 1: fitNULLGLMM
IsSparseKin <- TRUE
isCateVarianceRatio <- TRUE
if (opt$traitType == "quantitative") {
	invNormalize <- TRUE
} else {
	invNormalize <- FALSE
}

outputPrefix <- file.path(opt$outdir, opt$phenotype)
outputPrefix_varRatio <- file.path(opt$outdir, paste0(opt$phenotype, "_cate"))

# Now, let's fit the null model.
fitNULLGLMM(
    plinkFile = plinkFile,
    phenoFile = opt$phenofile,
    phenoCol = opt$phenotype,
    traitType = opt$traitType,
    invNormalize = invNormalize,
    covarColList = unlist(strsplit(opt$covarColList, split=",")),
    sampleIDColinphenoFile = opt$sampleIDColinphenoFile,
    LOCO = FALSE, # Note that we are not running LOCO - doing so will require 22 separate null models.
    outputPrefix = outputPrefix,
    outputPrefix_varRatio = outputPrefix_varRatio,
    IsOverwriteVarianceRatioFile = TRUE,
    IsSparseKin = TRUE,
    sparseGRMFile = sparse_GRM,
    sparseGRMSampleIDFile = sparse_GRM_SampleID,
    isCateVarianceRatio = TRUE,
    cateVarRatioMinMACVecExclude = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    cateVarRatioMaxMACVecInclude = c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    useSparseGRMtoFitNULL = TRUE)
