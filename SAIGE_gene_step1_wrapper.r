#!/usr/bin/env Rscript

# Ensure that the SAIGE environment has been activated:
# conda activate /well/lindgren/users/mmq446/conda/skylake/envs/RSAIGE

library(SAIGE, lib.loc='/well/lindgren/flassen/software/tmp/') 
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--phenofile", required=TRUE,
    help = "Which phenotype file to pass")
parser$add_argument("--phenotype", required=TRUE,
    help = "Which phenotype column to use in the phenotype file that was passed")
parser$add_argument("--covarColList",
    default = paste0(
        "age,age2,age3,Inferred.Gender,genotyping.array,",
        "PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21"),
    help = "Collection of covariates in the phenotype file to use")
parser$add_argument("--sampleIDColinphenoFile",
    default = "ID",
    help = "Sample ID column in the phenotype file")
parser$add_argument("--traitType",
    default = "quantitative",
    help = "What type of trait is the phenotype? Binary or quantitative?")
parser$add_argument("--outdir",
    default = "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k")
       
args <- parser$parse_args()

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
if (args$traitType == "quantitative") {
	invNormalize <- TRUE
} else {
	invNormalize <- FALSE
}

outputPrefix <- file.path(args$outdir, args$phenotype)
outputPrefix_varRatio <- file.path(args$outdir, paste0(args$phenotype, "_cate"))

# Now, let's fit the null model.
fitNULLGLMM(
    plinkFile = plinkFile,
    phenoFile = args$phenofile,
    phenoCol = args$phenotype,
    traitType = args$traitType,
    invNormalize = invNormalize,
    covarColList = args$covarColList,
    sampleIDColinphenoFile = args$sampleIDColinphenoFile,
    LOCO = FALSE, # Note that we are not running LOCO - doing so will require 22 separate null models.
    outputPrefix = outputPrefix,
    outputPrefix_varRatio = outputPrefix_varRatio,
    IsOverwriteVarianceRatioFile = TRUE,
    IsSparseKin = TRUE,
    sparseGRMFile = sparse_GRM,
    sparseGRMSampleIDFile = sparse_GRM_SampleID,
    isCateVarianceRatio = TRUE,
    cateVarRatioIndexVec = NULL,
    cateVarRatioMinMACVecExclude = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    cateVarRatioMaxMACVecInclude = c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    useSparseGRMtoFitNULL = TRUE)
