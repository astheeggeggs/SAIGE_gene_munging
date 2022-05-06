#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
library(SAIGE)

BLASctl_installed <- require(RhpcBLASctl)
library(optparse)
library(data.table)
library(methods)
print(sessionInfo())

option_list <- list(
    make_option("--vcfFile", type="character",
        default="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr18.vcf.bgz",
        help="Path to vcf file")
    make_option("--vcfFileIndex", type="character",
        default="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr18.vcf.bgz.csi",
        help="Path to index for vcf file by tabix, '.tbi' by 'tabix -p vcf file.vcf.gz'")
    make_option("--phenotype", type="character", default="Body_mass_index_BMI",
        help = cat("Phenotype name. Must match the phenotype name in the",
        "filepath to GMMATmodelFile, varianceRatioFile and sparseSigmaFile"))
    make_option("--groupFile", type="character",
        default="/well/lindgren/dpalmer/test_MC4R.tsv",
        help = cat("Path to the file containing the group information",
        "for gene-based tests. Each line is for one gene/set of",
        "variants. The first element is for gene/set name. The rest of",
        "the line is for variant ids included in this gene/set. For",
        "vcf/sav, the genetic marker ids are in the format",
        "chr:pos_ref/alt. For bgen, the genetic marker ids should",
        "match the ids in the bgen file. Each element in the line is",
        "seperated by tab."))
    make_option("--chr", type="character", default="chr18",
        help = "Which chromosome are you analysing?")
    make_option("--phenotypeFolder", type="character",
        default = "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k",
        help = cat("Absolute filepath to the folder containing the",
        "GMMATmodelFile and varianceRatioFile"))
    make_option("--minMAF", type="numeric", default = 0,
        help= cat("Maximum minor allele frequency of markers to test in",
        "group test."))
    make_option("--maxMAF_in_groupTest", type="character", default = "0.5",
        help= cat("Maximum minor allele frequency of markers to test in",
        "group test."))
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# Now, let's run the gene based tests and single variant analyses
# Inputs from step 1.
GMMATmodelFile <- file.path(opt$phenotypeFolder, paste0(opt$phenotype, ".rda"))
varianceRatioFile <- file.path(opt$phenotypeFolder, paste0(opt$phenotype, "_cate.varianceRatio.txt"))

sparse_GRM_folder <- "/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input"
sparse_GRM <- file.path(sparse_GRM_folder, "211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
sparse_GRM_SampleID <- paste0(sparse_GRM, ".sampleIDs.txt")

# Output locations
dir.create(file.path(opt$phenotypeFolder, "step2_results"), showWarnings=FALSE)
annotation <- gsub("^.*chr[0-9XY]+_(.*)_saige_gene.tsv.gz", "\\1", opt$groupFile)
SAIGEOutputFile <- file.path(
    opt$phenotypeFolder,
    "step2_results",
    paste0(opt$phenotype, "_", opt$chr, "_minMAF", opt$minMAF,
        "_maxMAFforGroupTest_", opt$maxMAFforGroupTest, "_", annotation)
    )

SPAGMMATtest(
    vcfFile = opt$vcfFile,
    vcfFileIndex = opt$vcfFileIndex,
    vcfField = "GT",
    chrom = opt$chr,
    min_MAF = opt$minMAF,
    maxMAF_in_groupTest = opt$maxMAFforGroupTest,
    min_Info = 0,
    GMMATmodelFile = GMMATmodelFile,
    varianceRatioFile = varianceRatioFile,
    SAIGEOutputFile = SAIGEOutputFile,
    is_output_moreDetails = TRUE,
    LOCO = FALSE,
    condition = "",
    sparseGRMFile = sparse_GRM,
    sparseGRMSampleIDFile = sparse_GRM_SampleID,
    groupFile = opt$groupFile,
    r.corr = 0,
    is_single_in_groupTest = TRUE,
    cateVarRatioMinMACVecExclude = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    cateVarRatioMaxMACVecInclude = c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    dosage_zerod_cutoff = 0.2,
    MACCutoff_to_CollapseUltraRare = 10,
)


