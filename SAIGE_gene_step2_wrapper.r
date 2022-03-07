#!/usr/bin/env Rscript

# Ensure that the SAIGE environment has been activated:
# conda activate /well/lindgren/users/mmq446/conda/skylake/envs/RSAIGE

library(SAIGE)#, lib.loc='/well/lindgren/flassen/software/tmp/') 
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--vcfFile",
    default="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr18.vcf.bgz",
    help="Path to vcf file")
parser$add_argument("--vcfFileIndex",
    default="/well/lindgren/UKBIOBANK/nbaya/wes_200k/ukb_wes_qc/data/filtered/ukb_wes_200k_filtered_chr18.vcf.bgz.csi",
    help="Path to index for vcf file by tabix, '.tbi' by 'tabix -p vcf file.vcf.gz'")
parser$add_argument("--phenotype", default="Body_mass_index_BMI",
    help = cat("Phenotype name. Must match the phenotype name in the",
    "filepath to GMMATmodelFile, varianceRatioFile and sparseSigmaFile"))
parser$add_argument("--groupFile",
    # default="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/ukb_wes_200k_filtered_chr22_pLoF_saige_gene.tsv.gz",
    default="/well/lindgren/dpalmer/test_MC4R.tsv",
    help = cat("Path to the file containing the group information",
    "for gene-based tests. Each line is for one gene/set of",
    "variants. The first element is for gene/set name. The rest of",
    "the line is for variant ids included in this gene/set. For",
    "vcf/sav, the genetic marker ids are in the format",
    "chr:pos_ref/alt. For bgen, the genetic marker ids should",
    "match the ids in the bgen file. Each element in the line is",
    "seperated by tab."))
parser$add_argument("--chr", default="chr18",
    help = "Which chromosome are you analysing?")
parser$add_argument("--phenotypeFolder",
    default = "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k",
    help = cat("Absolute filepath to the folder containing the",
    "GMMATmodelFile, varianceRatioFile, and sparseSigmaFile"))
parser$add_argument("--minMAF", default = 0,
    help= cat("Maximum minor allele frequency of markers to test in",
    "group test."))
parser$add_argument("--maxMAFforGroupTest", default = 0.5,
    help= cat("Maximum minor allele frequency of markers to test in",
    "group test."))

args <- parser$parse_args()

# Now, let's run the gene based tests and single variant analyses
# Inputs from step 1.
GMMATmodelFile <- file.path(args$phenotypeFolder, paste0(args$phenotype, ".rda"))
varianceRatioFile <- file.path(args$phenotypeFolder, paste0(args$phenotype, "_cate.varianceRatio.txt"))
sparseSigmaFile <- grep(paste0(args$phenotype, "_cate.varianceRatio.txt_"), dir(args$phenotypeFolder, full.names=TRUE), value=TRUE)

# Output locations
dir.create(file.path(args$phenotypeFolder, "step2_results"), showWarnings=FALSE)
annotation <- gsub("^.*chr[0-9XY]+_(.*)_saige_gene.tsv.gz", "\\1", args$groupFile)
SAIGEOutputFile <- file.path(
    args$phenotypeFolder,
    "step2_results",
    paste0(args$phenotype, "_", args$chr, "_minMAF", args$minMAF,
        "_maxMAFforGroupTest_", args$maxMAFforGroupTest, "_", annotation)
    )

SPAGMMATtest(
    vcfFile = args$vcfFile,
    vcfFileIndex = args$vcfFileIndex,
    vcfField = "GT",
    chrom = args$chr,
    IsDropMissingDosages = FALSE,
    minMAF = as.numeric(args$minMAF),
    maxMAFforGroupTest = as.numeric(args$maxMAFforGroupTest),
    minInfo = 0,
    GMMATmodelFile = GMMATmodelFile,
    varianceRatioFile = varianceRatioFile,
    SAIGEOutputFile = SAIGEOutputFile,
    IsSparse = TRUE,
    IsOutputAFinCaseCtrl = TRUE,
    IsOutputHetHomCountsinCaseCtrl = TRUE,
    IsOutputNinCaseCtrl = TRUE,
    IsOutputlogPforSingle = TRUE,
    LOCO = FALSE,
    condition = "",
    sparseSigmaFile = sparseSigmaFile,
    groupFile = args$groupFile,
    kernel = "linear.weighted",
    method = "optimal.adj",
    weights.beta.rare = c(1, 25),
    weights.beta.common = c(1, 25),
    weightMAFcutoff = 0.01,
    weightsIncludeinGroupFile = FALSE,
    weights_for_G2_cond = NULL,
    r.corr = 0,
    IsSingleVarinGroupTest = TRUE,
    cateVarRatioMinMACVecExclude = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    cateVarRatioMaxMACVecInclude = c(1.5, 2.5, 3.5, 4.5, 5.5, 10.5, 20.5),
    dosageZerodCutoff = 0.2,
    IsOutputPvalueNAinGroupTestforBinary = FALSE,
    IsAccountforCasecontrolImbalanceinGroupTest = TRUE,
    method_to_CollapseUltraRare = "absence_or_presence",
    MACCutoff_to_CollapseUltraRare = 10,
    DosageCutoff_for_UltraRarePresence = 0.5
)

