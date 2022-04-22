library(data.table)
library(dplyr)
library(latex2exp)

source("QC_scripts/utils/pretty_plotting.r")

get_SAIGE_output_path <- function(phenotype="coronary_artery_disease", chr=20,
	minMAF="0", maxMAF="0.01", variant_class = "pLoF",
	SAIGE_results_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/step2_results/", continue=FALSE) 
{
	file <- paste0(SAIGE_results_dir, phenotype, "_chr", chr, "_minMAF", minMAF, "_maxMAFforGroupTest_", maxMAF, "_", variant_class)
	if(file.exists(file)) {
		return(file)
	} else {
		if (!continue) {
			stop(paste("File:", file, "does not exist"))
		} else {
			warning(paste("File:", file, "does not exist"))
		}
	}
}

read_and_create_qq <- function(
	phenotype = "coronary_artery_disease",
	minMAFs=c("0"),#, "0", "0.01"),
	maxMAFs=c("0.01"),#, "0.5", "0.5"),
	variant_classes=c(
		"pLoF",
		"damaging_missense",
		"other_missense",
		"synonymous"),
	outdir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/plots/"
	)
{
	if(length(minMAFs) != length(maxMAFs)) { stop("minMAFs and maxMAFs must be the same length") }
	
	# Create the QQ plots
	pdf(file=paste0(outdir, phenotype, "_qq.pdf"), width=4.5, height=3.5)
	for (i in 1:length(minMAFs)) {
		for (variant_class in variant_classes) {
			variant_class_plot <- gsub("\\|", ", ", variant_class)
			variant_class_plot <- gsub("_", " ", variant_class_plot)
			variant_class_filename <- gsub("\\|", "_", variant_class)
			dt <- list()
			for (chr in seq(1,22)) {
				cat(paste0("chromosome ", chr, "..."))
				dt[[chr]] <- fread(get_SAIGE_output_path(phenotype, chr, minMAFs[i], maxMAFs[i], variant_class_filename))
			}
			cols <- c("Gene", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT")
			dt <- rbindlist(dt)[, ..cols]
			dt[, Burden := ifelse(is.na(Pvalue_Burden), -log10(Pvalue), -log10(Pvalue_Burden))]
			dt[, SKAT := ifelse(is.na(Pvalue_SKAT), -log10(Pvalue), -log10(Pvalue_SKAT))]
			dt[, `SKAT-O` := -log10(Pvalue)]
			dt[, GeneID := gsub("_.*", "", Gene)]
			dt[, GeneName := gsub("^[^_]*_([^_]*)_.*", "\\1", Gene)]
			ribbon_p <- 0.95
			dt <- melt(dt, id.vars=c("GeneID", "GeneName"), measure.vars=c("Burden", "SKAT", "SKAT-O"), variable.name = "test" , value.name="Pvalue")
			dt <- data.table(
				dt %>% arrange(test, desc(Pvalue)) %>% select(-test),
				dt %>% group_by(test) %>% arrange(desc(Pvalue)) %>% 
					summarize(
						Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
						clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
						cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
						)
					)
			cat(paste("\nPhenotype:", phenotype, "\n"))
			create_pretty_qq_plot(
				plot_title=paste0(
					gsub("_", " ", paste0(toupper(substring(phenotype, 1,1)), substring(phenotype, 2))), "\n",
					"MAF: ", minMAFs[i], "-", maxMAFs[i], "\n",
					variant_class_plot
					),
				cex_labels=2,
				dt %>% mutate(labels=GeneName), aes(x=Pvalue_expected, y=Pvalue, color=test),
				save_figure=FALSE, n_to_include=10,
				x_label=TeX("$-\\log_{10}(\\mathit{p}_{expected})$"), 
				y_label=TeX("$-\\log_{10}(\\mathit{p}_{observed})$"),
				key_cols=c("test", "Pvalue"),
				aes_ribbon = aes(ymin=clower, ymax=cupper),
				width=170, height=120
			)
		}
	}
	dev.off()
}

read_and_create_gene_manhattan <- function(
	phenotype = "coronary_artery_disease",
	minMAFs=c("0"),#, "0", "0.01"),
	maxMAFs=c("0.01"),#, "0.5", "0.5"),
	variant_classes=c(
		"pLoF",
		"damaging_missense",
		"other_missense",
		"synonymous"),
	outdir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/plots/",
	gene_mapping="/well/lindgren/dpalmer/SAIGE_gene_munging/data/gene_mapping.txt.gz",
	tests_to_plot=c("Burden"),
	threshold=4,
	significance_T=0.05/20000
	)
{
	if(length(minMAFs) != length(maxMAFs)) { stop("minMAFs and maxMAFs must be the same length") }
	# Mapping file to obtain gene-start positions for each gene.
	# The default mapping file is taken from Biomart, build 38 with stable ID, start and end, chr, and gene name.
	if (grepl("*.gz$", gene_mapping)) {
		dt_gene <- fread(cmd = paste("zcat", gene_mapping))
	} else {
		dt_gene <- fread(gene_mapping)
	}

	names(dt_gene) <- c("GeneID", "GeneIDversion", "start", "stop", "chr", "GeneName_biomart")
	setkey(dt_gene, "GeneID")
	
	# Create the QQ plots
	pdf(file=paste0(outdir, phenotype, "_manhattan.pdf"), width=(230/25.4), height=(100/25.4))
	for (i in 1:length(minMAFs)) {
		for (variant_class in variant_classes) {
			variant_class_plot <- gsub("\\|", ", ", variant_class)
			variant_class_plot <- gsub("_", " ", variant_class_plot)
			variant_class_filename <- gsub("\\|", "_", variant_class)
			dt <- list()
			for (chr in seq(1,22)) {
				cat(paste0("chromosome ", chr, "..."))
				dt[[chr]] <- fread(get_SAIGE_output_path(phenotype, chr, minMAFs[i], maxMAFs[i], variant_class_filename))
			}
			cols <- c("Gene", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT")
			dt <- rbindlist(dt)[, ..cols]
			dt[, Burden := ifelse(is.na(Pvalue_Burden), -log10(Pvalue), -log10(Pvalue_Burden))]
			dt[, SKAT := ifelse(is.na(Pvalue_SKAT), -log10(Pvalue), -log10(Pvalue_SKAT))]
			dt[, `SKAT-O` := -log10(Pvalue)]
			dt[, GeneID := gsub("_.*", "", Gene)]
			dt[, GeneName := gsub("^[^_]*_([^_]*)_.*", "\\1", Gene)]
			setkey(dt, "GeneID")
			dt <- merge(dt, dt_gene)
			cat(paste("\nPhenotype:", phenotype, "\n"))
			for (t in tests_to_plot) {
				dt_tmp <- dt %>% select(chr, start, Pvalue, GeneName, matches(paste0("^", t, "$")))
				make_manhattan_plot(
					dt_tmp$chr, dt_tmp$start, dt$Pvalue, labels=dt$GeneName,
					title=paste0(
						gsub("_", " ", paste0(toupper(substring(phenotype, 1,1)), substring(phenotype, 2))), "\n",
						"MAF: ", minMAFs[i], "-", maxMAFs[i], "\n",
						variant_class_plot
						),
					threshold=threshold,
					significance_T=significance_T,
					save_figure=FALSE,
					print_p=TRUE
				)
			}
		}
	}
	dev.off()
}


# cts_phenotypes <- c(
# 	"Alanine_aminotransferase",
#     "Albumin",
#     "Alkaline_phosphatase",
#     "Apolipoprotein_A",
#     "Apolipoprotein_B",
#     "Aspartate_aminotransferase",
#     "C_reactive_protein",
#     "Calcium",
#     "Cholesterol",
#     "Creatinine",
#     "Cystatin_C",
#     "Direct_bilirubin",
#     "Gamma_glutamyltransferase",
#     "Glucose",
#     "Glycated_haemoglobin_HbA1c",
#     "HDL_cholesterol",
#     "IGF_1",
#     "LDL_direct",
#     "Lipoprotein_A",
#     "Oestradiol",
#     "Phosphate",
#     "Rheumatoid_factor",
#     "SHBG",
#     "Testosterone",
#     "Total_bilirubin",
#     "Total_protein",
#     "Triglycerides",
#     "Urate",
#     "Urea",
#     "Vitamin_D",
#     "Body_mass_index_BMI",
#     "Hip_circumference",
#     "Standing_height",
#     "Waist_circumference",
#     "Body_fat_percentage",
#     "Visceral_adipose_tissue_volume_VAT",
#     "Total_adipose_tissue_volume",
#     "Abdominal_fat_ratio",
#     "Liver_proton_density_fat_fraction_AMRA"
#     )

# Phenotypes
cts_phenotypes <- c(
    "Alanine_aminotransferase_residual",
    "Albumin_residual",
    "Alkaline_Phosphatase_residual",
    "Apolipoprotein_A_residual",
    "Apolipoprotein_B_residual",
    "Aspartate_aminotransferase_residual",
    "C_reactive_Protein_residual",
    "Calcium_residual", 
    "Cholesterol_residual", # What's going on here?
    "Creatinine_Serum_residual",
    "Creatinine_Urine_residual",
    "Cystatin_C_Serum_residual",
    "Direct_Bilirubin_residual",
    "Gamma_glutamyltransferase_residual",
    "Glucose_residual",
    "HbA1c_residual",
    "HDL_Cholesterol_residual",
    "IGF_1_residual",
    "Lipoprotein_A_residual",
    "Microalbumin_Urine_residual",
    "Oestradiol_residual",
    "Phosphate_residual",
    "Potassium_Urine_residual",
    "Rheumatoid_factor_residual",
    "SHBG_residual",
    "Sodium_Urine_residual",
    "Testosterone_residual",
    "Total_Bilirubin_residual",
    "Total_Protein_residual",
    "Triglyceride_residual",
    "Urate_residual",
    "Urea_residual",
    "Vitamin_D_residual",
    # "Alanine_aminotransferase_M_residual",
    # "Albumin_M_residual",
    # "Alkaline_Phosphatase_M_residual",
    # "Apolipoprotein_A_M_residual",
    # "Apolipoprotein_B_M_residual",
    # "Aspartate_aminotransferase_M_residual",
    # "C_reactive_Protein_M_residual",
    # "Calcium_M_residual",
    # "Cholesterol_M_residual",
    # "Creatinine_Serum_M_residual",
    # "Creatinine_Urine_M_residual",
    # "Cystatin_C_Serum_M_residual",
    # "Direct_Bilirubin_M_residual",
    # "Gamma_glutamyltransferase_M_residual",
    # "Glucose_M_residual",
    # "HbA1c_M_residual",
    # "HDL_Cholesterol_M_residual",
    # "IGF_1_M_residual",
    # "Lipoprotein_A_M_residual",
    # "Microalbumin_Urine_M_residual",
    # "Oestradiol_M_residual",
    # "Phosphate_M_residual",
    # "Potassium_Urine_M_residual",
    # "Rheumatoid_factor_M_residual",
    # "SHBG_M_residual",
    # "Sodium_Urine_M_residual",
    # "Testosterone_M_residual",
    # "Total_Bilirubin_M_residual",
    # "Total_Protein_M_residual",
    # "Triglyceride_M_residual",
    # "Urate_M_residual",
    # "Urea_M_residual",
    # "Vitamin_D_M_residual",
    # "Alanine_aminotransferase_F_residual",
    # "Albumin_F_residual",
    # "Alkaline_Phosphatase_F_residual",
    # "Apolipoprotein_A_F_residual",
    # "Apolipoprotein_B_F_residual",
    # "Aspartate_aminotransferase_F_residual",
    # "C_reactive_Protein_F_residual",
    # "Calcium_F_residual",
    # "Cholesterol_F_residual",
    # "Creatinine_Serum_F_residual",
    # "Creatinine_Urine_F_residual",
    # "Cystatin_C_Serum_F_residual",
    # "Direct_Bilirubin_F_residual",
    # "Gamma_glutamyltransferase_F_residual",
    # "Glucose_F_residual",
    # "HbA1c_F_residual",
    # "HDL_Cholesterol_F_residual",
    # "IGF_1_F_residual",
    # "Lipoprotein_A_F_residual",
    # "Microalbumin_Urine_F_residual",
    # "Oestradiol_F_residual",
    # "Phosphate_F_residual",
    # "Potassium_Urine_F_residual",
    # "Rheumatoid_factor_F_residual",
    # "SHBG_F_residual",
    # "Sodium_Urine_F_residual",
    # "Testosterone_F_residual",
    # "Total_Bilirubin_F_residual",
    # "Total_Protein_F_residual",
    # "Triglyceride_F_residual",
    # "Urate_F_residual",
    # "Urea_F_residual",
    # "Vitamin_D_F_residual",
    "BMI", # What's going on here?
    "BMI_imp",
    # "standing_height", # This one didn't fit a null model
    "VAT",
    "AFR",
    "AMRA",
    "body_fat_percentage",
    "waist_circumference",
    "hip_circumference",
    "WHR",
    "WHR_adj_BMI" #,
    # "WHR_adj_BMI_M"
    # "WHR_adj_BMI_F"
    )


binary_phenotypes <- c(
    # "BC_combined",
    # "BC_combined_primary_care",
    # "CAD_combined",
    # "CAD_combined_primary_care",
    # "COPD_combined",
    # "COPD_combined_primary_care",
    # "CLD_combined",
    # "CLD_combined_primary_care",
    # "CC_combined",
    # "CC_combined_primary_care",
    # "DEM_combined",
    # "DEM_combined_primary_care",
    # "INF_combined",
    # "INF_combined_primary_care",
    # "LC_combined",
    # "LC_combined_primary_care",
    # "NAFLD_combined",
    # "NAFLD_combined_primary_care",
    # "RF_combined",
    # "RF_combined_primary_care",
    # "RF_acute_combined",
    # "RF_acute_combined_primary_care",
    # "RF_chronic_combined",
    # "RF_chronic_combined_primary_care",
    # "STR_combined",
    # "STR_combined_primary_care",
    # "STR_hem_combined",
    # "STR_hem_combined_primary_care",
    # "STR_isc_combined",
    # "STR_isc_combined_primary_care",
    # "ADHD_combined",
    # "ADHD_combined_primary_care",
    # "AD_combined", # What's going on here?
    # "AD_combined_primary_care", # What's going on here?
    # "DEP_combined",
    # "DEP_combined_primary_care",
    # "HAM_combined",
    # "HAM_combined_primary_care",
    # "AUT_combined",
    # "AUT_combined_primary_care",
    # "CIRR_combined",
    # "CIRR_combined_primary_care",
    # "CD_combined",
    # "CD_combined_primary_care",
    # "IBD_combined",
    # "IBD_combined_primary_care",
    # "NASH_combined",
    # "NASH_combined_primary_care",
    # "PSOR_combined",
    # "PSOR_combined_primary_care",
    # "HEM_combined",
    # "HEM_combined_primary_care",
    # "PRO_combined",
    # "PRO_combined_primary_care",
    # "OLI_combined", # What's going on here?
    # "OLI_combined_primary_care",
    # "PRE_combined",
    # "PRE_combined_primary_care",
    # "EP_combined", # What's going on here?
    # "EP_combined_primary_care", # What's going on here?
    "PKD_combined",
    "PKD_combined_primary_care",
    "HAB_combined",
    "HAB_combined_primary_care",
    "POI_combined",
    "POI_combined_primary_care",
    "HYP_combined",
    "HYP_combined_primary_care",
    "DM_T1D",
    "DM_T2D",
    "DM_GD"
    )

# binary_phenotypes <- c(
#     "colorectal_cancer", # - Problem with this one in chromosome 21 - what's going on?
#     "Trachea_bronchus_lung_cancer",
#     "breast_cancer",
#     "hypothalamic_amenorrhea",
#     "POI",
#     "dementia",
#     "Alzheimers_disease",
#     "depression",
#     "autism",
#     "ADHD",
#     "renal_failure", #- Problem with this one - empty...what's the case count?!
#     "coronary_artery_disease",
#     "ischaemic_heart_disease",
#     "stroke_hemorrhagic",
#     "stroke",# - Problem with this one - empty...what's the case count?!
#     "ischaemic_stroke",
#     "chronic_obstructive_pulmonary_disease",
#     "Crohns_disease",
#     "IBD",
#     "Cirrhosis",
#     "NASH",
#     "NAFLD", # - problem with NAFLD in chr21
#     "psoriasis",
#     "hyperandrogenism",
#     "hematuria",
#     "proteinuria",
#     "acute_renal_failure",
#     "chronic_kidney_disease",
#     "male_infertility", #- problem with this one - empty...what's the case count?!
#     "oligomenorrhea",
#     "habitual_aborter",
#     "female_infertility",
#     "ectopic_pregnancy",
#     "Preeclampsia",
#     "GDM",
#     "intrahepatic_cholestasis_in_pregnancy",
#     "polycystic_kidney_disease",
#     "T2D",
#     "T1D",
#     "GDM2",
#     "kallmann_syndrome",
#     "E230",
#     "PCOS1",
#     "PCOS2",
#     "PCOS3",
#     "PCOS4"
# )

for (phenotype in cts_phenotypes) {
	read_and_create_qq(phenotype=phenotype)
	read_and_create_gene_manhattan(phenotype=phenotype)
}

for (phenotype in binary_phenotypes) {
	read_and_create_qq(phenotype=phenotype)
	read_and_create_gene_manhattan(phenotype=phenotype)
}

checking_outputs <- function(
	phenotypes,
	minMAFs=c("0"),#, "0", "0.01"),
	maxMAFs=c("0.01"),#, "0.5", "0.5"),
	variant_classes=c(
		"pLoF",
		"damaging_missense",
		"other_missense",
		"synonymous")
	) {

	for(phenotype in phenotypes) {
		# Checking that all results files have been created
		for (i in 1:length(minMAFs)) {
			for (variant_class in variant_classes) {
				variant_class_plot <- gsub("\\|", ", ", variant_class)
				variant_class_plot <- gsub("_", " ", variant_class_plot)
				variant_class_filename <- gsub("\\|", "_", variant_class)
				dt <- list()
				for (chr in c(seq(1,22), "X")) {
					file <- get_SAIGE_output_path(phenotype, chr, minMAFs[i], maxMAFs[i], variant_class_filename, continue=TRUE)
					if( !file.exists(file)) {
						print(file)
					}
				}
			}
		}
	}
}

checking_outputs(binary_phenotypes)
checking_outputs(cts_phenotypes)
