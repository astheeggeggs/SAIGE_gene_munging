library(data.table)
library(dplyr)
library(latex2exp)

source("QC_scripts/utils/pretty_plotting.r")

get_SAIGE_output_path <- function(phenotype="coronary_artery_disease", chr=20,
	minMAF="0", maxMAF="0.01", variant_class = "pLoF",
	SAIGE_results_dir="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/step2_results/") 
{
	file <- paste0(SAIGE_results_dir, phenotype, "_chr", chr, "_minMAF", minMAF, "_maxMAFforGroupTest_", maxMAF, "_", variant_class)
	if(file.exists(file)) {
		return(file)
	} else {
		stop(paste("File:", file, "does not exist"))
	}
}

read_and_create_qq <- function(
	phenotype = "coronary_artery_disease",
	minMAFs=c("0", "0", "0.01"),
	maxMAFs=c("0.01", "0.5", "0.5"),
	variant_classes=c(
		# "pLoF",
		# "damaging_missense",
		"other_missense"),#,
		# "synonymous"),
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
			cat(paste("\nPhenotype:", phenotype))
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

cts_phenotypes <- c(
	"Alanine_aminotransferase",
    "Albumin",
    "Alkaline_phosphatase",
    "Apolipoprotein_A",
    "Apolipoprotein_B",
    "Aspartate_aminotransferase",
    "C_reactive_protein",
    "Calcium",
    "Cholesterol",
    "Creatinine",
    "Cystatin_C",
    "Direct_bilirubin",
    "Gamma_glutamyltransferase",
    "Glucose",
    "Glycated_haemoglobin_HbA1c",
    "HDL_cholesterol",
    "IGF_1",
    "LDL_direct",
    "Lipoprotein_A",
    "Oestradiol",
    "Phosphate",
    "Rheumatoid_factor",
    "SHBG",
    "Testosterone",
    "Total_bilirubin",
    "Total_protein",
    "Triglycerides",
    "Urate",
    "Urea",
    "Vitamin_D",
    "Body_mass_index_BMI",
    "Hip_circumference",
    "Standing_height",
    "Waist_circumference",
    "Body_fat_percentage",
    "Visceral_adipose_tissue_volume_VAT",
    "Total_adipose_tissue_volume",
    "Abdominal_fat_ratio",
    "Liver_proton_density_fat_fraction_AMRA"
    )

binary_phenotypes <- c(
    # "colorectal_cancer", - Problem with this one in chromosome 21 - what's going on?
    "Trachea_bronchus_lung_cancer",
    "breast_cancer",
    "hypothalamic_amenorrhea",
    "POI",
    "dementia",
    "Alzheimers_disease",
    "depression",
    "autism",
    "ADHD",
    # "renal_failure", - Problem with this one - empty...what's the case count?!
    "coronary_artery_disease",
    "ischaemic_heart_disease",
    "stroke_hemorrhagic",
    # "stroke", - Problem with this one - empty...what's the case count?!
    "ischaemic_stroke",
    "chronic_obstructive_pulmonary_disease",
    "Crohns_disease",
    "IBD",
    "Cirrhosis",
    "NASH",
    # "NAFLD", - problem with NAFLD in chr21
    "psoriasis",
    "hyperandrogenism",
    "hematuria",
    "proteinuria",
    "acute_renal_failure",
    "chronic_kidney_disease",
    # "male_infertility", - problem with this one - empty...what's the case count?!
    "oligomenorrhea",
    "habitual_aborter",
    "female_infertility",
    "ectopic_pregnancy",
    "Preeclampsia",
    "GDM",
    "intrahepatic_cholestasis_in_pregnancy",
    "polycystic_kidney_disease",
    "T2D",
    "T1D",
    "GDM2",
    "kallmann_syndrome",
    "E230",
    "PCOS1",
    "PCOS2",
    "PCOS3",
    "PCOS4"
)

for (phenotype in cts_phenotypes) {
	read_and_create_qq(phenotype=phenotype)
}

for (phenotype in binary_phenotypes) {
	read_and_create_qq(phenotype=phenotype)
}

