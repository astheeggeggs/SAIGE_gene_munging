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


minMAF = c(0, 0, 0.01),
    maxMAFforGroupTest = c(0.01, 0.5, 0.5)

read_and_create_qq <- function(
	phenotype = "coronary_artery_disease",
	minMAFs=c("0", "0", "0.01"),
	maxMAFs=c("0.01","0.5", "0.5"),
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

phenotypes <- c(
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
    "Body_fat_percentage"
    )

for (phenotype in phenotypes) {
	read_and_create_qq(phenotype=phenotype)
}
