library(data.table)
library(dplyr)

dt <- fread(cmd = "gzcat debugging.tsv.bgz") %>% 
select(
	-excessHet_vcf.filters,
	-Nik_QC.filters,
	-Nik_QC.info.ExcessHet,
	-split_multi.info.AC,
	-split_multi.info.AN
	)

dt <- dt %>% mutate(
	variant = paste0(
		locus, ":",
		gsub("^[^A-Z]*([A-Z]+)[^A-Z]*([A-Z]+)[^A-Z]+", "\\1:\\2", alleles)
		)
) %>% select(-locus, -alleles)

dt <- dt %>% mutate(
	excessHet_vcf.info.AC = as.integer(gsub("\\[([0-9]+)\\]","\\1", excessHet_vcf.info.AC)),
	excessHet_vcf.info.AQ = as.integer(gsub("\\[([0-9]+)\\]","\\1", excessHet_vcf.info.AQ)),
	excessHet_vcf.info.AF = gsub("\\[(.*)\\]","\\1", excessHet_vcf.info.AF),

	split_multi.filters = gsub("\\[\"(.*)\"\\]","\\1", split_multi.filters),
	split_multi.info.AQ = as.integer(gsub("\\[([0-9]+)\\]","\\1", split_multi.info.AQ)),
	split_multi.info.AF = gsub("\\[(.*)\\]","\\1", split_multi.info.AF),

	Nik_QC.info.AQ = as.integer(gsub("\\[([0-9]+)\\]","\\1", Nik_QC.info.AQ)),
	Nik_QC.info.AC = as.integer(gsub("\\[([0-9]+)\\]","\\1", Nik_QC.info.AC)),
	Nik_QC.info.AN = as.integer(gsub("\\[([0-9]+)\\]","\\1", Nik_QC.info.AN)),
	Nik_QC.info.AF = gsub("\\[(.*)\\]","\\1", Nik_QC.info.AF),

	split_multi.in_variant_list = ifelse(is.na(split_multi.in_variant_list), FALSE, split_multi.in_variant_list),
	genebass.in_variant_list = ifelse(is.na(genebass.in_variant_list), FALSE, genebass.in_variant_list),
	Nik_QC.in_variant_list = ifelse(is.na(Nik_QC.in_variant_list), FALSE, Nik_QC.in_variant_list),
	Duncan_QC.in_variant_list = ifelse(is.na(Duncan_QC.in_variant_list), FALSE, Duncan_QC.in_variant_list),
	) %>% mutate(
	excessHet_vcf.info.AF = as.numeric(gsub("E", "e", excessHet_vcf.info.AF)),
	split_multi.info.AF = as.numeric(gsub("E", "e", split_multi.info.AF)),
	Nik_QC.info.AF = as.numeric(gsub("E", "e", Nik_QC.info.AF))
	)

