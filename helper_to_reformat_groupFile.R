groupFile <- "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/ukb_wes_200k_filtered_chr20_pLoF_saige_gene.tsv.gz"
annotations <- "pLoF"

incon <- gzfile(groupFile)
lines <- readLines(incon)
lines <- rep(lines, each=2)

for (i in seq(2, length(lines), by=2)) {
	lines[i] <- gsub("[_/]", ":", lines[i])
	lines[i] <- gsub("\t[^\t]+", "\tpLoF", lines[i])
	lines[i] <- gsub("(^[^\t]+)", "\\1\tanno", lines[i])
}

for (i in seq(1, length(lines), by=2)) {
	lines[i] <- gsub("[_/]", ":", lines[i])
	lines[i] <- gsub("(^[^\t]+)", "\\1\tvar", lines[i])
}

writeLines(lines, "/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/tmp_groupFile.tsv")
