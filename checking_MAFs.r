library(data.table)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)

dt <- list()
for (chr in c(seq(1,22), "X")) {
	cat(paste0(chr, "... "))
	dt[[as.character(chr)]] <- fread(cmd = paste0("gzcat debugging_chr", chr, ".tsv.bgz"))
}
dt <- rbindlist(dt)
dt[, variant:=paste(locus, ref, alt, sep=":")]

# First, just those that are not present in Nik's initial QC
# Then, those that are not present after my additional QC

variants <- list(
	`Raw VCF` = dt$variant,
	`genebass` = (dt %>% filter(genebass_in_variant_list))$variant,
	`Initial QC` = (dt %>% filter(!is.na(Nik_QC_AF)))$variant,
	`Subsequent QC` = (dt %>% filter(Duncan_QC_in_variant_list))$variant
	)

venn <- Venn(variants)
data <- process_data(venn)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  # 4. region label layer
  geom_sf_label(aes(label = count), data = venn_region(data)) +
  theme_void()

dt_check <- dt %>% filter(genebass_in_variant_list, is.na(Nik_QC_AF))

variants_sub <- list(
	`in our LCR` = (dt_check %>% filter(!not_in_lcr))$variant,
	`not in target +50bp` = (dt_check %>% filter(!is_in_target_padded_50bp))$variant,
	`invariant` = (dt_check %>% filter(is.na(excessHet_vcf_AC) & is_in_target_padded_50bp & not_in_lcr))$variant,
	`genebass` = (dt_check %>% filter(genebass_in_variant_list))$variant,
	`Excesshet` = (dt_check %>% filter(!excesshet_leq_54.69))$variant
	)
venn_sub <- Venn(variants_sub)

data_sub <- process_data(venn_sub)
ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data_sub)) +
  # 3. set label layer
  geom_sf_text(aes(label = name), data = venn_setlabel(data_sub)) +
  # 4. region label layer
  geom_sf_label(aes(label = count), data = venn_region(data_sub)) +
  theme_void()
