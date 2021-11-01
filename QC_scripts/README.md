## 00_tidy_bim_bed_fam_determine_EUR.r
*Inputs*

- UKBB plink files from UK Biobank
- 1000G ped file
- output from `00_plink_file_preparation.r` (UKBB fam file in format which can be read by plink).

*Outputs*

- `ukb_autosomes_combined{.bim,.bed,.fam}`
- `ukb_autosomes_combined_EUR{.bim, .bed, .fam}`
- `1000G_phase3_common_norel_EUR.tsv` (European samples)
- `final_EUR_list.tsv` (Non-Finnish European samples)
- plots

This script runs a bunch of munging, estimates superpopulations of the individuals by projecting into the PC spaces defined by 1000 genomes samples. Note that ancestry is determined using the genotype data here.
We then create a random forest classifier and estimate the EUR samples. Following that restriction, we repeat the process, defining PCs in 1000G Europeans, to define the non-Finnish Europeans and remove any putatative Finns from the dataset. The output is a collection of 'Non-Finnish European' samples in the genotype data, which we will filter to as a later step in the QC pipeline.

## 00_create_phenotype_ht.py
*Inputs*

- `QC_phenotypes.tsv.gz` (output from `00_auxiliary_sample_information_preparation.r`)

*Outputs*

- `phenotypes.ht`

Reads in and merges the phenotype data into a hail table and writes a ht.

## 01_create_hardcalls_mt.py
*Inputs*

- `ukb_wes_{TRANCHE}_filtered_chr{CHR}.mt`

*Outputs*

- `ukb_wes_{TRANCHE}_filtered_hardcalls_chr{CHR}.mt`

This step reads in the mt for each chromosome, selects just the GT entries of the mt and write the hardcalls to disk after repartitioning the mt. This is to speed up subsequent steps where we don't require dosages etc.

## 02_prefilter_variants.py
*Inputs*

- `ukb_wes_{TRANCHE}_filtered_hardcalls_chr{CHR}.mt`

*Outputs*

- `02_prefilter_chr{CHR}.keep.variant.ht`

In this step, we remove invariant rows. Note that invariant rows (variants) will already have been removed from Nik's earlier QC pipeline, so this is likely redundant.

## 03_initial_sample_qc.py
*Inputs*

- `ukb_wes_{TRANCHE}_filtered_chr{CHR}.mt`
- `02_prefilter_chr{CHR}.keep.variant.ht`

*Outputs*

- `03_chr{CHR}_initial_sample_qc.tsv.bgz`

In this step, we run the sample QC function in hail to determine sample qc metrics specific to chromosomes. Note that when we move to the cloud, we can streamline this and the next step to avoid the need for a loop to integrate metrics over chromosomes. Note - autosomes only is fine here.

## 03_initial_sample_qc_plot.r
*Inputs*

- `03_chr{CHR}_initial_sample_qc.tsv.bgz`

*Outputs*

- `03_initial_sample_qc.tsv`
- plots

This function is called by 03_initial_sample_qc_filter.r. Firstly, the R script gathers all the sample QC metrics across the autosomes in the appropriate manner. These metrics are then wrtten to a new .tsv alongside the chromosome specific sample metrics. We then create a series of plots, and define an empirical hard cutoff for a series of the metrics.

## 03_initial_sample_qc_filter.r
*Inputs*

- `03_initial_sample_qc.tsv`

*Outputs*

- `03_initial_qc.keep.sample_list`
- `03_sample_count.tsv`

This function runs the previous plotting functions (03_initial_sample_qc_plot.r), and defines the passing subset of samples, and writes it. We also create and write a summary table of the number of samples removed due to filtering so far in the pipeline.

## 04_prune_genotyped_snps.sh

*Inputs*

- UKBB plink files from UK Biobank

*Outputs*

- `ukb_snp_chrX_pruned.prune.{in, out}`

This first step takes the X chromosome, and LD prunes to define a collection of pseudo-independent SNPs for subsequent F-statistic evaluation. We filter to the collection of samples with exome sequence data available to speed things.

## 04_impute_sex.py

*Inputs*

- UKBB plink files for X and Y
- `ukb_snp_chrX_pruned.prune.in`
- `QC_phenotypes.tsv.gz` (output from `00_auxiliary_sample_information_preparation.r`)
- `03_initial_qc.keep.sample_list`
- `phenotypes.ht`

*Outputs*

- `04_imputesex.ht`
- `04_imputesex.tsv.bgz`
- `04_ycalled.tsv.bgz`

Here, we go back to the genotype data and check for sex swaps. This is likely already done, but we're being ultra careful. To do this we read in the sex chromosome genotype data, and determine the F-statistic for samples using the X chromosome, and check the number of non-missing allele counts on the Y.

## 04_impute_sex_plot.r

*Inputs*

- `04_imputesex.tsv.bgz`
- `04_ycalled.tsv.bgz`

*Outputs*

- `04_sexcheck.remove.sample_list`
- plots

Plot the output of 04_impute_sex.py. We plot the distribution of the F-statistic on the X, and define a cutoff for sex labelling. We also plot the X F-statistic against the number of reads on the Y chromosome. After adding genetically defined sex, we compare to the self assigned sex in the phenotype file and remove mismatches.

## 05_annotate_variants_vep.py

*Inputs*

- `vep_env.json`
- `gnomad.exomes.r2.1.1.sites.{CHR}.liftover_grch38.vcf.bgz`
- `ukb_wes_{TRANCHE}_filtered_chr{CHR}.mt`

*Outputs*

- `ukb_wes_{TRANCHE}_filtered_chr{CHR}_vep_qc.ht`

We annotate the mt using VEP, and run the process_consequences function. We also annotate with whether the variant is present in gnomAD v2 (following liftover).

## 07_ultra_rare_counts.py

*Inputs*

- `ukb_wes_{TRANCHE}_filtered_chr{CHR}.mt`
- `ukb_wes_{TRANCHE}_filtered_chr{CHR}_vep_qc.ht`
- `phenotypes.ht`
- `02_prefilter_chr{CHR}.keep.variant.ht`
- `03_initial_qc.keep.sample_list`
- `04_sexcheck.remove.sample_list`

*Outputs*

- `07_URVs_chr{CHR}.tsv.bgz`

Read in and filter down the mt to the subset of samples we've restricted to so far, and annotate the mt with the VEP and gnomAD annotations from the previous step. Determine counts of variants in each sample the are singletons in the dataset and not present in gnomAD. Note again, this is run separately for each chromosome, we can streamline this and the next step to avoid the need for a loop to integrate metrics over chromosomes.

## 07_ultra_rare_counts_plot.r

*Inputs*

- `07_URVs_chr{CHR}.tsv.bgz`

*Outputs*

- `07_URVs.tsv`
- `07_URV.remove.sample_list`
- plots

Firstly, the R script gathers all the sample QC metrics across the autosomes in the appropriate manner. We then create a series of plots to look at the distribution of ultra-rare singleton counts, and assign a hard filter where appropriate. This threshold will be used in the next step to remove samples.

## 08_final_variant_qc.py

*Inputs*

- `ukb_wes_{TRANCHE}_filtered_hardcalls_chr{CHR}.mt`
- `phenotypes.ht`
- `02_prefilter_chr{CHR}.keep.variant.ht`
- `03_initial_qc.keep.sample_list`
- `04_sexcheck.remove.sample_list`
- `04_imputesex.ht`
- `07_URV.remove.sample_list`

*Outputs*

- `08_final_qc.variants_chr{CHR}.tsv.bgz`

Filter samples down to the collection of samples we've restriced to so far, and determine variant metrics. Note that we need to run this on the X. 

## 08_final_variant_qc_plot.r

*Inputs*

- `08_final_qc.variants_chr{CHR}.tsv.bgz`

*Outputs*

- `08_final_qc.variants.tsv`
- plots

Firstly, the R script gathers all the sample QC metrics across the chromosomes in the appropriate manner. We then plot the resultant variant metrics determined above - call rate and pHWE, and determine empirical cutoffs.

## 08_final_variant_qc_filter.r

*Inputs*

- `08_final_qc.variants.tsv`

*Outputs*

- `08_final_qc.keep.variant_list`
- `08_variant_count.tsv`

This function runs the previous plotting functions (08_final_variant_qc_plot.r). We write the resultant variants to file and also save a summary table of the number of variants removed due to filtering so far in the pipeline.

## 09_final_sample_qc.py

*Inputs*

- `ukb_wes_{TRANCHE}_filtered_hardcalls_chr{CHR}.mt`
- `phenotypes.ht`
- `02_prefilter_chr{CHR}.keep.variant.ht`
- `03_initial_qc.keep.sample_list`
- `04_sexcheck.remove.sample_list`
- `04_imputesex.ht`
- `07_URV.remove.sample_list`
- `08_final_qc.keep.variant_list`

*Outputs*

- `09_final_qc_chr{CHR}.before.samples.tsv.bgz`
- `09_final_qc_chr{CHR}.after.samples.tsv.bgz`

Here, we read in the mt and determine sample metrics before and after the variant filters in step 09. 

## 09_final_sample_qc_plot.r

*Inputs*

- `09_final_qc_chr{CHR}.{before,after}.samples.tsv.bgz`

*Outputs*

- `09_final_qc.{before,after}.samples.tsv`
- plots

After determining sample metrics in 09_final_sample.py for each chromosome we amalgamate in the appropriate manner, and plot the distribution of the metrics across the automsomes, and highlight samples that lie > 3 standard deviations away from the mean of that metric.

## 09_final_filter_sample_qc.r

*Inputs*

- `09_final_qc.{before,after}.samples.tsv`

*Outputs*

- `09_final_sample.summary.tsv`
- `09_final_qc.keep.sample_list`

Following the plotting above, we determine the collection of samples to be removed and write the final sample list to file. We also write a summary table which details the number of samples removed at each step.

## 10_create_qc_mt.py

*Inputs*

- `ukb_wes_{TRANCHE}_filtered_chr{CHR}.mt`
- `08_final_qc.keep.variant_list`
- `09_final_qc.keep.sample_list`

*Outputs*

- `10_european.strict_filtered_chr{CHR}.mt`

This is the final step! We filter down to the final mt by pruning out the appropriates samples and variants. We run PCA on the samples, so that we can check that the resultant samples are relatively homogeneous in PC space with no clear outliers, and write a full mt and hard callled mt to disk. Note here that we also need to output the vcfs as these are required by SAIGE-gene.

## 10_create_qc_vcf.py

*Inputs*

- `10_european.strict_filtered_chr{CHR}.mt`

*Outputs*

- `10_european.strict_filtered_chr{CHR}.vcf.bgz`

Read in the saved matrix table, and export a compressed vcf.

## 10_create_vcf_index.sh

*Inputs*

- `10_european.strict_filtered_chr{CHR}.vcf.bgz`

*Outputs*

- `10_european.strict_filtered_chr{CHR}.vcf.bgz.csi`

Index the vcf, ensuring that indexing is `.csi`, to ensure compatability with SAIGE.
