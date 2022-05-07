# This is an example of an interactive singularity container for step 2.

saige_singularity_dir="/well/lindgren/dpalmer/saige_singularity"
version="1.0.5"
singularity shell --bind /well/lindgren/UKBIOBANK,/well/lindgren/dpalmer ${saige_singularity_dir}/saige.${version}.sif

# Navigate to the extdata folder and run this:

chr="20"
vcfFile="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr${chr}.vcf.bgz"
vcfFileIndex="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr${chr}.vcf.bgz.csi"
annotation="pLoF"
groupFile="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/ukb_wes_200k_filtered_chr${chr}_${annotation}_saige_gene.tsv.gz"
phenotypeFolder="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k"
phenotype="Alanine_aminotransferase_residual"
GMMATmodelFile="${phenotypeFolder}/${phenotype}.rda"
varianceRatioFile="${phenotypeFolder}/${phenotype}_cate.varianceRatio.txt"
sparse_GRMFile="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
sparse_GRM_SampleID="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
outputfile="/well/lindgren/UKBIOBANK/dpalmer/${phenotype}_chr${chr}_${annotation}_test"

Rscript step2_SPAtests.R \
    --vcfFile ${vcfFile} \
    --vcfField GT \
    --vcfFileIndex ${vcfFileIndex} \
    --groupFile ${groupFile} \
    --chr "chr${chr}" \
    --minMAF 0 --maxMAF_in_groupTest 0.01 \
    --GMMATmodelFile ${GMMATmodelFile} \
    --varianceRatioFile ${varianceRatioFile} \
    --sparseGRMFile ${sparse_GRMFile} --sparseGRMSampleIDFile ${sparse_GRM_SampleID} \
    --is_output_moreDetails TRUE \
    --LOCO FALSE \
    --groupFile ${groupFile} \
    --is_single_in_groupTest TRUE \
    --cateVarRatioMinMACVecExclude "0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude "1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --SAIGEOutputFile ${outputfile}
