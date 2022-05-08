# This is an example of an interactive singularity container for step 2.

saige_singularity_dir="/well/lindgren/dpalmer/saige_singularity"
version="1.0.7"
singularity shell --bind /well/lindgren/UKBIOBANK,/well/lindgren/dpalmer ${saige_singularity_dir}/saige.${version}.sif

# Navigate to the extdata folder of the new version of SAIGE (https://github.com/saigegit/SAIGE) and run this:
cd /well/lindgren/dpalmer/SAIGE/extdata

chr="20"
vcfFile="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr${chr}.vcf.bgz"
vcfFileIndex="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr${chr}.vcf.bgz.csi"
annotation="pLoF"
groupFile="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/ukb_wes_200k_filtered_chr${chr}_${annotation}_saige_gene.tsv.gz"
# I created the groupFile below using the helper function.
groupFile="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/tmp_groupFile.tsv"
phenotypeFolder="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k"
phenotype="Alanine_aminotransferase_residual"
GMMATmodelFile="${phenotypeFolder}/${phenotype}.rda"
varianceRatioFile="${phenotypeFolder}/${phenotype}_cate.varianceRatio.txt"
sparse_GRMFile="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
sparse_GRM_SampleID="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
outputfile="/well/lindgren/UKBIOBANK/dpalmer/${phenotype}_chr${chr}_${annotation}_test"

Rscript create_tmp_groupFile.R --groupFile

Rscript step2_SPAtests.R \
    --vcfFile ${vcfFile} \
    --vcfField GT \
    --vcfFileIndex ${vcfFileIndex} \
    --groupFile ${groupFile} \
    --minMAF 0 \
    --maxMAF_in_groupTest 0.01 \
    --GMMATmodelFile ${GMMATmodelFile} \
    --varianceRatioFile ${varianceRatioFile} \
    --sparseGRMFile ${sparse_GRMFile} \
    --sparseGRMSampleIDFile ${sparse_GRM_SampleID} \
    --is_output_moreDetails TRUE \
    --LOCO FALSE \
    --is_single_in_groupTest TRUE \
    --cateVarRatioMinMACVecExclude "0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude "1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --SAIGEOutputFile ${outputfile} \
    --annotation_in_groupTest="pLoF"


saige_singularity_dir="/well/lindgren/dpalmer/saige_singularity"
version="0.45"
singularity shell --bind /well/lindgren/UKBIOBANK,/well/lindgren/dpalmer ${saige_singularity_dir}/saige.${version}.sif

# Navigate to the extdata folder of the old version of SAIGE (https://github.com/weizhouUMICH/SAIGE) and run this:
cd /well/lindgren/dpalmer/SAIGE_old/extdata

chr="20"
annotation="pLoF"
vcfFile="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr${chr}.vcf.bgz"
vcfFileIndex="/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/final_mt/10_european.strict_filtered_chr${chr}.vcf.bgz.csi"
groupFile="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_variants_vep/200k/SAIGE_gene_input/ukb_wes_200k_filtered_chr${chr}_${annotation}_saige_gene.tsv.gz"
phenotypeFolder="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k"
phenotype="Alanine_aminotransferase_residual"
GMMATmodelFile="${phenotypeFolder}/${phenotype}.rda"
varianceRatioFile="${phenotypeFolder}/${phenotype}_cate.varianceRatio.txt"
sparseSigmaFile="${phenotypeFolder}/${phenotype}_cate.varianceRatio.txt_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseSigma.mtx"
outputfile="/well/lindgren/UKBIOBANK/dpalmer/${phenotype}_chr${chr}_${annotation}_test"

Rscript step2_SPAtests.R \
    --vcfFile ${vcfFile} \
    --vcfField GT \
    --vcfFileIndex ${vcfFileIndex} \
    --groupFile ${groupFile} \
    --chrom chr${chr} \
    --minMAF 0 \
    --maxMAFforGroupTest 0.01 \
    --sparseSigmaFile ${sparseSigmaFile} \
    --GMMATmodelFile ${GMMATmodelFile} \
    --varianceRatioFile ${varianceRatioFile} \
    --LOCO FALSE \
    --IsSingleVarinGroupTest TRUE \
    --cateVarRatioMinMACVecExclude "0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude "1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --SAIGEOutputFile ${outputfile} \
    --IsSparse TRUE \
    --IsOutputAFinCaseCtrl TRUE \
    --IsOutputHetHomCountsinCaseCtrl TRUE \
    --IsOutputNinCaseCtrl TRUE \
    --IsOutputlogPforSingle TRUE \
    --IsSingleVarinGroupTest TRUE
