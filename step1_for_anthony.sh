# This is an example of an interactive singularity container for step 1.

saige_singularity_dir="/well/lindgren/dpalmer/saige_singularity"
version="1.0.5"
singularity shell --bind /well/lindgren/UKBIOBANK,/well/lindgren/dpalmer ${saige_singularity_dir}/saige.${version}.sif

# Navigate to the extdata folder of the new version of SAIGE (https://github.com/saigegit/SAIGE) and run this:
cd /well/lindgren/dpalmer/SAIGE/extdata

phenofile="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k/UKBB_WES200k_filtered_cts_dec2021_phenotypes.tsv.gz"
phenotype="Alanine_aminotransferase_residual"
covars="age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sex,sequencing.batch"
ID_col="eid"
trait_type="quantitative"

plinkFile="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_1_22_including_rare_saige_input"
outputPrefix="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/${phenotype}"
outputPrefix_varRatio="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/${phenotype}_cate"
sparse_GRMFile="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
sparse_GRM_SampleID="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

Rscript step1_fitNULLGLMM.R \
    --plinkFile ${plinkFile} \
    --phenoFile ${phenofile} \
    --phenoCol ${phenotype} \
    --covarColList ${covars} \
    --sampleIDColinphenoFile ${ID_col} \
    --traitType ${trait_type} \
    --invNormalize TRUE \
    --LOCO FALSE \
    --outputPrefix ${outputPrefix} --outputPrefix_varRatio ${outputPrefix_varRatio} \
    --IsOverwriteVarianceRatioFile TRUE \
    --sparseGRMFile ${sparse_GRMFile} --sparseGRMSampleIDFile ${sparse_GRM_SampleID} \
    --isCateVarianceRatio TRUE \
    --cateVarRatioMinMACVecExclude "0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude "1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --useSparseGRMtoFitNULL TRUE 

# This is an example of an interactive singularity container for step 1 - using the old version of SAIGE

saige_singularity_dir="/well/lindgren/dpalmer/saige_singularity"
version="0.45"
singularity shell --bind /well/lindgren/UKBIOBANK,/well/lindgren/dpalmer ${saige_singularity_dir}/saige.${version}.sif

# Navigate to the extdata folder of the old version of SAIGE (https://github.com/weizhouUMICH/SAIGE) and run this:
cd /well/lindgren/dpalmer/SAIGE_old/extdata

phenofile="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/200k/UKBB_WES200k_filtered_cts_dec2021_phenotypes.tsv.gz"
phenotype="Alanine_aminotransferase_residual"
covars="age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,sex,sequencing.batch"
ID_col="eid"
trait_type="quantitative"

plinkFile="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_1_22_including_rare_saige_input"
outputPrefix="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/${phenotype}"
outputPrefix_varRatio="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_SAIGE_output/200k/${phenotype}_cate"
sparse_GRMFile="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
sparse_GRM_SampleID="/well/lindgren/UKBIOBANK/flassen/projects/KO/wes_ko_ukbb/data/saige/grm/input/211101_long_ukb_wes_200k_sparse_autosomes_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

Rscript step1_fitNULLGLMM.R \
    --plinkFile ${plinkFile} \
    --phenoFile ${phenofile} \
    --phenoCol ${phenotype} \
    --covarColList ${covars} \
    --sampleIDColinphenoFile ${ID_col} \
    --traitType ${trait_type} \
    --invNormalize TRUE \
    --LOCO FALSE \
    --outputPrefix ${outputPrefix} --outputPrefix_varRatio ${outputPrefix_varRatio} \
    --IsOverwriteVarianceRatioFile TRUE \
    --sparseGRMFile ${sparse_GRMFile} \
    --sparseGRMSampleIDFile ${sparse_GRM_SampleID} \
    --isCateVarianceRatio TRUE \
    --cateVarRatioMinMACVecExclude "0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude "1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --useSparseGRMtoFitNULL TRUE 
