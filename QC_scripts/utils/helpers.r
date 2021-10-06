library(data.table)
library(dplyr)

create_pheno_dt <- function(TRANCHE)
{
    PHENOFILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/", TRANCHE, '/QC_phenotypes.tsv.gz')
    
    dt <- fread(cmd=paste('zcat', PHENOFILE)) %>% 
        mutate(s=ID) %>% 
        select(
            s, ukbb.centre, Submitted.Gender, Inferred.Gender, genotyping.array, sequencing.batch,
            self.report.ethnicity, genetic.eur.oct2021, genetic.eur.no.fin.oct2021,
            genetically_european, `in.white.British.ancestry.subset`)
    
    # Rename to obtain the non-Europeans
    dt <- dt %>% mutate(
        white.british = ifelse(`in.white.British.ancestry.subset` == 1, "White-British", "Non white-British"),
        genetic.eur.no.fin.oct2021 = ifelse(genetic.eur.no.fin.oct2021, "NFE", "Non-NFE"),
        genetic.eur.oct2021 = ifelse(genetic.eur.oct2021, "European", "Non-European"),
        sequencing.batch = ifelse(sequencing.batch == 1, "Batch 1", "Batch 2"),
        ukbb.centre = factor(ukbb.centre),
        self.report.ethnicity = factor(self.report.ethnicity)
        )

    # Re-level manually.
    dt$self.report.ethnicity <- recode_factor(
        dt$self.report.ethnicity,
        `1` = "White",
        `1001` = "British",
        `2001` = "White and Black Caribbean",
        `3001` = "Indian",
        `4001` = "Caribbean",
        `2` = " Mixed",
        `1002` = "Irish",
        `2002` = "White and Black African",
        `3002` = "Pakistani",
        `4002` = "African",
        `3` = " Asian or Asian British",
        `1003` = "Any other white background",
        `2003` = "White and Asian",
        `3003` = "Bangladeshi",
        `4003` = "Any other Black background",
        `4` = " Black or Black British",
        `2004` = "Any other mixed background",
        `3004` = "Any other Asian background",
        `5` = "Chinese",
        `6` = "Other ethnic group",
        `-1` = "Do not know",
        `-3` = "Prefer not to answer"
        )

    dt$ukbb.centre <- recode_factor(
        dt$ukbb.centre,
        `11012` =  "Barts",
        `11021` =  "Birmingham",
        `11011` =  "Bristol",
        `11008` =  "Bury",
        `11003` =  "Cardiff",
        `11024` =  "Cheadle (revisit)",
        `11020` =  "Croydon",
        `11005` =  "Edinburgh",
        `11004` =  "Glasgow",
        `11018` =  "Hounslow",
        `11010` =  "Leeds",
        `11016` =  "Liverpool",
        `11001` =  "Manchester",
        `11017` =  "Middlesborough",
        `11009` =  "Newcastle",
        `11013` =  "Nottingham",
        `11002` =  "Oxford",
        `11007` =  "Reading",
        `11014` =  "Sheffield",
        `10003` =  "Stockport (pilot)",
        `11006` =  "Stoke",
        `11022` =  "Swansea",
        `11023` =  "Wrexham",
        `11025` =  "Cheadle (imaging)",
        `11026` =  "Reading (imaging)",
        `11027` =  "Newcastle (imaging)",
        `11028` =  "Bristol (imaging)"
        )
    
    dt <- data.table(dt)
    setkey(dt, 's')
    return(dt)
}