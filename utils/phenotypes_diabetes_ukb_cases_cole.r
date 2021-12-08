############################################################################################################
#*# I received this script from Jenny Censin on 07 Feb 2020 (see email). Copied content of that script into this script, and then pasted in the following lines: 
#*# The rest of this script was left as is, unless indicated with a comment from me (Marked as #*#). Added some spacing to improve legibility
# UK Biobank - Curated Diabetes Definitions

# Developed for paper: Udler, M. S., et al. (2018). "Type 2 diabetes genetic loci informed by multi-trait associations point to disease mechanisms and subtypes: A soft clustering analysis." PLoS Medicine 15(9): e1002654.
# By Joanne B. Cole 2018 (jcole@broadinstitute.org)
# These R commands simplified and written-up for sharing Feb 6, 2020

# Based on algorithm in Eastwood, S. V., et al. (2016). "Algorithms for the Capture and Adjudication of Prevalent and Incident Diabetes in UK Biobank." PloS One 11(9): e0162388.
# Notes: 
# Includes cases reported at repeat visits from UK Biobank (instance 0, 1, 2)
# Note that comparison with other phenotypes may need to take this into account (i.e. age, BMI, biomarker levels, etc. may all change between visits)
# Includes field 2443 (touchsreen - any diabetes diagnosed by a doctor) in the pipeline, increases sample size, but also potentially noise.
# start with Europeans (based on genetic clustering on 1KGP reference, ~455K) and unrelateds (if interested)
# ends with 1 control (dm_unlikely) and several case definitions - need to be combined for case-control binary phenotype
# Other considerations not included here: HbA1c levels, ICD9/10codes

############################################################################################################
## Read in Phenotype File
library(data.table)
df = fread("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv", data.table=FALSE, sep = ",")
colnames(df) <- paste0("f.", gsub("-", ".", colnames(df)))

############################################################################################################
# Creating self report (sr) diabetes (any/GDM/T1DM/T2DM) Nurse Interview (ni):
# f.20002 = Non-cancer illness code, self-reported, nurse interview

# 1220 = diabetes
df$dm_any_sr_ni.0 <- 0
for (i in 0:28){ df$dm_any_sr_ni.0[df[[paste ('f.20002.0.', i, sep="")]] == 1220] <- 1 } 
table(df$dm_any_sr_ni.0, useNA="always")

df$dm_any_sr_ni.1 <- 0
for (i in 0:28){ df$dm_any_sr_ni.1[df[[paste ('f.20002.1.', i, sep="")]] == 1220] <- 1 } 
table(df$dm_any_sr_ni.1, useNA="always")

df$dm_any_sr_ni.2 <- 0
for (i in 0:28){ df$dm_any_sr_ni.2[df[[paste ('f.20002.2.', i, sep="")]] == 1220] <- 1 } 
table(df$dm_any_sr_ni.2, useNA="always")

# and include the self-diagnosis TS Field of 2443 
#*# ISSUE: this field only has 4 numerical values, hence == "Yes" is not possible.
#*# Changed == "Yes" to == 1, as per http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=100349.
#*# The -1 and -3 values also get allocated as 0, but don't use this anyway. 
df$dm_any_sr_ts.0 <- 0
df$dm_any_sr_ts.0 = ifelse(df$f.2443.0.0 == 1, 1, 0)
table(df$dm_any_sr_ts.0, useNA="always")

df$dm_any_sr_ts.1 <- 0
df$dm_any_sr_ts.1 = ifelse(df$f.2443.1.0 == 1, 1, 0)
table(df$dm_any_sr_ts.1, useNA="always")

df$dm_any_sr_ts.2 <- 0
df$dm_any_sr_ts.2 = ifelse(df$f.2443.2.0 == 1, 1, 0)
table(df$dm_any_sr_ts.2, useNA="always")

# If ever TS or NI diabetes
df$dm_any_sr_ni_ts = ifelse(
  (df$dm_any_sr_ni.0 ==1 | df$dm_any_sr_ni.1 ==1 | df$dm_any_sr_ni.2 ==1 | 
   df$dm_any_sr_ts.0 ==1 | df$dm_any_sr_ts.1 ==1 | df$dm_any_sr_ts.2 ==1), 1, 0)
table(df$dm_any_sr_ni_ts, useNA="always")

# 1221 = gestational diabetes
df$dm_gdm_sr_ni.0 <- 0
for (i in 0:28){ df$dm_gdm_sr_ni.0[df[[paste ('f.20002.0.', i, sep="")]] == 1221] <- 1 }
table(df$dm_gdm_sr_ni.0, useNA="always")

df$dm_gdm_sr_ni.1 <- 0
for (i in 0:28){ df$dm_gdm_sr_ni.1[df[[paste ('f.20002.1.', i, sep="")]] == 1221] <- 1 }
table(df$dm_gdm_sr_ni.1, useNA="always")

df$dm_gdm_sr_ni.2 <- 0
for (i in 0:28){ df$dm_gdm_sr_ni.2[df[[paste ('f.20002.2.', i, sep="")]] == 1221] <- 1 }
table(df$dm_gdm_sr_ni.2, useNA="always")

df$dm_gdm_sr_ni = ifelse((df$dm_gdm_sr_ni.0 ==1 | df$dm_gdm_sr_ni.1 ==1 | df$dm_gdm_sr_ni.2 ==1), 1, 0)
table(df$dm_gdm_sr_ni, useNA="always")

######################################################
# 1222 = T1D
df$dm_t1dm_sr_ni.0 <- 0
for (i in 0:28) { df$dm_t1dm_sr_ni.0[df[[paste ('f.20002.0.', i, sep="")]] == 1222] <- 1 }
table(df$dm_t1dm_sr_ni.0, useNA="always")

df$dm_t1dm_sr_ni.1 <- 0
for (i in 0:28) { df$dm_t1dm_sr_ni.1[df[[paste ('f.20002.1.', i, sep="")]] == 1222] <- 1 }
table(df$dm_t1dm_sr_ni.1, useNA="always")

df$dm_t1dm_sr_ni.2 <- 0
for (i in 0:28) { df$dm_t1dm_sr_ni.2[df[[paste ('f.20002.2.', i, sep="")]] == 1222] <- 1 }
table(df$dm_t1dm_sr_ni.2, useNA="always")

df$dm_t1dm_sr_ni = ifelse((df$dm_t1dm_sr_ni.0 ==1 | df$dm_t1dm_sr_ni.1 ==1 | df$dm_t1dm_sr_ni.2 ==1), 1, 0)
table(df$dm_t1dm_sr_ni, useNA="always")

######################################################
# 1223 = T2D
df$dm_t2dm_sr_ni.0 <- 0
for (i in 0:28){ df$dm_t2dm_sr_ni.0[df[[paste ('f.20002.0.', i, sep="")]] == 1223] <- 1 }
table(df$dm_t2dm_sr_ni.0, useNA="always")

df$dm_t2dm_sr_ni.1 <- 0
for (i in 0:28){ df$dm_t2dm_sr_ni.1[df[[paste ('f.20002.1.', i, sep="")]] == 1223] <- 1 }
table(df$dm_t2dm_sr_ni.1, useNA="always")

df$dm_t2dm_sr_ni.2 <- 0
for (i in 0:28){ df$dm_t2dm_sr_ni.2[df[[paste ('f.20002.2.', i, sep="")]] == 1223] <- 1 }
table(df$dm_t2dm_sr_ni.2, useNA="always")

df$dm_t2dm_sr_ni = ifelse((df$dm_t2dm_sr_ni.0 ==1 | df$dm_t2dm_sr_ni.1 ==1 | df$dm_t2dm_sr_ni.2 ==1), 1, 0)
table(df$dm_t2dm_sr_ni, useNA="always")

######################################################
# Replacing 10 people with both T1 and T2 diagnoses to non-specific
#*# JB comment: this is 11 in our dataset in February 2020.
# these people now have t1t2dm + dm, but have 0's for t1dm and t2dm specifically
df$t1t2dm_sr_ni <- ifelse((df$dm_t1dm_sr_ni == 1 & df$dm_t2dm_sr_ni == 1), 1, 0)
table(df$t1t2dm_sr_ni, useNA="always")

table(df$dm_any_sr_ni_ts, useNA="always")
df$dm_any_sr_ni_ts <- ifelse(df$t1t2dm_sr_ni==1, 1, df$dm_any_sr_ni_ts)
table(df$dm_any_sr_ni_ts, useNA="always")

table(df$dm_t1dm_sr_ni, useNA="always")
df$dm_t1dm_sr_ni <- ifelse(df$t1t2dm_sr_ni==1, 0, df$dm_t1dm_sr_ni)
table(df$dm_t1dm_sr_ni, useNA="always")

table(df$dm_t2dm_sr_ni, useNA="always")
df$dm_t2dm_sr_ni <- ifelse(df$t1t2dm_sr_ni==1, 0, df$dm_t2dm_sr_ni)
table(df$dm_t2dm_sr_ni, useNA="always")

############################################################################################################
# Creating TS GDM Status
# f.4041 = GDM, touchscreen
######################################################
# Creating TS GDM Status

# Did you only have diabetes during pregnancy?
table(df$f.4041.0.0, useNA="always")

# if they said NA or didn't say yes, then gdm is 0, otherwise 1
#*# ISSUE: this field only has 5 numerical values, hence == "Yes" is not possible. Changed != "Yes" to != 1, 
#*# as per http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=100349. 
#*# The other values also get allocated as 0, but don't use this anyway. 
df$dm_gdm_sr_ts.0 <- ifelse(df$f.4041.0.0 != 1 | is.na(df$f.4041.0.0), 0,1)
table(df$dm_gdm_sr_ts.0, useNA="always")
df$dm_gdm_sr_ts.1 <- ifelse(df$f.4041.1.0 != 1 | is.na(df$f.4041.1.0), 0,1)
table(df$dm_gdm_sr_ts.1, useNA="always")
df$dm_gdm_sr_ts.2 <- ifelse(df$f.4041.2.0 != 1 | is.na(df$f.4041.2.0), 0,1)
table(df$dm_gdm_sr_ts.2, useNA="always")

df$dm_gdm_sr_ts <- ifelse((df$dm_gdm_sr_ts.0 == 1 | df$dm_gdm_sr_ts.0 == 1 | df$dm_gdm_sr_ts.0 == 1), 1, 0)
table(df$dm_gdm_sr_ts, useNA="always")

############################################################################################################
# Diagnosis Age
# f.20009 - NI - Interpolated Age of participant when non-cancer illness (f.20002) first diagnosed
# f.2876 - TS - age diabetes diagnosed 
############################################################################################################

# Pull age of dx across visits for different NI SR diabetes codes
df$age_dm_any_dx_sr.0 <- -1
for (i in 0:28) {df$age_dm_any_dx_sr.0 <- ifelse(is.na(df[[paste ('f.20002.0.', i, sep="")]]) == T | df[[paste ('f.20002.0.', i, sep="")]] != 1220, df$age_dm_any_dx_sr.0, df[[paste ('f.20009.0.', i, sep="")]])}
df$age_dm_any_dx_sr.0 <- ifelse(df$age_dm_any_dx_sr.0<0 | is.na(df$age_dm_any_dx_sr.0)==T , NA, df$age_dm_any_dx_sr.0)
table(df$age_dm_any_dx_sr.0, useNA="always")

df$age_dm_gdm_dx_sr.0 <- -1
for (i in 0:28) {df$age_dm_gdm_dx_sr.0 <- ifelse(is.na(df[[paste ('f.20002.0.', i, sep="")]]) == T | df[[paste ('f.20002.0.', i, sep="")]] != 1221, df$age_dm_gdm_dx_sr.0, df[[paste ('f.20009.0.', i, sep="")]])}
df$age_dm_gdm_dx_sr.0 <- ifelse(df$age_dm_gdm_dx_sr.0<0 | is.na(df$age_dm_gdm_dx_sr.0)==T , NA, df$age_dm_gdm_dx_sr.0)
table(df$age_dm_gdm_dx_sr.0, useNA="always")

df$age_dm_t1dm_dx_sr.0 <- -1
for (i in 0:28) {df$age_dm_t1dm_dx_sr.0 <- ifelse(is.na(df[[paste ('f.20002.0.', i, sep="")]]) == T | df[[paste ('f.20002.0.', i, sep="")]] != 1222, df$age_dm_t1dm_dx_sr.0, df[[paste ('f.20009.0.', i, sep="")]])}
df$age_dm_t1dm_dx_sr.0 <- ifelse(df$age_dm_t1dm_dx_sr.0 <0 | is.na(df$age_dm_t1dm_dx_sr.0)==T , NA, df$age_dm_t1dm_dx_sr.0)
table(df$age_dm_t1dm_dx_sr.0, useNA="always")

df$age_dm_t2dm_dx_sr.0 <- -1
for (i in 0:28) {df$age_dm_t2dm_dx_sr.0 <- ifelse(is.na(df[[paste ('f.20002.0.', i, sep="")]]) == T | df[[paste ('f.20002.0.', i, sep="")]] != 1223, df$age_dm_t2dm_dx_sr.0, df[[paste ('f.20009.0.', i, sep="")]])}
df$age_dm_t2dm_dx_sr.0 <- ifelse(df$age_dm_t2dm_dx_sr.0 <0 | is.na(df$age_dm_t2dm_dx_sr.0)==T , NA, df$age_dm_t2dm_dx_sr.0)
table(df$age_dm_t2dm_dx_sr.0, useNA="always")


df$age_dm_any_dx_sr.1 <- -1
for (i in 0:28) {df$age_dm_any_dx_sr.1 <- ifelse(is.na(df[[paste ('f.20002.1.', i, sep="")]]) == T | df[[paste ('f.20002.1.', i, sep="")]] != 1220, df$age_dm_any_dx_sr.1, df[[paste ('f.20009.1.', i, sep="")]])}
df$age_dm_any_dx_sr.1 <- ifelse(df$age_dm_any_dx_sr.1<0 | is.na(df$age_dm_any_dx_sr.1)==T , NA, df$age_dm_any_dx_sr.1)
table(df$age_dm_any_dx_sr.1, useNA="always")

df$age_dm_gdm_dx_sr.1 <- -1
for (i in 0:28) {df$age_dm_gdm_dx_sr.1 <- ifelse(is.na(df[[paste ('f.20002.1.', i, sep="")]]) == T | df[[paste ('f.20002.1.', i, sep="")]] != 1221, df$age_dm_gdm_dx_sr.1, df[[paste ('f.20009.1.', i, sep="")]])}
df$age_dm_gdm_dx_sr.1 <- ifelse(df$age_dm_gdm_dx_sr.1<0 | is.na(df$age_dm_gdm_dx_sr.1)==T , NA, df$age_dm_gdm_dx_sr.1)
table(df$age_dm_gdm_dx_sr.1, useNA="always")

df$age_dm_t1dm_dx_sr.1 <- -1
for (i in 0:28) {df$age_dm_t1dm_dx_sr.1 <- ifelse(is.na(df[[paste ('f.20002.1.', i, sep="")]]) == T | df[[paste ('f.20002.1.', i, sep="")]] != 1222, df$age_dm_t1dm_dx_sr.1, df[[paste ('f.20009.1.', i, sep="")]])}
df$age_dm_t1dm_dx_sr.1 <- ifelse(df$age_dm_t1dm_dx_sr.1 <0 | is.na(df$age_dm_t1dm_dx_sr.1)==T , NA, df$age_dm_t1dm_dx_sr.1)
table(df$age_dm_t1dm_dx_sr.1, useNA="always")

df$age_dm_t2dm_dx_sr.1 <- -1
for (i in 0:28) {df$age_dm_t2dm_dx_sr.1 <- ifelse(is.na(df[[paste ('f.20002.1.', i, sep="")]]) == T | df[[paste ('f.20002.1.', i, sep="")]] != 1223, df$age_dm_t2dm_dx_sr.1, df[[paste ('f.20009.1.', i, sep="")]])}
df$age_dm_t2dm_dx_sr.1 <- ifelse(df$age_dm_t2dm_dx_sr.1 <0 | is.na(df$age_dm_t2dm_dx_sr.1)==T , NA, df$age_dm_t2dm_dx_sr.1)
table(df$age_dm_t2dm_dx_sr.1, useNA="always")


df$age_dm_any_dx_sr.2 <- -1
for (i in 0:28) {df$age_dm_any_dx_sr.2 <- ifelse(is.na(df[[paste ('f.20002.2.', i, sep="")]]) == T | df[[paste ('f.20002.2.', i, sep="")]] != 1220, df$age_dm_any_dx_sr.2, df[[paste ('f.20009.2.', i, sep="")]])}
df$age_dm_any_dx_sr.2 <- ifelse(df$age_dm_any_dx_sr.2<0 | is.na(df$age_dm_any_dx_sr.2)==T , NA, df$age_dm_any_dx_sr.2)
table(df$age_dm_any_dx_sr.2, useNA="always")

df$age_dm_gdm_dx_sr.2 <- -1
for (i in 0:28) {df$age_dm_gdm_dx_sr.2 <- ifelse(is.na(df[[paste ('f.20002.2.', i, sep="")]]) == T | df[[paste ('f.20002.2.', i, sep="")]] != 1221, df$age_dm_gdm_dx_sr.2, df[[paste ('f.20009.2.', i, sep="")]])}
df$age_dm_gdm_dx_sr.2 <- ifelse(df$age_dm_gdm_dx_sr.2<0 | is.na(df$age_dm_gdm_dx_sr.2)==T , NA, df$age_dm_gdm_dx_sr.2)
table(df$age_dm_gdm_dx_sr.2, useNA="always")

df$age_dm_t1dm_dx_sr.2 <- -1
for (i in 0:28) {df$age_dm_t1dm_dx_sr.2 <- ifelse(is.na(df[[paste ('f.20002.2.', i, sep="")]]) == T | df[[paste ('f.20002.2.', i, sep="")]] != 1222, df$age_dm_t1dm_dx_sr.2, df[[paste ('f.20009.2.', i, sep="")]])}
df$age_dm_t1dm_dx_sr.2 <- ifelse(df$age_dm_t1dm_dx_sr.2 <0 | is.na(df$age_dm_t1dm_dx_sr.2)==T , NA, df$age_dm_t1dm_dx_sr.2)
table(df$age_dm_t1dm_dx_sr.2, useNA="always")

df$age_dm_t2dm_dx_sr.2 <- -1
for (i in 0:28) {df$age_dm_t2dm_dx_sr.2 <- ifelse(is.na(df[[paste ('f.20002.2.', i, sep="")]]) == T | df[[paste ('f.20002.2.', i, sep="")]] != 1223, df$age_dm_t2dm_dx_sr.2, df[[paste ('f.20009.2.', i, sep="")]])}
df$age_dm_t2dm_dx_sr.2 <- ifelse(df$age_dm_t2dm_dx_sr.2 <0 | is.na(df$age_dm_t2dm_dx_sr.2)==T , NA, df$age_dm_t2dm_dx_sr.2)
table(df$age_dm_t2dm_dx_sr.2, useNA="always")

# Average over visits 
df$age_dm_any_dx_sr = rowMeans(df[,c("age_dm_any_dx_sr.0", "age_dm_any_dx_sr.1", "age_dm_any_dx_sr.2")],na.rm = TRUE)
table(df$age_dm_any_dx_sr, useNA="always")
df$age_dm_any_dx_sr[is.nan(df$age_dm_any_dx_sr)] <- NA
table(df$age_dm_any_dx_sr, useNA="always")

df$age_dm_gdm_dx_sr = rowMeans(df[,c("age_dm_gdm_dx_sr.0", "age_dm_gdm_dx_sr.1", "age_dm_gdm_dx_sr.2")],na.rm = TRUE)
table(df$age_dm_gdm_dx_sr, useNA="always")
df$age_dm_gdm_dx_sr[is.nan(df$age_dm_gdm_dx_sr)] <- NA
table(df$age_dm_gdm_dx_sr, useNA="always")

df$age_dm_t1dm_dx_sr = rowMeans(df[,c("age_dm_t1dm_dx_sr.0", "age_dm_t1dm_dx_sr.1", "age_dm_t1dm_dx_sr.2")],na.rm = TRUE)
table(df$age_dm_t1dm_dx_sr, useNA="always")
df$age_dm_t1dm_dx_sr[is.nan(df$age_dm_t1dm_dx_sr)] <- NA
table(df$age_dm_t1dm_dx_sr, useNA="always")

df$age_dm_t2dm_dx_sr = rowMeans(df[,c("age_dm_t2dm_dx_sr.0", "age_dm_t2dm_dx_sr.1", "age_dm_t2dm_dx_sr.2")],na.rm = TRUE)
table(df$age_dm_t2dm_dx_sr, useNA="always")
df$age_dm_t2dm_dx_sr[is.nan(df$age_dm_t2dm_dx_sr)] <- NA
table(df$age_dm_t2dm_dx_sr, useNA="always")

# f.2976 = age diabetes diagnosed - TS
# set -1 and -3 to NA
df$f.2976.0.0 = ifelse(df$f.2976.0.0 < 0, NA, df$f.2976.0.0)
table(df$f.2976.0.0, useNA="always")
df$f.2976.1.0 = ifelse(df$f.2976.1.0 < 0, NA, df$f.2976.1.0)
table(df$f.2976.1.0, useNA="always")
df$f.2976.2.0 = ifelse(df$f.2976.2.0 < 0, NA, df$f.2976.2.0)
table(df$f.2976.2.0, useNA="always")

# average TS dx age
df$f.2976 = rowMeans(df[,c("f.2976.0.0", "f.2976.1.0", "f.2976.2.0")],na.rm = TRUE)
df$f.2976[is.nan(df$f.2976)] <- NA
table(df$f.2976, useNA="always")

# if you have T1DM -> age at T1DM
# if you have T2DM -> age at T2DM
# if age at T1DM, T2DM, and any DM is NA, AND you do not have gestational diabetes -> TS age at dx (f.2976)
# this is because "Field 2976 was collected from men who indicated that a doctor had told them 
# they have diabetes, as defined by their answers to Field 2443 and all women except those who 
# indicated they had diabetes only during pregnancy, as defined by their answers to Field 2443
# order of operations: T2DM age, T1DM age, then if no T1 or T2 or GDM, then 2976, then any dm age

df$agedm_ts_or_ni <- ifelse(df$dm_t2dm_sr_ni == 1, df$age_dm_t2dm_dx_sr, ifelse(df$dm_t1dm_sr_ni == 1, df$age_dm_t1dm_dx_sr, ifelse(is.na(df$age_dm_any_dx_sr)==T & is.na(df$age_dm_t1dm_dx_sr)==T & is.na(df$age_dm_t2dm_dx_sr)==T & df$dm_gdm_sr_ts==0, df$f.2976, df$age_dm_any_dx_sr)))
table(df$agedm_ts_or_ni, useNA="always")
df$agedm_ts_or_ni <- ifelse(is.na(df$f.2976) ==T & is.na(df$age_dm_any_dx_sr)==T & is.na(df$age_dm_t1dm_dx_sr)==T & is.na(df$age_dm_t2dm_dx_sr)==T, NA, df$agedm_ts_or_ni)
table(df$agedm_ts_or_ni, useNA="always")
df$agedm_ts_or_ni <- ifelse(df$agedm_ts_or_ni < 0, NA, df$agedm_ts_or_ni)
table(df$agedm_ts_or_ni, useNA="always")

############################################################################################################
# Insulin & Medication Use
# f.6177 and f.6153 for males and females - medication use - TS
# f.20003 - SR NI medication use 
############################################################################################################

# Insulin - SR TS
#*# ISSUE: fields 6177 and 6153 only have numerical values, hence == "Insulin" is not possible. 
#*# Changed == "Insulin" to = 3 for both fields, as per 
#*# http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=100625 and http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=100626.

df$dm_insulin_sr_ts <- 0 
df$dm_insulin_sr_ts.0 <- 0
for (i in 0:2){df$dm_insulin_sr_ts.0[df[[paste ('f.6177.0.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.0, useNA="always")
for (i in 0:3){df$dm_insulin_sr_ts.0[df[[paste ('f.6153.0.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.0, useNA="always")

df$dm_insulin_sr_ts.1 <- 0
for (i in 0:2){df$dm_insulin_sr_ts.1[df[[paste ('f.6177.1.', i, sep="")]] == 3] <- 1 }
table(df$dm_insulin_sr_ts.1, useNA="always")
for (i in 0:3){df$dm_insulin_sr_ts.1[df[[paste ('f.6153.1.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.1, useNA="always")

df$dm_insulin_sr_ts.2 <- 0
for (i in 0:2){df$dm_insulin_sr_ts.2[df[[paste ('f.6177.2.', i, sep="")]] == 3] <- 1 }
table(df$dm_insulin_sr_ts.2, useNA="always")
for (i in 0:3){df$dm_insulin_sr_ts.2[df[[paste ('f.6153.2.', i, sep="")]] == 3] <- 1 } 
table(df$dm_insulin_sr_ts.2, useNA="always")

df$dm_insulin_sr_ts = ifelse((df$dm_insulin_sr_ts.0 == 1 | df$dm_insulin_sr_ts.1 == 1 | df$dm_insulin_sr_ts.2 == 1), 1, 0)
table(df$dm_insulin_sr_ts, useNA="always")

# Other medication: insulin, metformin, non-metformin OAD
df$meds_insulin_sr_ni.0 <- 0
for (i in 0:47){ df$meds_insulin_sr_ni.0[df[[paste ('f.20003.0.', i, sep="")]] == 1140883066] <- 1 }
df$meds_metformin_sr_ni.0 <- 0
table(df$meds_insulin_sr_ni.0, useNA="always")

v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47){ df$meds_metformin_sr_ni.0[df[[paste ('f.20003.0.', i, sep="")]] %in% v] <- 1 }
table(df$meds_metformin_sr_ni.0, useNA="always")

df$meds_nonmet_oad_sr_ni.0 <- 0
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47){ df$meds_nonmet_oad_sr_ni.0[df[[paste ('f.20003.0.', i, sep="")]] %in% v] <- 1 }
table(df$meds_nonmet_oad_sr_ni.0, useNA="always")

df$meds_insulin_sr_ni.1 <- 0
for (i in 0:47){ df$meds_insulin_sr_ni.1[df[[paste ('f.20003.1.', i, sep="")]] == 1140883066] <- 1 }
df$meds_metformin_sr_ni.1 <- 0
table(df$meds_insulin_sr_ni.1, useNA="always")

v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47){ df$meds_metformin_sr_ni.1[df[[paste ('f.20003.1.', i, sep="")]] %in% v] <- 1 }
table(df$meds_metformin_sr_ni.1, useNA="always")

df$meds_nonmet_oad_sr_ni.1 <- 0
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47){ df$meds_nonmet_oad_sr_ni.1[df[[paste ('f.20003.1.', i, sep="")]] %in% v] <- 1 }
table(df$meds_nonmet_oad_sr_ni.1, useNA="always")

df$meds_insulin_sr_ni.2 <- 0
for (i in 0:47){ df$meds_insulin_sr_ni.2[df[[paste ('f.20003.2.', i, sep="")]] == 1140883066] <- 1 }

df$meds_metformin_sr_ni.2 <- 0
table(df$meds_insulin_sr_ni.2, useNA="always")

v <- c(1140884600, 1140874686, 1141189090)
for (i in 0:47){ df$meds_metformin_sr_ni.2[df[[paste ('f.20003.2.', i, sep="")]] %in% v] <- 1 }
table(df$meds_metformin_sr_ni.2, useNA="always")

df$meds_nonmet_oad_sr_ni.2 <- 0
v <- c(1140874718, 1140874744, 1140874746, 1141152590, 1141156984, 1140874646, 1141157284, 1140874652, 1140874674, 1140874728, 1140868902, 1140868908, 1140857508, 1141173882, 1141173786, 1141168660, 1141171646, 1141171652, 1141153254, 1141177600, 1141177606)
for (i in 0:47){ df$meds_nonmet_oad_sr_ni.2[df[[paste ('f.20003.2.', i, sep="")]] %in% v] <- 1 }
table(df$meds_nonmet_oad_sr_ni.2, useNA="always")

df$meds_insulin_sr_ni = ifelse((df$meds_insulin_sr_ni.0 == 1 | df$meds_insulin_sr_ni.1 == 1 | df$meds_insulin_sr_ni.2 == 1), 1, 0)
table(df$meds_insulin_sr_ni, useNA="always")
df$meds_metformin_sr_ni = ifelse((df$meds_metformin_sr_ni.0 == 1 | df$meds_metformin_sr_ni.1 == 1 | df$meds_metformin_sr_ni.2 == 1), 1, 0)
table(df$meds_metformin_sr_ni, useNA="always")
df$meds_nonmet_oad_sr_ni = ifelse((df$meds_nonmet_oad_sr_ni.0 == 1 | df$meds_nonmet_oad_sr_ni.1 == 1 | df$meds_nonmet_oad_sr_ni.2 == 1), 1, 0)
table(df$meds_nonmet_oad_sr_ni, useNA="always")

# Combination medication variables
df$meds_any_sr_ni <- 0 # this single variable captures all 3 categories of DM medications
df$meds_any_sr_ni[(df$meds_insulin_sr_ni ==1 | df$meds_metformin_sr_ni ==1 | df$meds_nonmet_oad_sr_ni ==1)] <- 1
table(df$meds_any_sr_ni, useNA="always")
table(df$meds_metformin_sr_ni, useNA="always")
table(df$meds_nonmet_oad_sr_ni, useNA="always") 

df$meds_any_sr_ni_ts <- 0 # this single variable captures all 3 categories of DM medications AND the insulin TS
df$meds_any_sr_ni_ts[(df$meds_insulin_sr_ni ==1 | df$meds_metformin_sr_ni ==1 | df$meds_nonmet_oad_sr_ni ==1 | df$dm_insulin_sr_ts ==1)] <- 1
table(df$dm_insulin_sr_ts, useNA="always")

table(df$meds_any_sr_ni, useNA="always")
table(df$meds_any_sr_ni_ts, useNA="always")

############################################################################################################
# PLOS ONE Flowchart to create cases and controls 
# Flowchart... 

# 1.1
df$dm_unlikely <- 1
# any SR diabetes or T1DM or T2DM from NI -> they have 0 for diabetes unlikely
df$dm_unlikely[(df$dm_any_sr_ni_ts == 1 | df$dm_t1dm_sr_ni == 1 | df$dm_t2dm_sr_ni == 1)] <- 0
# any on diabetes medications -> they have 0 for diabetes unlikely
df$dm_unlikely[(df$meds_any_sr_ni == 1 | df$dm_insulin_sr_ts == 1)] <- 0
# any female with gestational diabetes -> they have 0 for diabetes unlikely
df$dm_unlikely <- ifelse((df$f.31.0.0==0 & (df$dm_gdm_sr_ni==1 | df$dm_gdm_sr_ts==1)), 0, df$dm_unlikely)
# create a field if diabetes is NOT unlikely -> pass on step 1.2
df$tostep1.2 <- ifelse(df$dm_unlikely == 1, 0, 1)

# 1.2
df$possible_gdm <- 0
# if they are NOT unlikely for diabetes, and they have gestational diabets (TS) and no medication and no insulin and no T1D or T2D diagnosis -> then they have possible gestational diabetes only
df$possible_gdm[(df$tostep1.2==1 & df$dm_gdm_sr_ts == 1 & df$meds_any_sr_ni == 0 & df$dm_insulin_sr_ts == 0 & df$dm_t1dm_sr_ni == 0 & df$dm_t2dm_sr_ni == 0)] <- 1
# if they are NOT unlikely for diabetes, and they have gestational diabetes (NI) and diagnosis of any diabetes is < 50, and they have no medication and no insulin and and no T1D or T2D diagnosis -> then they have possible gestational diabetes only
df$possible_gdm[(df$tostep1.2==1 & df$dm_gdm_sr_ni == 1 & df$age_dm_gdm_dx_sr < 50 & df$meds_any_sr_ni == 0 & df$dm_insulin_sr_ts == 0 & df$dm_t1dm_sr_ni == 0 & df$dm_t2dm_sr_ni == 0)] <- 1
df$tostep1.3 <- 0
# create a field if diabetes is NOT unlikely and possible gestational diabetes is 0 -> pass on to step 1.3
# 1.3 = likely diabetes, NOT gestational 
df$tostep1.3[(df$tostep1.2 == 1 & df$possible_gdm == 0)] <- 1

# 1.3
df$possible_t2dm_a <- 0
# if they might have diabetes and might not have gestational diabets, AND non-metformin OAD
df$possible_t2dm_a[(df$tostep1.3==1 & df$meds_nonmet_oad_sr_ni == 1)] <- 1
df$tostep1.4 <- 0
# 1.4 = likely diabetes, not gestational, and NOT non-metformin OAD 
df$tostep1.4[(df$tostep1.3 == 1 & df$possible_t2dm_a == 0)] <- 1

# 1.4
# if they are european (all in my dataset) or other and age of diagnosis is between 0 and 36. This equals 1. 
df$europeans_36 <- ifelse((df$agedm_ts_or_ni > 0 & df$agedm_ts_or_ni < 37), 1, 0)

df$tostep1.5 <- 0
# likely diabetes, not gestational, NOT non-metformim OAD, EUR and diagnosed 0-36
df$tostep1.5[(df$tostep1.4==1 & df$europeans_36 == 1)] <- 1
df$possible_t2dm_b <- 0
# likely diabetes, not gestational, NOT non-metformin OAD, EUR, diagnosed ABOVE 36 
df$possible_t2dm_b[(df$tostep1.4 ==1 & df$tostep1.5 == 0)] <- 1

#1.5
df$possible_t1dm_temp <- 0
# likely diabetes, young dx, and insulin and/or insulin within 1year of diagnosis
#*# ISSUE: field 2986 only has numerical values, hence == "Yes" is not possible. 
#*# Changed == "Yes" to == 1, as per http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=100349 

df$possible_t1dm_temp[(df$tostep1.5 ==1 & df$dm_insulin_sr_ts == 1 & df$meds_insulin_sr_ni ==1)] <-1
df$possible_t1dm_temp[(df$tostep1.5 ==1 & (df$f.2986.0.0 == 1 | df$f.2986.1.0 == 1 | df$f.2986.2.0 == 1))] <- 1
# likely diabetes, young dx, self reported T1D
df$possible_t1dm_temp[(df$tostep1.5 ==1 & df$dm_t1dm_sr_ni ==1)] <- 1
table(df$possible_t1dm_temp, useNA="always")
df$possible_t2dm_c <- 0
# likely diabetes, young dx ages, and NOT self-reported T1D
df$possible_t2dm_c[(df$tostep1.5 ==1 & df$possible_t1dm_temp == 0)] <- 1
df$possible_t2dm_temp <- 0
df$possible_t2dm_temp[(df$possible_t2dm_a == 1 | df$possible_t2dm_b == 1 | df$possible_t2dm_c == 1)] <- 1

###################################################### 
# Summary
# dm_unlikely = SR, diabetes medication, or gestational only (AKA this field ==0 is "likely diabetes" below)
# possible_gdm = likely diabetes, gestational diabetes, diagnosis < 50, no medication, no insulin, no T1D or T2D diagnosis
# 1.3 = likely diabetes, not gestational
# possible_t2dm_a = likely diabetes (notgbm), and on non-metformin OAD medication
# 1.4 = likely diabetes (notgbm) - all dx ages - and NOT on non-metformin OAD medication
# 1.5 = likely diabetes (notgbm) - young dx ages - and NOT on non-metformin OAD medication
# possible_t2dm_b = likely diabetes (notgbm) - old dx ages
# possible_t1dm_temp = likely diabetes (notgbm) - young dx ages - and NOT on non-metformin OAD medication, and insulin and/or insulin within 1year of diagnosis, and self reported T1D
# possible_t2dm_c = likely diabetes (notgbm) - young dx ages - and NOT on non-metformin OAD medication, and NOT self-reported T1D
# possible_t2dm_temp = any of the following:
# likely diabetes (notgbm), and on non-metformin OAD medication OR 
# likely diabetes (notgbm) - old dx ages OR
# likely diabetes (notgbm) - young dx ages - and NOT on non-metformin OAD medication, and NOT self-reported T1D

# Figure 1b. Prevalence algorithm 2: Adjudicating type 1 diabetes and classification into probable and possible categories. 

# 2.1
df$probable_t1dm_a <- 0
df$probable_t1dm_a[df$possible_t1dm_temp == 1 & df$dm_t1dm_sr_ni == 1] <- 1
df$tostep2.2 <- 0
df$tostep2.2[df$possible_t1dm_temp ==1 & df$probable_t1dm_a == 0] <- 1

table(df$possible_t1dm_temp, useNA="always")
table(df$probable_t1dm_a, useNA="always")
table(df$tostep2.2, useNA="always")

# 2.2
df$probable_t1dm_b <- 0
df$probable_t1dm_b[(df$tostep2.2 ==1 & df$dm_insulin_sr_ts == 1 & df$f.2986.0.0 == 1)] <- 1
df$probable_t1dm_b[(df$tostep2.2 ==1 & df$meds_insulin_sr_ni == 1 & df$f.2986.0.0 == 1)] <- 1

df$probable_t1dm <- 0
df$probable_t1dm[df$probable_t1dm_a == 1 | df$probable_t1dm_b ==1] <- 1

df$possible_t1dm <- 0
df$possible_t1dm[df$tostep2.2 ==1 & df$probable_t1dm_b == 0] <- 1
table(df$possible_t1dm, useNA="always")

# probable_t1dm_a = possible_t1dm_temp = likely diabetes (notgbm) - young dx ages - and NOT on non-metformin OAD medication, and insulin and/or insulin within 1year of diagnosis, and self reported T1D
# 2.2 = likely diabetes (notgbm), young dx ages, NOT on non-metformin OAD medication, and insulin and/or insulin within 1year of diagnosis, and self reported T1D
# probable_t1dm_b = 2.2 + insulin
# TF: probable_t1dm = a or b
# TF: possible_t1dm = small subset
# possible_t1dm = a or b:

# Figure 1c. Prevalence algorithm 3: Adjudicating type 2 diabetes and classification into probable and possible categories. 

# 3.1
df$tostep3.2 <- 0
df$tostep3.2[(df$possible_t2dm_temp == 1 & df$meds_metformin_sr_ni == 1 & df$meds_insulin_sr_ni ==0 & df$meds_nonmet_oad_sr_ni ==0 & df$dm_insulin_sr_ts ==0)] <- 1

df$tostep3.3_a <- 0
df$tostep3.3_a[(df$possible_t2dm_temp == 1 & df$tostep3.2 == 0)] <- 1

# 3.2
df$dm_unlikely_b <- 0
df$dm_unlikely_b[(df$tostep3.2 ==1 & df$dm_any_sr_ni_ts == 0 & df$dm_gdm_sr_ni == 0 & df$dm_t1dm_sr_ni == 0 & df$dm_t2dm_sr_ni == 0)] <- 1
df$dm_unlikely[df$dm_unlikely_b==1] <-1

df$tostep3.3_b <- 0
df$tostep3.3_b[(df$tostep3.2 ==1 & df$dm_unlikely_b == 0)] <- 1

df$tostep3.3 <- 0
df$tostep3.3[(df$tostep3.3_a ==1 | df$tostep3.3_b ==1)] <- 1

# 3.3
df$probable_t2dm_a <- 0
df$probable_t2dm_a[(df$tostep3.3 == 1 & df$meds_nonmet_oad_sr_ni ==1)] <- 1
df$tostep3.4 <- 0
df$tostep3.4[(df$tostep3.3 == 1 & df$probable_t2dm_a ==0)] <- 1

# 3.4
df$probable_t2dm_b <- 0
df$probable_t2dm_b[(df$tostep3.4 == 1 & df$dm_insulin_sr_ts == 0 & df$meds_insulin_sr_ni == 0)] <- 1

df$probable_t2dm <- 0
df$probable_t2dm[(df$probable_t2dm_a == 1 | df$probable_t2dm_b ==1)] <- 1

df$tostep3.5 <- 0
df$tostep3.5[(df$tostep3.4 == 1 & df$probable_t2dm == 0)] <- 1

# 3.5
df$possible_t2dm <- 0
df$possible_t2dm[(df$tostep3.5 ==1 & df$dm_t1dm_sr_ni == 0)] <- 1

df$probable_t1dm_c <- 0
df$probable_t1dm_c[(df$tostep3.5 ==1 & df$dm_t1dm_sr_ni == 1)] <- 1
df$probable_t1dm[df$probable_t1dm_c ==1] <- 1

# Create an age discrepancy flag field between visits and TS and NI
# This is greater than 10 years
df$age_discrepancy = NULL
df$age_discrepancy = ifelse(
  (
    abs(df$age_dm_any_dx_sr.0 - df$age_dm_any_dx_sr.1) > 10 |
    abs(df$age_dm_any_dx_sr.0 - df$age_dm_any_dx_sr.2) > 10 |
    abs(df$age_dm_any_dx_sr.1 - df$age_dm_any_dx_sr.2) > 10 |
    abs(df$age_dm_gdm_dx_sr.0 - df$age_dm_gdm_dx_sr.1) > 10 |
    abs(df$age_dm_gdm_dx_sr.0 - df$age_dm_gdm_dx_sr.2) > 10 |
    abs(df$age_dm_gdm_dx_sr.1 - df$age_dm_gdm_dx_sr.2) > 10 |
    abs(df$age_dm_t2dm_dx_sr.0 - df$age_dm_t2dm_dx_sr.1) > 10 |
    abs(df$age_dm_t2dm_dx_sr.0 - df$age_dm_t2dm_dx_sr.2) > 10 |
    abs(df$age_dm_t2dm_dx_sr.1 - df$age_dm_t2dm_dx_sr.2) > 10 |
    abs(df$age_dm_t1dm_dx_sr.0 - df$age_dm_t1dm_dx_sr.1) > 10 |
    abs(df$age_dm_t1dm_dx_sr.0 - df$age_dm_t1dm_dx_sr.2) > 10 |
    abs(df$age_dm_t1dm_dx_sr.1 - df$age_dm_t1dm_dx_sr.2) > 10 |
    abs(df$f.2976.0.0 - df$f.2976.1.0) > 10 |
    abs(df$f.2976.0.0 - df$f.2976.2.0) > 10 |
    abs(df$f.2976.1.0 - df$f.2976.2.0) > 10 |
    abs(df$age_dm_any_dx_sr - df$f.2976) > 10
    ), 1, NA)
df$age_discrepancy[is.na(df$age_discrepancy)] <- 0
table(df$age_discrepancy, useNA="always")

############################################################################################################
# Simplify
#*# JB comment: change the ID column to "f.eid".

UKB_diabetes <- subset(df, select=c(
  "f.eid",
  "age_discrepancy",
  "agedm_ts_or_ni", # combined field for age of diagnosis: agedm_ts_or_ni
  "dm_insulin_sr_ts",
  "meds_insulin_sr_ni",
  "meds_metformin_sr_ni",
  "meds_nonmet_oad_sr_ni",
  "meds_any_sr_ni", # insulin nurse interview, metformin, nonmet_oad
  "meds_any_sr_ni_ts", # insulin nurse interview, metformin, nonmet_oad, AND insulin touchscreen
  'dm_unlikely', 'possible_t1dm', 'probable_t1dm', 'possible_t2dm', 'probable_t2dm', 'possible_gdm'
))

#*# JB: added lines to save the final extracted table: 
file <- "/well/lindgren/UKBIOBANK/dpalmer/ukbb_diabetes_EastwooddmAlgorithm.tsv"
fwrite(UKB_diabetes, file, sep="\t")
