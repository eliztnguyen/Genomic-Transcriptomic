#written summer 2018 by ETN for Solomon Lab
#this scipt is written to be interactive and user friendly for non-coders
#this script is meant to break down the Allen Brain Aging TBI dataset
#Can be used for each ROI, OR, for DonorInformation
#it breaks it down by different AD criterias (all vs DSM4 vs NINCDS)
#it further breaks down the dataset into male and female subsets to facilitate working with the data on an exploratory level

###running this with files following the 'AllenBrainAgingTBIgenes.R' script will further subset the data accordingly.

#you will need to install the dplyr package prior to accessing it
#the dplyr package only needs to be installed once
#for reference, the code for installing the package is commented out in the next line
#>>install.packages("dplyr") #uncomment and install if needed

#load library
library(dplyr)

#select expression file with donor information for relevant ROI
#you can choose expression data with ALL participants, or EXCLUDING TBI
expressionDonors <- read.csv(file.choose())

#if you want a summary count of participants in each group
summarize(group_by(expressionDonors, sex, act_demented), count=n())
#summary with TBI broken up....
summarize(group_by(expressionDonors, ever_tbi_w_loc, sex, act_demented), count=n())

#choose folder where you want data to be deposited
setwd(choose.dir())

#filter out data based on sex
m_all <- filter(expressionDonors, sex == "M")
f_all <- filter(expressionDonors, sex == "F")
#arrange data by dementia status
m_all <- arrange(m_all, desc(act_demented))
f_all <- arrange(f_all, desc(act_demented))

#write csv files
write.csv(m_all, file = "m_all.csv", row.names = FALSE)
write.csv(f_all, file = "f_all.csv", row.names = FALSE)

########################################################################

#filter out data based on DSM4 -- ONLY No Dementia and AD
mf_DSM4 <- filter(expressionDonors, dsm_iv_clinical_diagnosis == "No Dementia" |
                    dsm_iv_clinical_diagnosis == "Alzheimer's Disease Type")

#if you want a summary count of participants in each group
summarize(group_by(mf_DSM4, sex, dsm_iv_clinical_diagnosis), count=n())
#summary with TBI broken up....
summarize(group_by(mf_DSM4, ever_tbi_w_loc, sex, dsm_iv_clinical_diagnosis), count=n())


#filter out data based on sex
m_DSM4 <- filter(mf_DSM4, sex=="M")
f_DSM4 <- filter(mf_DSM4, sex=="F")
#arrange data by dementia status
m_DSM4 <- arrange(m_DSM4, desc(dsm_iv_clinical_diagnosis))
f_DSM4 <- arrange(f_DSM4, desc(dsm_iv_clinical_diagnosis))

#write csv files
write.csv(m_DSM4, file = "m_DSM4.csv", row.names = FALSE)
write.csv(f_DSM4, file = "f_DSM4.csv", row.names = FALSE)

########################################################################

#filter out data based on NINCDS-ARDA -- ONLY No Dementia, Probable, Possible
mf_NINCDS <- filter(expressionDonors, nincds_arda_diagnosis == "No Dementia" |
                      nincds_arda_diagnosis == "Probable Alzheimer'S Disease" |
                      nincds_arda_diagnosis == "Possible Alzheimer'S Disease")
#recodes NINCDS "probable" to AD
mf_NINCDS_combined <- mutate(mf_NINCDS, nincds_AD_prep = dplyr::recode(mf_NINCDS$nincds_arda_diagnosis, "Probable Alzheimer'S Disease" = "AD"))
#recodes NINCDS "possible" to AD
mf_NINCDS_combined <- mutate(mf_NINCDS_combined, nincds_AD = dplyr::recode(mf_NINCDS_combined$nincds_AD, "Possible Alzheimer'S Disease" = "AD"))
#cleans up data so mf_NINCDS_combined has the in between step removed
mf_NINCDS_combined <- select(mf_NINCDS_combined, -nincds_AD_prep)

#if you want a summary count of participants in each group
summarize(group_by(mf_NINCDS_combined, sex, nincds_AD), count=n())
#summary with TBI broken up....
summarize(group_by(mf_NINCDS_combined, ever_tbi_w_loc, sex, nincds_AD), count=n())

#filter out data based on sex
m_NINCDS_combined <- filter(mf_NINCDS_combined, sex=="M")
f_NINCDS_combined <- filter(mf_NINCDS_combined, sex=="F")
#arrange data by AD status
m_NINCDS_combined <- arrange(m_NINCDS_combined, nincds_AD)
f_NINCDS_combined <- arrange(f_NINCDS_combined, nincds_AD)

#write csv files
write.csv(m_NINCDS_combined, file = "m_NINCDS_combined.csv", row.names = FALSE)
write.csv(f_NINCDS_combined, file = "f_NINCDS_combined.csv", row.names = FALSE)