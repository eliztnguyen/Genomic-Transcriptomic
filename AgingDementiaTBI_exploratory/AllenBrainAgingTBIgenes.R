#written summer 2018 by ETN for Solomon Lab -- with initial guidance/encouragement from her friend Andrew Aycoth
#this script is written to be interactive and user friendly for non-coders
#This script is meant to organize RNAseq data downloaded for any one gene from the
#Allen Brain Aging and TBI Database
#This script will work best if all data for the given gene is collected into one folder/directory
#Naming the folder with the gene name will probably be helpful
#Set the working directory (setwd) to said folder when prompted
#Then you can easily access the requested files (expression, columns) for said gene
#Donor information is downloaded independent from the gene data
#An expression datasheet for each ROI (Hipp, PCx, TCx, FWM) will be written separately to this folder

#you will need to install the dplyr package prior to accessing it
#the dplyr package only needs to be installed once
#for reference, the code for installing the package is commented out in the next line
#install.packages("dplyr") #uncomment and install if needed

#load library
library(dplyr)

#choose folder where your gene data files are
setwd(choose.dir())

#read relevant files into R
expression <- read.csv(file.choose(),head=FALSE)
#expression file has no header, so we need to note this when reading it in
columns <- read.csv(file.choose())
donors <- read.csv(file.choose())

#make gene ID name of row **found to be NOT necessary**
#rownames(expression) <- expression[,1]

#remove gene name from expression to enable joining
subset_expression <- expression[-1]
#transpose data and form into dataframe (table) rather than a matrix
subset_expression_trans <- as.data.frame(t(subset_expression))

#add expression values as a column(c) to the donor details
columnsE <- cbind(subset_expression_trans, columns)
#merge donor list and expression donor details into one data frame
donorColumnsE <- merge(donors, columnsE, by="donor_id")
#select out only donor info and relevant ROI and expression info
min_donorColumnsE <- select(donorColumnsE, V1, structure_abbreviation, donor_id:nia_reagan)

#filter out data based on ROI
HIP <- filter(min_donorColumnsE, structure_abbreviation == "HIP")
PCx <- filter(min_donorColumnsE, structure_abbreviation == "PCx")
TCx <- filter(min_donorColumnsE, structure_abbreviation == "TCx")
FWM <- filter(min_donorColumnsE, structure_abbreviation == "FWM")

#choose folder where you want data for all participants to be deposited
setwd(choose.dir())

#write data to csv files based on ROI
write.csv(HIP, file = "HIP.csv", row.names = FALSE)
write.csv(PCx, file = "PCx.csv", row.names = FALSE)
write.csv(TCx, file = "TCx.csv", row.names = FALSE)
write.csv(FWM, file = "FWM.csv", row.names = FALSE)

#filter out data for only participants WITHOUT TBI history
HIP_noTBI <- filter(HIP, ever_tbi_w_loc == "N")
PCx_noTBI <- filter(PCx, ever_tbi_w_loc == "N")
TCx_noTBI <- filter(TCx, ever_tbi_w_loc == "N")
FWM_noTBI <- filter(FWM, ever_tbi_w_loc == "N")

#choose folder where you want data for noTBI participants to be deposited
setwd(choose.dir())

#write data to csv files based on ROI excluding those with TBI history
write.csv(HIP_noTBI, file = "HIP_noTBI.csv", row.names = FALSE)
write.csv(PCx_noTBI, file = "PCx_noTBI.csv", row.names = FALSE)
write.csv(TCx_noTBI, file = "TCx_noTBI.csv", row.names = FALSE)
write.csv(FWM_noTBI, file = "FWM_noTBI.csv", row.names = FALSE)
