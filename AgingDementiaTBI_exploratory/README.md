# Aging Dementia TBI - exploratory
These scripts were written to be interactive and user friendly for non-coders.
This script is meant to organize RNAseq data downloaded for any one gene from the Allen Brain Aging and TBI Database.
https://aging.brain-map.org/overview/explore

Scripts:
AllenBrainAgingTBIgenes.R = Associates RNAseq data with the relevant donor ID/info and breaks the data by region of interest (ROI);
AllenBrain_ROIexpression_byGroup.R = Takes the ROI data and subsets it by Alzheimer's status and sex.

Files for sample run/reference:
Contents.txt = summary of data downloaded from database;
Expression.csv = RNA-sequencing expression data for the NR3C1 gene (glucocorticoid receptor);
Columns.csv = donor ID and sample information associated with the RNA-sequencing expression data for the relevant gene;
DonorInformation.csv = demographic and clinical information associated with each donor ID;




