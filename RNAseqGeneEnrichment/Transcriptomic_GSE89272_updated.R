### GERO 599
### Elizabeth Nguyen
### Transcriptomics assignment for GSE89272

###################################################################################################################
# Dataset for Assigment: Diverse interventions that extend mouse lifespan
# Publicly available from GEO Datasets repository
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89272
### Series GSE89272 (Public on Feb 16, 2017)
### Diverse interventions that extend mouse lifespan suppress shared age-associated epigenetic changes at critical gene regulatory regions (RNA-Seq)
### Platform: GPL19057 	Illumina NextSeq 500 (Mus musculus)
###################################################################################################################

# sessionInfo()
#R version 3.6.2 (2019-12-12)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 18362)

# Load required libraries
library('DESeq2')              # DESeq2_1.26.0
library('pvclust')             # pvclust_2.2-0
library('clusterProfiler')     # clusterProfiler_3.14.3
library('org.Mm.eg.db')        # org.Mm.eg.db_3.10.0
library('pheatmap')            # pheatmap_1.0.12
library('Vennerable')          # Vennerable_3.1.0.9000


###################################################################################################################


# Set the working directory
### for Elizabeth's computer
setwd('/Users/nguyenet/Documents/Genomics for Biologists/Transcriptomic Assignment/')
options(stringsAsFactors = F)


###################################################################################################################
### 1. Import and clean RNA-seq data in preparation for DESeq2 analysis
###################################################################################################################

# read in data file as data frame and assign to variable
Hepa.Data <- read.delim('ME_2018-03-26_GSE89272_Hepatocytes_kallisto_mapping.txt')
### view data frame
#View(Hepa.Data)
### currently have gene symbols in first column, followed by Ensemble GeneID, then Ensemble Transcript ID


# clean column names to make more readable
colnames(Hepa.Data) <- sub(".kallisto_res.abundance.tsv", "", colnames(Hepa.Data))   # removes '.kallisto_res.abundance.tsv' from column names
colnames(Hepa.Data) <- sub("Hepatocytesâ..", "", colnames(Hepa.Data))   # removes 'Hepatocytesâ..' from column names


# Get total read count per gene
### file is from kallisto, which maps the transcriptome, not the genome (row/line = one transcript)
### sum all the transcripts of each gene to get a total read count per gene 
Hepa.Data.totalCount <- aggregate(Hepa.Data[,-c(1:6)],
                                 by = list(c(Hepa.Data$GeneSymbol)),
                                 FUN = 'sum')


# Make GeneSymbol row names for DESeq2 to work with later as 'countData'
### rename GeneSymbol (first column of Hepa.Data.totalCount) as row name
rownames(Hepa.Data.totalCount) <- Hepa.Data.totalCount[,1]


# Get read count into integer form since DESeq2 expects count data in the form of a matrix of integer values
### rounds the count per gene to the closest integer
### command below also removes the first column
integer.Hepa.Data.totalCount <- round(Hepa.Data.totalCount[,-1])


# filter low coverage genes
### Judgement call = want to filter for genes expressed in all samples. The actual cutoff depends on sequencing depth and a number of replicates.
### Berenice recommends retaining 10,000-15,000 genes, not more. Less would be too stringent.
### for each row, if expression across all samples is not greater than 50, row removed
filter.integer.Hepa.Data.totalCount <- integer.Hepa.Data.totalCount[rowSums(integer.Hepa.Data.totalCount) >= 50,]


# Write cleaned-up & filtered integer total count data to external file
my.output <- paste(Sys.Date(),"GSE89272_filter.integer.Hepa.Data.totalCount.txt",sep="_")
write.table(filter.integer.Hepa.Data.totalCount, file = my.output, sep = "\t", eol="\r", na = "NA", row.names = T)


# Subset 'filtered.integer.Hepa.Data.totalCount' to include 2 months OR 22 months ONLY
### select for only '2_months' data
Hepa.Data.2Months <- filter.integer.Hepa.Data.totalCount[, 1:8]
### select for only '22_months' data
Hepa.Data.22Months <- filter.integer.Hepa.Data.totalCount[, 9:16]


# Write cleaned-up & filtered integer total count data for 2 months AND 22 months to external file
### 2 months
my.output <- paste(Sys.Date(),"GSE89272_2Months.filter.integer.Hepa.Data.totalCount.txt",sep="_")
write.table(Hepa.Data.2Months, file = my.output, sep = "\t", eol="\r", na = "NA", row.names = T)
### 22 months
my.output <- paste(Sys.Date(),"GSE89272_22Months.filter.integer.Hepa.Data.totalCount.txt",sep="_")
write.table(Hepa.Data.22Months, file = my.output, sep = "\t", eol="\r", na = "NA", row.names = T)


###################################################################################################################
### 2. Make Sample Information Table in preparation for DESeq2 analysis
###################################################################################################################

# Make dataframe of sample information
### makes an object with genotypes in the order they appear in 'filter.integer.Hepa.Data.totalCount'
genotype <- as.factor(c(rep("Dwarf_Homo",4),rep("WT_Het",4),rep("Dwarf_Homo",4),rep("WT_Het",4)))
### makes an object with ages in the order they appear in 'filter.integer.Hepa.Data.totalCount'
ages <- as.factor(c(rep("2_months",8),rep("22_months",8)))
### creates dataframe with genotype & ages corresponding to columns in 'filter.integer.Hepa.Data.totalCount'
Hepa.annotation <- data.frame(genotype, ages)
### create new column with sample names, as they appear in 'filter.integer.Hepa.Data.totalCount'
Hepa.annotation$Samples <- colnames(filter.integer.Hepa.Data.totalCount)
### make Samples into row names
rownames(Hepa.annotation) <- Hepa.annotation$Sample
### Simplify Hepa.annotation for later use as 'colData'
Hepa.annotation <- Hepa.annotation[,c("genotype","ages")]


# Write Hepa annotation to external file
my.output <- paste(Sys.Date(),"GSE89272_Hepa.annotation.txt",sep="_")
write.table(Hepa.annotation, file = my.output, sep = "\t", eol="\r", na = "NA", row.names = T, col.names = T)


# Make dataframe of sample information for 2 months AND 22 months
### 2 months
Hepa.annotation.2Months <- Hepa.annotation[-c(9:16),]
### 22 months
Hepa.annotation.22Months <- Hepa.annotation[-c(1:8),]


# Write Hepa annotation for 2 months AND 22 months to external file
### 2 months
my.output <- paste(Sys.Date(),"GSE89272_2Months.Hepa.annotation.txt",sep="_")
write.table(Hepa.annotation.2Months, file = my.output, sep = "\t", eol="\r", na = "NA", row.names = T, col.names = T)
### 22 months
my.output <- paste(Sys.Date(),"GSE89272_22Months.Hepa.annotation.txt",sep="_")
write.table(Hepa.annotation.22Months, file = my.output, sep = "\t", eol="\r", na = "NA", row.names = T, col.names = T)


###################################################################################################################
### 3. Differential expression analysis using DESeq2 & normalization of data
###################################################################################################################

################### A. Prep data for DESeq2 analysis

# Double check that sample orders are consistent between count table and annotation dataframe
all(rownames(Hepa.annotation) == colnames(filter.integer.Hepa.Data.totalCount))
all(rownames(Hepa.annotation.2Months) == colnames(Hepa.Data.2Months))
all(rownames(Hepa.annotation.22Months) == colnames(Hepa.Data.22Months))


# Create DESeq2 object(s)
### NOTE: put the variable of interest at the end of the formula and make sure the control level is the first level
### Object for all data
ddsAll <- DESeqDataSetFromMatrix(countData = filter.integer.Hepa.Data.totalCount, 
                              colData = Hepa.annotation, 
                              design = ~ ages + genotype)
ddsAll$ages <- relevel(ddsAll$ages,"2_months")   # set '2_months' as the reference for 'ages', otherwise, DESeq will use whichever comes first in alphabetical order


### Object for 2 months
dds2Months <- DESeqDataSetFromMatrix(countData = Hepa.Data.2Months, 
                              colData = Hepa.annotation.2Months, 
                              design = ~ genotype)
### Object for 22 months
dds22Months <- DESeqDataSetFromMatrix(countData = Hepa.Data.22Months, 
                                     colData = Hepa.annotation.22Months, 
                                     design = ~ genotype)


# check to make sure you have the right samples
#as.data.frame(colData(ddsAll))
#as.data.frame(colData(dds2Months))
#as.data.frame(colData(dds22Months))


################### B. Run DESeq2, extract results, and normalize counts

# run DESeq2 algorithm(s)
ddsAll <- DESeq(ddsAll)
dds2Months <- DESeq(dds2Months)
dds22Months <- DESeq(dds22Months)


# extract DESeq2 results for 2 months AND 22 months separately
### the contrast arguments make the estimates the logarithmic fold change log2('Dwarf_Homo'/'WT_Het'); reference = 'WT_Het'
### NOTE: using contrast will additionally set to 0 the estimated LFC in a comparison of two groups, where all of the counts in the two groups are equal to 0
resAll <- results(ddsAll, contrast = c("genotype", "Dwarf_Homo", "WT_Het"))
### 2 Months
res2Months <- results(dds2Months, contrast = c("genotype", "Dwarf_Homo", "WT_Het"))
### 22 Months
res22Months <- results(dds22Months, contrast = c("genotype", "Dwarf_Homo", "WT_Het"))


# Write DESeq2 results to external file
### All
my.output <- paste(Sys.Date(),"GSE89272_ALL.Hepa.DwHomoVsWTHet.DESeq2.results.txt",sep="_")
write.table(resAll, file = my.output, sep = "\t", quote = F)
### 2 Months
my.output <- paste(Sys.Date(),"GSE89272_2Months.Hepa.DwHomoVsWTHet.DESeq2.results.txt",sep="_")
write.table(res2Months, file = my.output, sep = "\t", quote = F)
### 22 Months
my.output <- paste(Sys.Date(),"GSE89272_22Months.Hepa.DwHomoVsWTHet.DESeq2.results.txt",sep="_")
write.table(res22Months, file = my.output, sep = "\t", quote = F)


# normalize RNA-seq dataset for graphical purposes
### NOTE: 'getVarianceStabilizedData' requires first running 'estimateSizeFactors' AND 'estimateDispersions'; 'DESeq' does both steps, in addition to Negative Binomial GLM fitting and Wald statistics
### Running 'getVarianceStabilizedData' before or after running the 'DESeq' function does not make a difference in the count matrix produced
norm.ddsAll <- getVarianceStabilizedData(ddsAll)
norm.dds2Months <- getVarianceStabilizedData(dds2Months)
norm.dds22Months <- getVarianceStabilizedData(dds22Months)


################### C. Diagnoistic Plots

# Make MA plot: scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
### If data seems flat, you should NOT really trust the significance analysis
### All data
#plotMA(resAll, main = "All Data")
### 2 months
#plotMA(res2Months, main = "2 Months")
### 22 months
#plotMA(res22Months, main = "22 Months")


###################################################################################################################
### 4. PVclust - hierarchical clustering with bootstrap resampling
###################################################################################################################


# call clustering function, enter desired clustering parameters
# For REAL analysis: 100 or more bootstraps!!!!
result.pvclust <- pvclust(na.omit(norm.ddsAll), nboot=100)


# plot results to external pdf
pdf(paste(Sys.Date(),"GSE89272_Hepa_Data.PVclust.pdf",sep="_"))
plot(result.pvclust)
dev.off()


###################################################################################################################
### 5. Multidimensional scaling analysis
###################################################################################################################

# do MDS analysis
mds.result <- cmdscale(1-cor(norm.ddsAll,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
# exgtract first and second columns of the MDS analysis, x & y for each of the samples
x <- mds.result[, 1]
y <- mds.result[, 2]


# check to see if all subjects kept and in same order as in norm.ddsAll
### check dimension
dim(mds.result)
#16  2    ### all 16 subjects present
### check order
all(rownames(mds.result) == colnames(norm.ddsAll))
#TRUE     ### all subjects in the same order


# double check details for genotype and age
### check for 2 months
mds2m <- grep("_2m", rownames(mds.result))
mds2m
#1  2  3  4  5  6  7  8
length(mds2m)
#8
### check for 22 months
mds22m <- grep("_22m", rownames(mds.result))
mds22m
#9 10 11 12 13 14 15 16
length(mds22m)
#8
### check for homo
grep("Hom", rownames(mds.result))
#1  2  3  4  9 10 11 12
### check for het
grep("Het", rownames(mds.result))
#5  6  7  8 13 14 15 16


# specify details for graphs
# specifying that for the genotype, the first 4  will be pch 16, the next 4 pch 18
my.pch.genotype <- c(rep(16,4), rep(17,4))
# specifying that for age, the firt 8 will be coral, the next 5 darkblue
my.colors.age <- c(rep("red",8), rep("darkblue",8))


# create MDS plots and write to PDF
### Specify PDF name
pdf(paste(Sys.Date(),"GSE89272_Hepa_Data.MDSplot.pdf",sep="_"))
### make plot
plot(x, y,
     pch = my.pch.genotype, col = my.colors.age,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",cex=2)
legend("bottom", c("WT_Heterozygous","Dwarf_Homozygous"), col = "grey", pch = c(17,16), bty = 'n', pt.cex = 2)
legend("bottomright", c("2m","22m"), col = c("red","darkblue"), pch = 15, bty = 'n', pt.cex = 2)
### close out PDF file
dev.off()


###################################################################################################################
#### 6. Heatmap of ALL differentially expressed genes
##################################################################################################################

# Used to plot genes that are signficnatly changed

# remove NAs
resAll <- resAll[!is.na(resAll$padj),]

# make list of genes significnatly changed at FDR5
sig.genes <- rownames(resAll)[resAll$padj < 0.05] # FDR < 5%
# counts how many genes are significnatly changed at FDR5
my.num.genes <- length(sig.genes)


# heatmap drawing 
my.heatmap.out <- paste(Sys.Date(),"GSE89272_ALL_DWARFvsWT_Heatmap_significant_genes.pdf",sep="_")

pdf(my.heatmap.out, height = 20, width = 10, onefile = F)
my.heatmap.title <- paste("Significant Genotype Changes (FDR < 5%), ",my.num.genes, " genes", sep="")
pheatmap(norm.ddsAll[sig.genes,],   # takes only rows from 'norm.ddsAll' that are listed in 'sig.genes'
         cluster_cols = F,
         cluster_rows = T,  # so genes that are upregulated cluster together and genes that are downregualted cluster together
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 30)
dev.off()


# if you get an error 4, run "dev.off()" again, and then rerun PDF plot lines through "dev.off()"
# this is something that has been happening with the newest version of pheatmap


###################################################################################################################
#### 7. Prep data for Cluster Profiler
###################################################################################################################

################### A. 2 months

# remove NAs
res2Months <- res2Months[!is.na(res2Months$padj),]

# genes up/downregulated with FDR < 5% 
### 2 months
Months2.up <- intersect(rownames(res2Months)[res2Months$padj < 0.05], rownames(res2Months)[res2Months$log2FoldChange > 0])
Months2.down <- intersect(rownames(res2Months)[res2Months$padj < 0.05], rownames(res2Months)[res2Months$log2FoldChange < 0])
Months2.bg <- rownames(res2Months)


# count # of genes significantly up/down regulated and in background
length(Months2.up)
length(Months2.down)
length(Months2.bg)


# convert gene symbol to EntrezID
ids.Months2.up <- bitr(Months2.up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ids.Months2.down <- bitr(Months2.down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ids.Months2.bg <- bitr(Months2.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


################### B. 22 months

# Prep data for Cluster Profiler

# remove NAs
res22Months <- res22Months[!is.na(res22Months$padj),]


# genes up/downregulated with FDR < 5% 
### 22 months
Months22.up <- intersect(rownames(res22Months)[res22Months$padj < 0.05], rownames(res22Months)[res22Months$log2FoldChange > 0])
Months22.down <- intersect(rownames(res22Months)[res22Months$padj < 0.05], rownames(res22Months)[res22Months$log2FoldChange < 0])
Months22.bg <- rownames(res22Months)


# count # of genes significantly up/down regulated and in background
length(Months22.up)
length(Months22.down)
length(Months22.bg)


# convert gene symbol to EntrezID
ids.Months22.up <- bitr(Months22.up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ids.Months22.down <- bitr(Months22.down, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ids.Months22.bg <- bitr(Months22.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


################### B. Intersection

# genes flipping up/downregulated with age with FDR < 5% & b/g
upThenDownGenes <- intersect(Months2.up, Months22.down)
upThenDownGenes.bg <- c(Months2.up, Months22.down)

downThenUpGenes <- intersect(Months2.down, Months22.up)
downThenUpGenes.bg <- c(Months2.down, Months22.up)


# count # of genes flipped
length(upThenDownGenes)
length(downThenUpGenes)


# convert gene symbol to EntrezID
ids.upThenDownGenes <- bitr(upThenDownGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ids.upThenDownGenes.bg <- bitr(upThenDownGenes.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

ids.downThenUpGenes <- bitr(downThenUpGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
ids.downThenUpGenes.bg <- bitr(downThenUpGenes.bg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


###################################################################################################################
#### 8. Functional Enrichment with Cluster Profiler for FLIPPED genes
###################################################################################################################

################### A. GO BP over-representation test
# for upThenDownGenes 
ego.bp.upThenDown <- enrichGO(gene  = ids.upThenDownGenes$ENTREZID,
                              universe      = ids.upThenDownGenes.bg$ENTREZID,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = 'ENTREZID',
                              ont           = "BP", # ontogology we are testing here is BP: biological process
                              pAdjustMethod = "BH", # type of correction method we are using; normal way to calculate an FDR
                              pvalueCutoff  = 0.01, # not necessary, but a secondary precaution; before corrections, p-value needs to be below 0.01
                              qvalueCutoff  = 0.05, # q-values, tells FDR has to be below 0.05
                              readable      = TRUE)
# for downThenUpGenes
ego.bp.downThenUp <- enrichGO(gene  = ids.downThenUpGenes$ENTREZID,
                              universe      = ids.downThenUpGenes.bg$ENTREZID,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = 'ENTREZID',
                              ont           = "BP", # ontogology we are testing here is BP: biological process
                              pAdjustMethod = "BH", # type of correction method we are using; normal way to calculate an FDR
                              pvalueCutoff  = 0.01, # not necessary, but a secondary precaution; before corrections, p-value needs to be below 0.01
                              qvalueCutoff  = 0.05, # q-values, tells FDR has to be below 0.05
                              readable      = TRUE)


# write results to file
write.table(ego.bp.upThenDown@result,  file = paste(Sys.Date(),"GSE89272_upThenDown_GO_BP_UP_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
write.table(ego.bp.downThenUp@result, file = paste(Sys.Date(),"GSE89272_downThenUp_GO_BP_DWN_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# check to see if any significant GO terms
#count(ego.bp.upThenDown@result$p.adjust < 0.05)
#0
#count(ego.bp.downThenUp@result$p.adjust < 0.05)
#0


###################################################################################################################
#### 9. Functional Enrichment with Cluster Profiler for 2 Months (YOUNG)
###################################################################################################################

################### A. GO BP over-representation test
# for up-regulated genes
ego.bp.up.2Months <- enrichGO(gene  = ids.Months2.up$ENTREZID,
                      universe      = ids.Months2.bg$ENTREZID,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = 'ENTREZID',
                      ont           = "BP", # ontogology we are testing here is BP: biological process
                      pAdjustMethod = "BH", # type of correction method we are using; normal way to calculate an FDR
                      pvalueCutoff  = 0.01, # not necessary, but a secondary precaution; before corrections, p-value needs to be below 0.01
                      qvalueCutoff  = 0.05, # q-values, tells FDR has to be below 0.05
                      readable      = TRUE)
# for down regulated genes
ego.bp.dwn.2Months <- enrichGO(gene          = ids.Months2.down$ENTREZID,
                       universe      = ids.Months2.bg$ENTREZID,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)


# look at file
#head(ego.bp.up.2Months@result)
#head(ego.bp.dwn.2Months@result)


# write results to file
write.table(ego.bp.up.2Months@result,  file = paste(Sys.Date(),"GSE89272_2Months_GO_BP_UP_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
write.table(ego.bp.dwn.2Months@result, file = paste(Sys.Date(),"GSE89272_2Months_GO_BP_DWN_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# make some plots
pdf(paste(Sys.Date(),"GSE89272_2Months_GO_BP_dotplot_FDR5.pdf", sep = "_"))
dotplot(ego.bp.up.2Months,  x = "Count", title = "Upregulated GO BP Terms")
dotplot(ego.bp.dwn.2Months, x = "Count", title = "Downregulated GO BP Terms")
dev.off()

# gene concept visual
pdf(paste(Sys.Date(),"GSE89272_2Months_GO_BP_GeneConcept_Network_FDR5.pdf", sep = "_"))
cnetplot(ego.bp.up.2Months,  categorySize="pvalue")
cnetplot(ego.bp.dwn.2Months, categorySize="pvalue")
dev.off()

# cluster enrichment terms
pdf(paste(Sys.Date(),"GSE89272_2Months_GO_BP_Enrichment_Map_FDR5.pdf", sep = "_"))
emapplot(ego.bp.up.2Months)
emapplot(ego.bp.dwn.2Months)
dev.off()


################### B. KEGG over-representation test
# for up-regulated genes
kk.up.2Months  <- enrichKEGG(gene  = ids.Months2.up$ENTREZID,
                     universe      = ids.Months2.bg$ENTREZID,
                     keyType       = 'ncbi-geneid',
                     organism      = 'mmu',
                     pvalueCutoff  = 0.01)
# for down regulated genes
kk.dwn.2Months <- enrichKEGG(gene  = ids.Months2.down$ENTREZID,
                     universe      = ids.Months2.bg$ENTREZID,
                     keyType       = 'ncbi-geneid',
                     organism      = 'mmu',
                     pvalueCutoff  = 0.01)


# write results to file
write.table(kk.up.2Months@result,  file = paste(Sys.Date(),"GSE89272_2Months_KEGG_UP_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
write.table(kk.dwn.2Months@result, file = paste(Sys.Date(),"GSE89272_2Months_KEGG_DWN_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# make some plots
pdf(paste(Sys.Date(),"GSE89272_2Months_KEGG_dotplot_FDR5.pdf", sep = "_"))
dotplot(kk.up.2Months,  title = "Upregulated KEGG Pathways")
dotplot(kk.dwn.2Months, title = "Downregulated KEGG Pathways")
dev.off()

pdf(paste(Sys.Date(),"GSE89272_2Months_KEGG_GeneConcept_Network_FDR5.pdf", sep = "_"))
cnetplot(kk.up.2Months,  categorySize="pvalue")
cnetplot(kk.dwn.2Months, categorySize="pvalue")
dev.off()


pdf(paste(Sys.Date(),"GSE89272_2Months_KEGG_Enrichment_Map_FDR5.pdf", sep = "_"))
emapplot(kk.up.2Months)
emapplot(kk.dwn.2Months)
dev.off()


################### C. GO Gene Set Enrichment Analysis


# Prepare ranked GeneList using DESeq2 Wald statistic to rank genes
### Necessary for GSEA type analysis
res2Months$symbol <- rownames(res2Months)   # adds a column with Gene symbols taken from rownames
res2Months.Entrez <- merge(data.frame(res2Months), ids.Months2.bg, by.x = "symbol", by.y = "SYMBOL")   # adds Entrez ID to res2Months by matching the "symbol" column in res2Months with the "SYMBOL column in ids.Months2.bg 
Months2.geneList = res2Months.Entrez$stat   # makes object consisting of Wald statistic from DESeq analysis
names(Months2.geneList) = res2Months.Entrez$ENTREZID # identifies each Wald statistic with the ENTREZID from DESeq analysis with ENRTEZ column
Months2.geneList = sort(Months2.geneList, decreasing = TRUE)    # sorts gene list in decreasing order


# run upregulated and downregulated together
### cares about rank, rather than relative up and down regulation
### can detect the combined effect of many genes being regulated to a smaller amount consistently
go.bp.gsea.2Months <- gseGO(geneList     = Months2.geneList,
                      OrgDb              = org.Mm.eg.db,
                      keyType            = "ENTREZID",
                      ont                = "BP",
                      nPerm              = 1000,
                      minGSSize          = 100,
                      maxGSSize          = 500,
                      pvalueCutoff       = 0.05,
                      verbose            = FALSE)


# write results to file
write.table(go.bp.gsea.2Months@result, file = paste(Sys.Date(),"GSE89272_2Months_GO_BP_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
head(go.bp.gsea.2Months)   # to identify GO ID for interesting pathways to plot


###################  D. KEGG Gene Set Enrichment Analysis
kegg.gsea.2Months <- gseKEGG(geneList     = Months2.geneList,
                     organism             = 'mmu',
                     keyType              = 'ncbi-geneid',
                     nPerm                = 1000,
                     pvalueCutoff         = 0.05,
                     verbose              = FALSE)

# write results to file
write.table(kegg.gsea.2Months@result, file = paste(Sys.Date(),"GSE89272_2Months_KEGG_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
head(kegg.gsea.2Months)   # to identify KEGG ID for interesting pathways to plot


# get how many significant KEGG terms
nrow(kegg.gsea.2Months)


# GSEA plot of interesting pathways
pdf(paste(Sys.Date(),"GSE89272_2Months_GSEA_KEGG.pdf", sep = "_"))
gseaplot(kegg.gsea.2Months, geneSetID = "mmu04141", title = "Protein processing in endoplasmic reticulum")
gseaplot(kegg.gsea.2Months, geneSetID = "mmu04610", title = "Complement and coagulation cascades")
gseaplot(kegg.gsea.2Months, geneSetID = "mmu00983", title = "Drug metabolism - other enzymes")
dev.off()


###################################################################################################################
#### 10. Functional Enrichment with Cluster Profiler for 22 Months (OLD)
###################################################################################################################

################### A. GO BP over-representation test
# for up-regulated genes
ego.bp.up.22Months <- enrichGO(gene  = ids.Months22.up$ENTREZID,
                              universe      = ids.Months22.bg$ENTREZID,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = 'ENTREZID',
                              ont           = "BP", # ontogology we are testing here is BP: biological process
                              pAdjustMethod = "BH", # type of correct method we are using; normal way to calculate an FDR
                              pvalueCutoff  = 0.01, # not necessary, but a secondary precaution; before corrections, p-value needs to be below 0.01
                              qvalueCutoff  = 0.05, # q-values, tells FDR has to be below 0.05
                              readable      = TRUE)
# for down regulated genes
ego.bp.dwn.22Months <- enrichGO(gene         = ids.Months22.down$ENTREZID,
                               universe      = ids.Months22.bg$ENTREZID,
                               OrgDb         = org.Mm.eg.db,
                               keyType       = 'ENTREZID',
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE)


# look at file
#head(ego.bp.up.22Months@result)
#head(ego.bp.dwn.22Months@result)


# write results to file
write.table(ego.bp.up.22Months@result,  file = paste(Sys.Date(),"GSE89272_22Months_GO_BP_UP_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
write.table(ego.bp.dwn.22Months@result, file = paste(Sys.Date(),"GSE89272_22Months_GO_BP_DWN_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# make some plots
pdf(paste(Sys.Date(),"GSE89272_22Months_GO_BP_dotplot_FDR5.pdf", sep = "_"))
dotplot(ego.bp.up.22Months,  x = "Count", title = "Upregulated GO BP Terms")
dotplot(ego.bp.dwn.22Months, x = "Count", title = "Downregulated GO BP Terms")
dev.off()

# gene concept visual
pdf(paste(Sys.Date(),"GSE89272_22Months_GO_BP_GeneConcept_Network_FDR5.pdf", sep = "_"))
cnetplot(ego.bp.up.22Months,  categorySize="pvalue")
cnetplot(ego.bp.dwn.22Months, categorySize="pvalue")
dev.off()

# cluster enrichment terms
pdf(paste(Sys.Date(),"GSE89272_22Months_GO_BP_Enrichment_Map_FDR5.pdf", sep = "_"))
emapplot(ego.bp.up.22Months)
emapplot(ego.bp.dwn.22Months)
dev.off()


################### B. KEGG over-representation test
# for up-regulated genes
kk.up.22Months  <- enrichKEGG(gene  = ids.Months22.up$ENTREZID,
                             universe      = ids.Months22.bg$ENTREZID,
                             keyType       = 'ncbi-geneid',
                             organism      = 'mmu',
                             pvalueCutoff  = 0.01)
# for down regulated genes
kk.dwn.22Months <- enrichKEGG(gene  = ids.Months22.down$ENTREZID,
                             universe      = ids.Months22.bg$ENTREZID,
                             keyType       = 'ncbi-geneid',
                             organism      = 'mmu',
                             pvalueCutoff  = 0.01)


# write results to file
write.table(kk.up.22Months@result,  file = paste(Sys.Date(),"GSE89272_22Months_KEGG_UP_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
write.table(kk.dwn.22Months@result, file = paste(Sys.Date(),"GSE89272_22Months_KEGG_DWN_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")


# make some plots
pdf(paste(Sys.Date(),"GSE89272_22Months_KEGG_dotplot_FDR5.pdf", sep = "_"))
dotplot(kk.up.22Months,  title = "Upregulated KEGG Pathways")
dotplot(kk.dwn.22Months, title = "Downregulated KEGG Pathways")
dev.off()

pdf(paste(Sys.Date(),"GSE89272_22Months_KEGG_GeneConcept_Network_FDR5.pdf", sep = "_"))
cnetplot(kk.up.22Months,  categorySize="pvalue")
cnetplot(kk.dwn.22Months, categorySize="pvalue")
dev.off()


pdf(paste(Sys.Date(),"GSE89272_22Months_KEGG_Enrichment_Map_FDR5.pdf", sep = "_"))
emapplot(kk.up.22Months)
emapplot(kk.dwn.22Months)
dev.off()


################### C. GO Gene Set Enrichment Analysis


# Prepare ranked GeneList using DESeq2 Wald statistic to rank genes
### Necessary for GSEA type analysis
res22Months$symbol <- rownames(res22Months)   # adds a column with Gene symbols taken from rownames
res22Months.Entrez <- merge(data.frame(res22Months), ids.Months22.bg, by.x = "symbol", by.y = "SYMBOL")   # adds Entrez ID to res22Months by matching the "symbol" column in res22Months with the "SYMBOL column in ids.Months22.bg 
Months22.geneList = res22Months.Entrez$stat   # makes object consisting of Wald statistic from DESeq analysis
names(Months22.geneList) = res22Months.Entrez$ENTREZID # identifies each Wald statistic with the ENTREZID from DESeq analysis with ENRTEZ column
Months22.geneList = sort(Months22.geneList, decreasing = TRUE)    # sorts gene list in decreasing order


# run upregulated and downregulated together
### cares about rank, rather than relative up and down regulation
### can detect the combined effect of many genes being regulated to a smaller amount consistently
go.bp.gsea.22Months <- gseGO(geneList     = Months22.geneList,
                            OrgDb              = org.Mm.eg.db,
                            keyType            = "ENTREZID",
                            ont                = "BP",
                            nPerm              = 1000,
                            minGSSize          = 100,
                            maxGSSize          = 500,
                            pvalueCutoff       = 0.05,
                            verbose            = FALSE)


# write results to file
write.table(go.bp.gsea.22Months@result, file = paste(Sys.Date(),"GSE89272_22Months_GO_BP_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
head(go.bp.gsea.22Months)   # to identify GO ID for interesting pathways to plot


# get how many significant GO terms
nrow(go.bp.gsea.22Months)


# GSEA plot of interesting pathways
pdf(paste(Sys.Date(),"GSE89272_22Months_GSEA_GO.pdf", sep = "_"))
gseaplot(go.bp.gsea.22Months, geneSetID = "GO:0016570", title = "histone modification")
gseaplot(go.bp.gsea.22Months, geneSetID = "GO:0016571", title = "histone methylation")
dev.off()


###################  D. KEGG Gene Set Enrichment Analysis
kegg.gsea.22Months <- gseKEGG(geneList     = Months22.geneList,
                             organism             = 'mmu',
                             keyType              = 'ncbi-geneid',
                             nPerm                = 1000,
                             pvalueCutoff         = 0.05,
                             verbose              = FALSE)

# write results to file
write.table(kegg.gsea.22Months@result, file = paste(Sys.Date(),"GSE89272_22Months_KEGG_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
head(kegg.gsea.22Months)   # to identify GO ID for interesting pathways to plot


# get how many significant KEGG terms
nrow(kegg.gsea.22Months)


# GSEA plot of interesting pathways
pdf(paste(Sys.Date(),"GSE89272_22Months_GSEA_KEGG.pdf", sep = "_"))
gseaplot(kegg.gsea.22Months, geneSetID = "mmu04141", title = "Protein processing in endoplasmic reticulum")
gseaplot(kegg.gsea.22Months, geneSetID = "mmu04610", title = "Complement and coagulation cascades")
gseaplot(kegg.gsea.22Months, geneSetID = "mmu00140", title = "Steroid hormone biosynthesis")
dev.off()


###################################################################################################################
#### 11. Venn diagrams comparing DE genes in Dwarf_Homo at 2 months vs 22 months
###################################################################################################################

################### A. Upregulated Genes

# makes a list of all genes upregulated at 2 months, and all genes upregulated at 22 months
my.up.genes <- list("2 Months"            = Months2.up,
                    "22 Months"           = Months22.up)
my.Venn.up <- Venn(my.up.genes)  

my.Venn.out <- paste(Sys.Date(),"GSE89272_2Months_vs_22Months_Dwarf_UpVenn.pdf",sep="_")
pdf(my.Venn.out, height = 20, width = 10, onefile = F)
plot(my.Venn.up, doWeights= T, show = list(Universe = FALSE))
dev.off()

# is the overlap significant?

# get the "universe/background" ==> ALL genes you are able to detect in your experiment (upregulated, downregulated, neither)
my.Venn.bg <- intersect(rownames(res2Months),rownames(res22Months))

# makes a contingency table of #upregulated@both, #upregulated@2MonthsOnly, #upregulated@22MonthsOnly, allDetectableGenes-allGenesUpregulatedInEither
my.up.contingency.mat <- matrix(c(173,889,330,length(my.Venn.bg)-173-889-330),2,2)

# perform Fisher exact/hypergeometric test
fisher.test(my.up.contingency.mat, alternative = "greater")
# p-value < 2.2e-16
# overlap is greater than expected by chance


################### B. Down-regulated Genes

# makes a list of all genes upregulated at 2 months, and all genes upregulated at 22 months
my.down.genes <- list("2 Months"            = Months2.down,
                    "22 Months"           = Months22.down)
my.Venn.down <- Venn(my.down.genes)  

my.Venn.out <- paste(Sys.Date(),"GSE89272_2Months_vs_22Months_Dwarf_DownVenn.pdf",sep="_")
pdf(my.Venn.out, height = 20, width = 10, onefile = F)
plot(my.Venn.down, doWeights= T, show = list(Universe = FALSE))
dev.off()

# is the overlap significant?

# get the "universe/background" ==> ALL genes you are able to detect in your experiment (upregulated, downregulated, neither)
my.Venn.bg <- intersect(rownames(res2Months),rownames(res22Months))

# makes a contingency table of #upregulated@both, #upregulated@2MonthsOnly, #upregulated@22MonthsOnly, allDetectableGenes-allGenesUpregulatedInEither
my.down.contingency.mat <- matrix(c(333,903,290,length(my.Venn.bg)-333-903-290),2,2)

# perform Fisher exact/hypergeometric test
fisher.test(my.down.contingency.mat, alternative = "greater")
# p-value < 2.2e-16
# overlap is greater than expected by chance


################### C. All regulated Genes

# makes a list of all genes upregulated at 2 months, and all genes upregulated at 22 months
my.all.genes <- list("2 Months Up"            = Months2.up,
                     "22 Month Ups"            = Months22.up,
                     "2 Months Down"           = Months2.down,
                     "22 Months Down"          = Months22.down)
my.Venn.all <- Venn(my.all.genes)  

my.Venn.out <- paste(Sys.Date(),"GSE89272_2Months_vs_22Months_Dwarf_AllVenn.pdf",sep="_")
pdf(my.Venn.out, height = 20, width = 10, onefile = F)
plot(my.Venn.all, doWeights= T, show = list(Universe = FALSE))
dev.off()


