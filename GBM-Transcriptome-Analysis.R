setwd("C:/Users/HP/Documents/VisualizationPractice/")

#install libraries
##TCGA Biolinks
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

##edgeR
BiocManager::install("edgeR")
##limma
BiocManager::install("limma", force = TRUE)
##EDAseq

BiocManager::install("EDASeq")
install.packages("gplots")

library("TCGAbiolinks")
library(SummarizedExperiment) #installed alongside TCGAbiolinks
library(biomaRt)
library("limma")
library("edgeR")
library("EDASeq")
library("gplots")

#project information
getProjectSummary("TCGA-GBM")
?GDCquery

#dOWNLOAD and preprocess data
gbmQ <- GDCquery(project = "TCGA-GBM",
                 data.category = "Transcriptome Profiling",
                 data.type = "Gene Expression Quantification")

GDCdownload(gbmQ)
gbm.data <- GDCprepare(gbmQ)
View(gbm.data)

#explore some metadata information

gbm.data$race
gbm.data$tumor_descriptor
gbm.data$barcode

#We can create a simple metadata for our use case

SimpleMeta <- data.frame("barcode" = gbm.data$barcode,
                         "race" = gbm.data$race,
                         "tumor_type" = gbm.data$tumor_descriptor)


#=========================002=======================================

#select unstranded dataset
gbm.raw.data <- assays(gbm.data)
dim(gbm.raw.data$unstranded) #datasize

#lets downsize our data to 5 primary and 5 recurring data alone.

selectedBarcodes <- c(subset(SimpleMeta, tumor_type == "Recurrence")$barcode[c(1:5)],
                      subset(SimpleMeta, tumor_type == "Primary")$barcode[c(1:5)])

selectedData <- gbm.raw.data$unstranded[, c(selectedBarcodes)]
dim(selectedData)

#data normalization and filtering
normData <- TCGAanalyze_Normalization(tabDF = selectedData,
                                      geneInfo = geneInfoHT,
                                      method = "geneLength")

filtData <- TCGAanalyze_Filtering(tabDF = normData,
                                  method = "quantile",
                                  qnt.cut = 0.25)


#=========================003=======================================

#we can decide to visualise it as a heatmap, but 6000 is too much and a waste of time
#we will pick the top n (50, 100, etc) differentially expressed genes

#Differentially expressed analysis(DEA)
selectResults <- 
  TCGAanalyze_DEA(mat1 = filtData[, c(selectedBarcodes)[1:5]],
                  mat2 = filtData[, c(selectedBarcodes)[6:10]],
                  Cond1type = 'Recurrence',
                  Cond2type = "Primary",
                  pipeline = "edgeR", #can be "limma,
                  fdr.cut = 0.01,
                  logFC.cut = 2)

plot(selectResults$logFC, -log10(selectResults$FDR))

#differentially expressed analysis with treatment levels

selectResults.Level <- TCGAanalyze_LevelTab(selectResults, "Recurrence", "primary",
                                            filtData[, c(selectedBarcodes)[1:5]],
                                            filtData[, c(selectedBarcodes)[6:10]])

#=========================004=======================================

#how can we visualize with a heatmap
head(selectResults.Level)
dim(selectResults.Level)

heat.data <- filtData[rownames(selectResults.Level),]

#color the plot by the kind of tumor

head(heat.data)
cancer.type <- c(rep("Recurrence", 5), rep('Primary', 5))

ccodes <- c()

for (i in cancer.type) {
  if (i == "Recurrence"){
    ccodes <- c(ccodes, "red")
  }else {
    ccodes <- c(ccodes, "blue")
  }
}

heatmap.2(x = as.matrix(heat.data),
          col = hcl.colors(10, palette = 'Blue-Red 2'), #search hcl colors in r
          Rowv = F, Colv = T,
          scale = 'row', 
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.5, cexCol = 1,
          main = "Heatmap of Recurrence vs Primary",
          na.color = 'black',
          ColSideColors = ccodes)



#=========================005=======================================

#Principal Component Analysis

pca <- TCGAvisualize_PCA(filtData,
                         selectResults.Level,
                         ntopgenes = 100,
                         selectedBarcodes[1:5],
                         selectedBarcodes[6:10])
          
#=========================005=======================================

#Enrichment Analysis
#View the volcano plot first

plot(x = selectResults.Level$logFC, y = -log10(selectResults.Level$FDR))

upreg.genes <- rownames(subset(selectResults.Level, logFC >2))
dnreg.genes <- rownames(subset(selectResults.Level, logFC < -2))

#convert ensemble IDs to gene IDs using biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(mart)
head(attributes)

upreg.genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                     filters = 'ensembl_gene_id', 
                     values = upreg.genes,
                     mart = mart)$hgnc_symbol

dnreg.genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                     filters = 'ensembl_gene_id', 
                     values = dnreg.genes,
                     quote = FALSE,
                     mart = mart)$hgnc_symbol

#now we can perform enrinchment analysis for both
up.EA <- TCGAanalyze_EAcomplete(TFname="UPREGULATED", upreg.genes)
dn.EA <- TCGAanalyze_EAcomplete(TFname="DOWNREGULATED", dnreg.genes)

#now we can visullize our results
TCGAvisualize_EAbarplot(tf = rownames())


#=========================007=======================================
#MACHINE LEARNING
#change of data
table(gbm.data$tumor_descriptor)


#we will be using this low grade glioblasstoma dataset

lgg.data <- readRDS('LGGrnaseq.rds')
meta <- readRDS("patient2LGGsubtypes.rds")

#=========================007=======================================
#predicting tumor status as either methylated or non-methylated

install.packages("caret")
install.packages("DALEX")
install.packages("pROC")

library(caret)
library(pROC)
library(DALEX)
set.seed(34567)

#UNSUPERVISED MACHINE LEARNING
##k- nearest neighbour
#how many samples do we have under each dataset

table(meta$subtype)

