library(scater)
library(scran)
library(limma)
library(BiocParallel)
library(BiocNeighbors)

xinQC <- readRDS("QC/xinQC.rds")
baronQC <- readRDS("QC/baronQC.rds")
segQC <- readRDS("QC/segQC.rds")
muraroQC <- readRDS("QC/muraroQC.rds")

keep_genes <- intersect(intersect(intersect(rownames(xinQC),rownames(baronQC)),rownames(segQC)),rownames(muraroQC))

xinQC <- xinQC[match(keep_genes, rownames(xinQC)), ]
baronQC <- baronQC[match(keep_genes, rownames(baronQC)), ]
segQC <- segQC[match(keep_genes, rownames(segQC)), ]
muraroQC <- muraroQC[match(keep_genes, rownames(muraroQC)), ]

# Combine Datasets

counts_panc <- cbind(counts(xinQC), counts(baronQC))
counts_panc <- cbind(counts_panc, counts(segQC))
counts_panc <- cbind(counts_panc, counts(muraroQC))

label.xin <- as.data.frame(xinQC$cell.type)
colnames(label.xin) <-"cell.type"
label.baron <- as.data.frame(baronQC$label)
colnames(label.baron) <-"cell.type"
label.seg <- as.data.frame(segQC$CellType)
colnames(label.seg) <-"cell.type"
label.muraro <- as.data.frame(muraroQC$label)
colnames(label.muraro) <-"cell.type"

label_panc <- rbind(label.xin,label.baron)
label_panc <- rbind(label_panc,label.seg)
label_panc <- rbind(label_panc,label.muraro)


colData(xinQC)$age <- NULL
colData(xinQC)$Sample.name <- NULL
colData(xinQC)$condition <- NULL
colData(xinQC)$ethnicity <- NULL
colData(xinQC)$gender <- NULL

names(colData(xinQC))

names(colData(baronQC))[1] ="donor.id"
names(colData(baronQC))[2] ="cell.type"

colData(segQC)$Disease <- NULL
colData(segQC)$Quality <- NULL

names(colData(segQC))[1] <- "cell.type"
names(colData(segQC))[2] <- "donor.id"

colData(muraroQC)$plate <- NULL

names(colData(muraroQC))[1] <- "cell.type"
names(colData(muraroQC))[2] <- "donor.id"


merged.sce <- SingleCellExperiment( 
  assays = list(counts = counts_panc),  
  rowData = rowData(xinQC), # same as rowData(pbmc4k) 
  colData = rbind(colData(xinQC), colData(baronQC),colData(segQC),colData(muraroQC)) 
) 

saveRDS(merged.sce, "merged_sce.rds")


# Normalization

merged_sce <- read_csv("eliminatedClasses/merged_sce.csv")
merged_sce<- as.data.frame(t(merged_sce))

merged.sce.eliminated <- SingleCellExperiment( 
  assays = list(counts = merged_sce),  
  rowData = rowData(merged.sce), # same as rowData(pbmc4k) 
  colData = donor_celltype_dataset 
)

merged.sce <- merged.sce.eliminated

library(scran)
clusters <- quickCluster(merged.sce)
merged.sce <- computeSumFactors(merged.sce, clusters=clusters)
merged.sce <- logNormCounts(merged.sce) 
summary(sizeFactors(merged.sce))

# Batch Effect Removal
colData(merged.sce)$batch <- donor.id$batch
limma_corrected <- limma::removeBatchEffect(logcounts(merged.sce), batch = merged.sce$batch)
assay(merged.sce, "logcounts_limma") <- limma_corrected ## add new assay


# Feature Selection

merged.sce.fs <- modelGeneVar(merged.sce)
chosen <- getTopHVGs(merged.sce.fs )#prop=0.1
str(chosen)

merged.sce.hvg <- merged.sce[chosen,]
dim(merged.sce.hvg)

saveRDS(merged.sce.hvg,"merged.sce.hvg.rds")


library(rhdf5)
library(scrattch.io)
write_dgCMatrix_csv(t(counts(merged.sce.hvg)), "merged.sce.hvg.csv", col1_name = "gene",
                    chunk_size = 1000)


## Remove 'unclear' labels from dataset

corrected_logcounts_limma <- read_csv("corrected_logcounts_limma.csv")
donor_celltype_dataset <- read_delim("donor-celltype-dataset.csv", 
                                     +     delim = ";", escape_double = FALSE, trim_ws = TRUE)











