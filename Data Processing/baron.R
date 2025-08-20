### DATA
# human1
h1 <- read.csv("GSM2230757_human1_umifm_counts.csv", header = T)
rownames(h1) <- h1[,1]
labels_h1 <- as.character(h1$assigned_cluster)
h1 <- h1[,4:ncol(h1)]
h1 <- t(h1)
# human2
h2 <- read.csv("GSM2230758_human2_umifm_counts.csv", header = T)
rownames(h2) <- h2[,1]
labels_h2 <- as.character(h2$assigned_cluster)
h2 <- h2[,4:ncol(h2)]
h2 <- t(h2)
# human3
h3 <- read.csv("GSM2230759_human3_umifm_counts.csv", header = T)
rownames(h3) <- h3[,1]
labels_h3 <- as.character(h3$assigned_cluster)
h3 <- h3[,4:ncol(h3)]
h3 <- t(h3)
# human4
h4 <- read.csv("GSM2230760_human4_umifm_counts.csv", header = T)
rownames(h4) <- h4[,1]
labels_h4 <- as.character(h4$assigned_cluster)
h4 <- h4[,4:ncol(h4)]
h4 <- t(h4)

# merge data
h <- cbind(h1, h2, h3, h4)

### ANNOTATIONS
# human
h_ann <- data.frame(
  human = c(
    rep(1, length(labels_h1)), 
    rep(2, length(labels_h2)), 
    rep(3, length(labels_h3)), 
    rep(4, length(labels_h4))
  ),
  cell_type1 = c(labels_h1, labels_h2, labels_h3, labels_h4))
rownames(h_ann) <- colnames(h)


### SINGLECELLEXPERIMENT
source("../utils/create_sce.R")
h_sceset <- create_sce_from_counts(h, h_ann)
saveRDS(h_sceset, "baron-human.rds")



############################################### SCE ###########

sce.human <- BaronPancreasData()

library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
gene.symb <- sub("__chr.*$", "", rownames(sce.human))
gene.ids <- mapIds(edb, keys=gene.symb, 
                   keytype="SYMBOL", column="GENEID")

# Removing duplicated genes or genes without Ensembl IDs.
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sce.human <- sce.human[keep,]
rownames(sce.human) <- gene.ids[keep]

saveRDS(sce.human, "baron.rds")


#############Sce QC #####

dataset <- "baron"
unfiltered <- sce.human

is.mito <- which(gene.ids == "MT")

library(scater)
stats <- perCellQCMetrics(sce.human)
qc <- quickPerCellQC(stats)
sce.human <- sce.human[,!qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

svg(paste0("QC/qcMetrics",dataset,".svg"))
gridExtra::grid.arrange(
  plotColData(unfiltered, x="donor", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, x="donor", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  #plotColData(unfiltered, x="donor", y="altexps_ERCC_percent",
   #           colour_by="discard") + ggtitle("ERCC percent"),
  ncol=2
)
dev.off()

colSums(as.matrix(qc))

saveRDS(sce.human, "QC/baronQC.rds")

########












