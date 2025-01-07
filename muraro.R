#muraro
library(scRNAseq)

sceM <- MuraroPancreasData()

df <- data.frame(counts(sceM))
df_label<- sceM$label

#sceM <- sceM[,!is.na(sceM$label)]
#sceM <- sceM[,-c(which(sceM$label == "unclear"))]

#counts <- assay(sceM, "counts")
#libsizes <- colSums(counts)
#size.factors <- libsizes/mean(libsizes)
#logcounts(sceM) <- log2(t(t(counts)/size.factors) + 1)

#saveRDS(sceM, "muraro.rds")

#####################
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
gene.symb <- sub("__chr.*$", "", rownames(sceM))
gene.ids <- mapIds(edb, keys=gene.symb, 
                   keytype="SYMBOL", column="GENEID")

# Removing duplicated genes or genes without Ensembl IDs.
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sceM <- sceM[keep,]
rownames(sceM) <- gene.ids[keep]

saveRDS(sceM, "muraro.rds")

#######################################SCE QC##########
sce<-readRDS("muraro.rds")

dataset <- "muraro"
unfiltered <- sce

library(scater)
stats <- perCellQCMetrics(sce)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
                     batch=sce$donor, subset=sce$donor!="D28")
sce <- sce[,!qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

svg(paste0("QC/qcMetrics",dataset,".svg"))
gridExtra::grid.arrange(
  plotColData(unfiltered, x="donor", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, x="donor", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, x="donor", y="altexps_ERCC_percent",
              colour_by="discard") + ggtitle("ERCC percent"),
  ncol=2
)
dev.off()

colSums(as.matrix(qc))

saveRDS(sceM, "QC/muraroQC.rds")

########
