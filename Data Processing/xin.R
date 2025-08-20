library(scRNAseq)
sce <- XinPancreasData()

library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce)$symbol
gene.ids <- mapIds(edb, keys=symbols, 
                   keytype="SYMBOL", column="GENEID")

# Removing duplicated genes or genes without Ensembl IDs.
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
sce <- sce[keep,]
rownames(sce) <- gene.ids[keep]

saveRDS(sce, "xin.rds")

#######################################SCE QC##########
sce<-readRDS("xin.rds")

dataset <- "xin"
unfiltered <- sce

library(scater)
names(sce@assays)<- "counts"
stats <- perCellQCMetrics(sce)
qc <- quickPerCellQC(stats)
sce <- sce[,!qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

svg(paste0("QC/qcMetrics",dataset,".svg"))
gridExtra::grid.arrange(
  plotColData(unfiltered, x="donor.id", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, x="donor.id", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  #plotColData(unfiltered, x="donor", y="altexps_ERCC_percent",
   #           colour_by="discard") + ggtitle("ERCC percent"),
  ncol=2
)
dev.off()

colSums(as.matrix(qc))

saveRDS(sce, "QC/xinQC.rds")

########
