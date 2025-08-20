library(scRNAseq)

sce.seger <- SegerstolpePancreasData()

library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce.seger)$symbol
ens.id <- mapIds(edb, keys=symbols, keytype="SYMBOL", column="GENEID")
ens.id <- ifelse(is.na(ens.id), symbols, ens.id)

# Removing duplicated rows.
keep <- !duplicated(ens.id)
sce.seger <- sce.seger[keep,]
rownames(sce.seger) <- ens.id[keep]

emtab.meta <- colData(sce.seger)[,c("cell type", "disease",
                                    "individual", "single cell well quality")]
colnames(emtab.meta) <- c("CellType", "Disease", "Donor", "Quality")
colData(sce.seger) <- emtab.meta

sce.seger$CellType <- gsub(" cell", "", sce.seger$CellType)
sce.seger$CellType <- paste0(
  toupper(substr(sce.seger$CellType, 1, 1)),
  substring(sce.seger$CellType, 2))

saveRDS(sce.seger, "seg.rds")

#######################################SCE QC##########
dataset<- "seg"
unfiltered <- sce.seger

low.qual <- sce.seger$Quality == "low quality cell"

library(scater)
stats <- perCellQCMetrics(sce.seger)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent",
                     batch=sce.seger$Donor,
                     subset=!sce.seger$Donor %in% c("HP1504901", "HP1509101"))

sce.seger <- sce.seger[,!(qc$discard | low.qual)]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

svg(paste0("QC/qcMetrics",dataset,".svg"))
gridExtra::grid.arrange(
  plotColData(unfiltered, x="Donor", y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count") +
    theme(axis.text.x = element_text(angle = 90)),
  plotColData(unfiltered, x="Donor", y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features") +
    theme(axis.text.x = element_text(angle = 90)),
  plotColData(unfiltered, x="Donor", y="altexps_ERCC_percent",
              colour_by="discard") + ggtitle("ERCC percent") +
    theme(axis.text.x = element_text(angle = 90)),
  ncol=2
)
dev.off()

colSums(as.matrix(qc))

saveRDS(sce.seger, "QC/segQC.rds")











