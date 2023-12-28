library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(limma)
library(DESeq2)
library(scater)
library(harmony)
library(DoubletFinder)
library(future)
library(tidyverse)

data_dir <- 'mm1_filtered_feature_bc_matrix'
list.files(data_dir)
data <- Read10X(data.dir = data_dir)
mm1 <- CreateSeuratObject(counts = data, project = "mm1", min.cells = 1)

sample_list <- c('gom', 'pmm', 'sat', 'lung', 'spleen', 'kidney', 'mm1')
mm1 = merge(x = gom, y = c(pmm, sat, lung, spleen, kidney, mm1),
            add.cell.ids = sample_list, project = 'mm1')

# Calculate the percentage of mitochondrial genes, 13 mitochondrial genes in total
mm1[["ND6"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ND6")
mm1[["ND5"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ND5")
mm1[["ND4L"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ND4L")
mm1[["ND4"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ND4")
mm1[["ND3"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ND3")
mm1[["ND2"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ND2")
mm1[["ND1"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ND1")
mm1[["CYTB"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "CYTB")
mm1[["COX3"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "COX3")
mm1[["COX2"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "COX2")
mm1[["COX1"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "COX1")
mm1[["ATP8"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ATP8")
mm1[["ATP6"]] <- PercentageFeatureSet(mm1, assay = "RNA", feature = "ATP6")
#or 
MT.genes <- rownames(m1)[grep("^MT-",rownames(m1))]
m1[["percent.MT"]] <- PercentageFeatureSet(m1, pattern = "^MT-")
p1 <- VlnPlot(m1, features = c("nCount_RNA", "nFeature_RNA", "percent.MT"), ncol = 3, pt.size = 0.1) + NoLegend() 
ggsave("m1-rawdata-volin.pdf", plot=p1, width = 20, height =8)
p2 <- FeatureScatter(m1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("m1-rawdatacounts_gene.pdf", plot=p2, width = 8, height =8)

# Calculate and add total mitochondrial gene percentage
head(mm1@meta.data, 5)  # View the first 5 rows of metadata
mm1[["percent.mt"]] <- apply(mm1@meta.data[, 5:17], 1, sum)  # Sum to get the total mitochondrial percentage
mm1@meta.data <- mm1@meta.data[, -c(5:17)]  # Remove intermediate data
head(mm1@meta.data, 5)  # View the updated metadata
write.csv(mm1@meta.data, file = "muscle-cells-metadata-rawdate.csv")

# Create a feature scatter plot
p2 <- FeatureScatter(mm1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("mm1-beforecounts_gene-dimplot.pdf", plot = p2, width = 8, height = 8)

min.cells <- 10
counts <- GetAssayData(mm1, layer = "counts")
num.cells <- rowSums(counts > 0)  # Sum of cells expressing each gene, returns a matrix
genes.use <- names(num.cells[which(num.cells >= min.cells)])  # Identify genes expressed in at least 10 cells
mykeepgene <- c(1:nrow(mm1))[rownames(mm1) %in% as.character(genes.use)]
mm2 <- subset(mm1, features = mykeepgene)
mm2

mm3 <- subset(mm2, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & nCount_RNA > 2000 & nCount_RNA < 50000 & percent.mt < 10)
mm3
saveRDS(mm3, file = "mm3.rds")
mm3 <- readRDS("/media/disk/chenlong/bamapig/all_data/mm3.rds")
write.csv(mm3@meta.data, file = "muscle-cells-metadata-fliter.csv")

# Create violin plots
p1 <- VlnPlot(mm3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3) + NoLegend()
ggsave("muscle-fliter-volin.pdf", plot = p1, width = 15, height = 8)

seurat_phase <- NormalizeData(mm3)
seurat_phase <- NormalizeData(mm3)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(seurat_phase), 10) ### Select the top 10 highly variable genes for plotting
plot1 <- VariableFeaturePlot(seurat_phase)    ### Create a plot of highly variable genes without labels
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size = 2.5) ### Add the top 10 genes to the plot
plot <- CombinePlots(plots = list(plot1, plot2), legend = "bottom") 
ggsave("VariableFeatures-plot.pdf", plot = plot, width = 12, height = 8)
seurat_phase <- ScaleData(seurat_phase, features = rownames(seurat_phase))
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase))
p1 <- ElbowPlot(seurat_phase, ndims = 50, reduction = "pca") ### Choose an appropriate number of PCs, usually at the inflection point.
ggsave("muscle-PCA-ElbowPlot.pdf", plot = p1, width = 5, height = 5)
seurat_phase <- RunHarmony(seurat_phase, group.by.vars = "orig.ident")
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, reduction = "harmony")
seurat_phase <- RunUMAP(seurat_phase, dims = 1:30, reduction = "harmony")
seurat_phase <- FindClusters(object = seurat_phase, resolution = 0.4)
saveRDS(seurat_phase, file = "yunashi-cellcycle.rds")
exp.mat <- read.table(file = "/media/18t/chenlong/bamapig/raw data/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes
g2m.genes
counts <- seurat_phase
cell_cycle_genes <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
                      "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
                      "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
                      "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
                      "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
                      "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
                      "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
                      "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
                      "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
                      "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
                      "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
                      "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
                      "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
                      "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

# Use the cell cycle gene list to select relevant genes
test <- counts[cell_cycle_genes, ]
rowSums(test)
seurat_phase <- CellCycleScoring(seurat_phase, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(seurat_phase@meta.data, 3)
write.csv(seurat_phase@meta.data, file = "seurat_phase@meta.data-cellscore.csv")
p1 <- DimPlot(seurat_phase, reduction = "umap", group.by = "Phase", split.by = "orig.ident")
ggsave("bamapig-PCA-cell-cycle-SCThou-cellcycle1.pdf", plot = p1, width = 25, height = 5)
# s.genes: "SLBP", "NASP"    g2m.genes: "TUBB4B", "NUSAP1"
p1 <- RidgePlot(mergehs, features = c("SLBP", "NASP", "TUBB4B", "NUSAP1"), ncol = 3)
ggsave("liver-cell-cycle-marker-heatmap.pdf", plot = p1, width = 8, height = 4)

gc()
# pK Identification (no ground-truth)
sweep.res.list <- paramSweep_v3(seurat_phase, PCs = 1:30, num.cores = 10)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
## Homotypic Doublet Proportion Estimate
homotypic.prop <- modelHomotypic(seurat_phase@active.ident)
nExp_poi <- round(0.07 * length(colnames(seurat_phase)))   ## Assuming 7% doublet formation rate
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
seurat_phase <- doubletFinder_v3(seurat_phase, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE)
# head(test@meta.data)
doublet <- grep("DF.classifications", names(seurat_phase@meta.data), value = TRUE)
table(seurat_phase[[doublet]])
# plot
p1 <- DimPlot(seurat_phase, reduction = "umap", group.by = doublet)
ggsave("spleen-DoubletFinder.pdf", plot = p1, width = 5, height = 5)
#
# Save clean singlet data
keep <- c(1:nrow(seurat_phase@meta.data))[seurat_phase[[doublet]] == "Singlet"]
seurat_phase2 <- subset(seurat_phase, cells = c(keep))
saveRDS(seurat_phase2, file = "celltype-cleandata-singlet.rds")

cluster <- table(seurat_phase@active.ident)
write.csv(cluster, file="bamapig-cluster-r0.4-cell-number.csv")
write.csv(seurat_phase@meta.data, file = "bamapig-clean-data-cluster.csv")
p1 <- VlnPlot(seurat_phase, features="nCount_RNA", group.by="seurat_clusters", same.y.lims=TRUE, fill.by="seurat_clusters", pt.size = 0.1) + NoLegend()
ggsave("muscle-cluster-nCount_SCT-VlnPlot.pdf", plot=p1, width = 20, height =5)
p1 <- VlnPlot(seurat_phase, features="nFeature_RNA", group.by="seurat_clusters", same.y.lims=TRUE, fill.by="seurat_clusters", pt.size = 0.1) + NoLegend()
ggsave("muscle-cluster-nFeature_SCT-VlnPlot.pdf", plot=p1, width = 20, height =5)

p1 <- DimPlot(seurat_phase, reduction = "umap", label = TRUE, repel = TRUE)  
ggsave("bamapig-r0.4-dims20-UMAP.pdf", plot=p1, width = 6, height =5)
p1 <- DimPlot(seurat_phase, reduction = "umap", group.by="orig.ident",label = TRUE, repel = TRUE)  
ggsave("bamapig-r0.4-dims20-UMAP-orig.pdf", plot=p1, width = 6, height =5)

seurat_phase <- readRDS(file="spleen-cleandata-singlet.rds")
plan(multisession, workers=4)           #The server uses 4 cores
difgene1 <- FindAllMarkers(object=seurat_phase, slot= "counts", min.pct=0.25, logfc.threshold=0.25, only.pos=TRUE)
difgene1 <- subset(difgene1, p_val_adj < 0.05)
write.csv(difgene1, file="bamapig-r0.4-cluster-marker.csv")

difgene2 <- subset(difgene1, p_val_adj < 0.05 & avg_log2FC > 2 & pct.1 >0.6)
write.csv(difgene2, file="bamapig-r0.4-cluster-marker-2-0.6.csv")

Idents(seurat_phase) <- seurat_phase$celltype
aa <- seurat_phase@meta.data
aa$celltype %>% table
#seurat_phase <- subset(seurat_phase, idents = c("24","25","26"), invert = TRUE)
#old
new.cluster.ids<-c("0"="T cell IL7R+","1"="T cell CRTAM+","2"="T cell GNLY+","3"="Fibroblasts124","4"="Proximal tubule cells","5"="T cell CRTAM+","6"="T cell/NK cell KLRB1+",
                   "7"="Macrophage","8"="II myofiber","9"="I myofiber","10"="T cell IL7R+","11"="Endothelial_cells12","12"="Fibroblasts123",
                   "13"="B_cells1","14"="Macrophage","15"="Distal tubule cells","16"="T cell CRTAM+","17"="B_cells2","18"="B_cells1",
                   "19"="Fibroblasts3","20"="Macrophage","21"="Endothelial_cells2","22"="Proximal tubule cells","23"="Macrophage")
#new
#seurat_phase[["celltype"]] <- NULL
new.cluster.ids<-c("0"="T cells GNLY+","1"="T cells CRTAM+","2"="T cells IL7R+","3"="Fibroblasts1","4"="myofiber",
                   "5"="T cells IL7R+","6"="Proximal tubule cells","7"="Macrophage","8"="Endothelial cells",
                   "9"="T cells IL7R+","10"="T cells KLRB1+","11"="Fibroblasts2","12"="B cells1",
                   "13"="Distal tubule cells","14"="Macrophage","15"="Adipocytes","16"="T cells CRTAM+",
                   "17"="B cells2","18"="Erythroid cells","19"="B cells1","20"="SMLC","21"="Endothelial cells",
                   "22"="T cells CRTAM+")
seurat_phase <- RenameIdents(seurat_phase, new.cluster.ids)
seurat_phase$celltype <- seurat_phase@active.ident 
head(seurat_phase@active.ident,5)
write.csv(seurat_phase@meta.data, file = "bamapig-cells-metadata-celltypehou.csv")
DimPlot(seurat_phase,group.by = "Phase", reduction = "umap")
p1 <- DimPlot(seurat_phase,group.by = "celltype", reduction = "umap", label = TRUE, repel = TRUE)
ggsave("r0.4-dims10-UMAP-jiandinghou_new.pdf", plot=p1, width =10, height =8)
p1 <- DimPlot(seurat_phase, reduction="umap", group.by="celltype", split.by="orig.ident", label = TRUE, repel = TRUE)
ggsave("fat-r0.1-dims10-UMAP-7apart.pdf", plot=p1, width = 42, height =5)

p1 <- DotPlot(object =seurat_phase,group.by = "celltype",features =c("GNLY","CRTAM","IL7R","ABCA6","ABCA8","TNNT3","TTLL11",
                                                                     "ALDOB","GPX3","RBPJ","FRMD4B","PREX2","PECAM1","KLRB1",
                                                                     "FN1","CDON","SMC6","CD79B","TMEM213","WFDC2","ADIPOR2","PLIN1",
                                                                     "JCHAIN","MZB1","AHSP","ALAS2","SYNPO2","DAAM2"))+    RotatedAxis()
ggsave("bamapig-r0.4-qipaotu-adi-new.pdf", plot=p1, width = 15, height =5)


p1 <- DotPlot(object =seurat_phase,group.by = "celltype",features =c("IL7R","EVL","CRTAM","MLF1","GNLY","KLRB1","NCR1","KLRF1",
                                                                     "SMC6","CD79B","CD79A","MS4A1","JCHAIN","MZB1","TXNDC5",
                                                                     "POU2AF1","DERL3",
                                                                     "TNNT3","MYBPC2","TTLL11","ACTN3","GADL1","ATP2A1",
                                                                     "ATP2A2","MYH7","ARPP21","TECRL","ESRRG","MYOZ2","TNNT1","MYH7B",
                                                                     "ABHD18","MYL3",
                                                                     "ABCA6","ABCA8","SPARCL1","KAZN","HMCN2","DKK2","NID1","BICC1","COL6A3","FBN1",
                                                                     "SLIT3","FGF2","HMCN1","PID1","DCN","PDGFRB","ACTA2","CACNB2","FN1",
                                                                     "GLI3","GRK5","CDON","COL1A1","NEGR1","COL1A2","COL3A1","FLT1",
                                                                     "ALDOB","GPX3","ASS1","PCK1","BHMT","FBP1","SLC13A3","SLC13A1","DCXR",
                                                                     "SLC34A1","LRP2","UMOD","DEFB1","CA2","WFDC2","ATP6V0D2",
                                                                     "ATP6V1C2","ATP6V1G3","TMEM213",
                                                                     "MECOM","ADAMTS9","PTPRB","PREX2","CYYR1","CDH5","PECAM1","ERG",
                                                                     "VWF","KDR","CLEC14A","ID1","ENG","PROX1",
                                                                     "RBPJ","FRMD4B","MTSS1","CD163","MRC1","CD14","CSF1R","PYCARD","CD86","LYZ",
                                                                     "ADIPOR2","PLIN1","ADIPOQ","ADAMTS12","PLIN4","LIPE"))+    RotatedAxis()
ggsave("bamapig-r0.5-qipaotu-adi.pdf", plot=p1, width = 33, height =5)
#T CELLS
P1 <- FeaturePlot(seurat_phase,features = c("CRTAM","MLF1","GNLY","KLRB1","NCR1","KLRF1","IL7R","CCR7","EVL"),reduction = 'umap')
ggsave(filename = 'T cells-maker-features.pdf', plot=P1, width = 15, height = 15)
#B CELLS
P1 <- FeaturePlot(seurat_phase,features = c("SMC6","CD79B","CD79A","MS4A1","JCHAIN","MZB1","TXNDC5","POU2AF1","DERL3"),reduction = 'umap')
ggsave(filename = 'B cells-maker-features.pdf', plot=P1, width = 15, height = 15)
#Ⅱ myofibers
P1 <- FeaturePlot(seurat_phase,features = c("TNNT3","MYBPC2","TTLL11","ACTN3","GADL1","ATP2A1"),reduction = 'umap')
ggsave(filename = 'II myofiber-maker-features.pdf', plot=P1, width = 15, height = 15)
#Ⅰ myofibers
P1 <- FeaturePlot(seurat_phase,features = c("ATP2A2","MYH7","ARPP21","TECRL","ESRRG","MYOZ2","TNNT1","MYH7B","ABHD18","MYL3"),reduction = 'umap')
ggsave(filename = 'I myofiber-maker-features.pdf', plot=P1, width = 15, height = 15)
#Fibroblasts1
P1 <- FeaturePlot(seurat_phase,features = c("ABCA6","ABCA8","SPARCL1","KAZN","HMCN2","DKK2","NID1","BICC1","COL6A3","FBN1"),reduction = 'umap')
ggsave(filename = 'Fibroblasts1-maker-features.pdf', plot=P1, width = 15, height = 15)
#Fibroblasts2
P1 <- FeaturePlot(seurat_phase,features = c("SLIT3","FGF2","HMCN1","PID1","DCN","PDGFRB"),reduction = 'umap')
ggsave(filename = 'Fibroblasts2-maker-features.pdf', plot=P1, width = 15, height = 15)
#Fibroblasts3
P1 <- FeaturePlot(seurat_phase,features = c("ACTA2","CACNB2","FN1","GLI3","GRK5","CDON"),reduction = 'umap')
ggsave(filename = 'Fibroblasts3-maker-features.pdf', plot=P1, width = 15, height = 15)
#Fibroblasts4
P1 <- FeaturePlot(seurat_phase,features = c("COL1A1","NEGR1","COL1A2","COL3A1"),reduction = 'umap')
ggsave(filename = 'Fibroblasts4-maker-features.pdf', plot=P1, width = 15, height = 15)
#Proximal tubule cells
P1 <- FeaturePlot(seurat_phase,features = c("ALDOB","GPX3","ASS1","PCK1","BHMT","FBP1","SLC13A3","SLC13A1","DCXR","SLC34A1","LRP2"),reduction = 'umap')
ggsave(filename = 'Proximal tubule cells-maker-features.pdf', plot=P1, width = 15, height = 15)
#Distal tubule cells 
P1 <- FeaturePlot(seurat_phase,features = c("UMOD","DEFB1","CA2","WFDC2","ATP6V0D2","ATP6V1C2","ATP6V1G3","TMEM213"),reduction = 'umap')
ggsave(filename = 'Distal tubule cells-maker-features.pdf', plot=P1, width = 15, height = 15)
#Endothelial_cells 
P1 <- FeaturePlot(seurat_phase,features = c("MECOM","ADAMTS9","PTPRB","PREX2","CYYR1","CDH5","PECAM1","ERG","VWF","KDR","CLEC14A","ID1","ENG","PROX1"),reduction = 'umap')
ggsave(filename = 'Endothelial-maker-features.pdf', plot=P1, width = 15, height = 15)
# Macrophage 
P1 <- FeaturePlot(seurat_phase,features = c("RBPJ","FRMD4B","MTSS1","CD163","MRC1","CD14","CSF1R","PYCARD","CD86","LYZ"),reduction = 'umap')
ggsave(filename = 'Macrophage-maker-features.pdf', plot=P1, width = 15, height = 15)
#Adipocytes
P1 <- FeaturePlot(seurat_phase,features = c("ADIPOR2","PLIN1","ADIPOQ","ADAMTS12","PLIN4","LIPE"),reduction = 'umap')
ggsave(filename = 'Adipocytes-maker-features.pdf', plot=P1, width = 15, height = 15)
#Smooth muscle like cells
P1 <- FeaturePlot(seurat_phase,features = c("SYNPO2","PDE3A","P3H2","DAAM2"),reduction = 'umap')
ggsave(filename = 'Smooth muscle like cells-maker-features.pdf', plot=P1, width = 15, height = 15)
#total
p1<-FeaturePlot(seurat_phase, features = c("GNLY","CRTAM","IL7R","ABCA6","TNNT3",
                                     "ALDOB","RBPJ","PECAM1","KLRB1",
                                     "FN1","SMC6","TMEM213","ADIPOR2",
                                     "JCHAIN","AHSP","SYNPO2"),cols = c("gray", "red"), pt.size = 1,combine = TRUE)+  
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
ggsave(filename = 'allcelltype-maker-features.pdf', plot=p1, width = 12, height = 10)


p1<-FeaturePlot(seurat_phase, features = c("SLBP", "NASP", "TUBB4B","NUSAP1"),cols = c("gray", "red"), pt.size = 1,combine = TRUE)+  
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))
ggsave(filename = 'yunshiS-G2M-maker-features.pdf', plot=p1, width = 12, height = 10)
