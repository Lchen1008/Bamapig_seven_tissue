#安装marker基因统计高效程序包
install.packages('BiocManager')
BiocManager::install('limma')
BiocManager::install("DESeq2")
BiocManager::install("dplyr")
#安装其他需要的软件包
install.packages('XML')
install.packages('locfit')

#加载一些分析需要的软件包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(limma)
library(DESeq2)
library(scater)
library(SingleR)
library(scRNAseq)
#首先R工作目录改变到10X文件所在目录
dir()
data_dir <- 'filtered_feature_bc_matrix'
#导入count matrix data
list.files(data_dir)
data <- Read10X(data.dir = data_dir)
#将data转化成seurat格式
liver <- CreateSeuratObject(counts = data, project = "liver", min.cells = 1)
#注意这里的 min.cells =1 参数使用也是一种质控条件，意思是确保feature至少在1个cell里面表达，去掉了那些counts数为0的feature
#计算线粒体基因占比，13个线粒体基因每个都算
liver[["ND6"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ND6")
liver[["ND5"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ND5")
liver[["ND4L"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ND4L")
liver[["ND4"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ND4")
liver[["ND3"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ND3")
liver[["ND2"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ND2")
liver[["ND1"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ND1")
liver[["CYTB"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "CYTB")
liver[["COX3"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "COX3")
liver[["COX2"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "COX2")
liver[["COX1"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "COX1")
liver[["ATP8"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ATP8")
liver[["ATP6"]] <- PercentageFeatureSet(liver, assay = "RNA", feature = "ATP6")
#计算并加入线粒体基因总占比
head(liver@meta.data, 2)    #查看metadata里面前5行数据
liver[["percent.mt"]] <- apply(liver@meta.data[,5:17], 1, sum)  #求和得到总的线粒体占比
liver@meta.data <- liver@meta.data[,-c(5:17)]  #删除中间数据
#
#计算核糖体基因占比
#查看RPS和RPL开头的基因名称
rb.genes <- rownames(liver)[grep("^RP[SL]",rownames(liver))]
#
liver[["percent.rb"]] <- PercentageFeatureSet(liver, pattern = "^RP[SL]")
head(liver@meta.data, 3)
#导出原始的cell数量信息，可用excel打开手动统计平均值
write.csv(liver@meta.data, file = "liver-cells-metadata-质控前的原始数量信息.csv")
#
#导出原始的feature counts表（很大一般不做这一步）
##counts <- as.data.frame(as.matrix(liver@assays$RNA@counts))
##write.csv(counts, file="liver-原始基因counts.csv")
#保存原始数据
saveRDS(liver, file = "lung-rawdata.rds")
?VlnPlot
p1 <- VlnPlot(liver, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb" ), ncol = 4, pt.size = 0.1) + NoLegend() 
ggsave("liver-原始数据-小提琴图.pdf", plot=p1, width = 20, height =8)
#单独作图并调整Y轴范围便于观察 y.max = 25000
p1 <- VlnPlot(liver, features ="nCount_RNA", pt.size = 0.1, y.max = 25000) + NoLegend()
ggsave("liver-原始counts-小提琴图.pdf", plot=p1, width = 8, height =8)
#单独做gene， 
p1 <- VlnPlot(liver, features ="nFeature_RNA", pt.size = 0.1, y.max = 3000) + NoLegend()
ggsave("liver-原始gene数-小提琴图.pdf", plot=p1, width = 8, height =8)
#单独做线粒体占比图 
p1 <- VlnPlot(liver, features ="percent.mt", pt.size = 0.1, y.max = 25) + NoLegend()
ggsave("liver-原始线粒体占比-小提琴图.pdf", plot=p1, width = 8, height =8)
#基因数与counts数散点图
p2 <- FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("liver-原始数据counts_gene-散点图.pdf", plot=p2, width = 8, height =8)
#要求单个基因至少在10个cells中有表达
min.cells <- 10
num.cells <- rowSums(as.matrix(liver@assays$RNA@counts) > 0)  #有表达的基因对应的细胞数求和，返回矩阵
genes.use <- names(num.cells[which(num.cells >= min.cells)])  #判断获取表达细胞数大于等于10个的基因ID 
mykeepgene <- c(1:nrow(liver))[rownames(liver)%in%as.character(genes.use)] #
liver_f0 <- subset(liver, features=mykeepgene)   #subset函数十分重要，需要仔细研究参数用法
liver_f0
#过滤掉血液的污染基因 (这里参考的是一篇cell文章中的42个污染基因）
filter_genes <-c("MALAT1", "SLC4A1", "KDM5D", "ANK1", "DDX3Y", "EIF2AK1", "HBQ1", "FTL", "GATA1", "KLF1", "USP9Y", "NFE2", "MT1G", "RPS4Y1", "HBZ", "GYPC", "HEMGN", "SLC25A37", "ALAS2", "EPB41", "AHSP", "GYPA", "UTY", "HBA2", "HBG2", "EIF1AY", "HBA1", "HBM", "HBE1", "HBG1", "MTRNR2L4", "HBB", "MTRNR2L5", "MTRNR2L8", "MTRNR2L10", "MTRNR2L3", "MTRNR2L1", "MTRNR2L7", "MTRNR2L12", "MTRNR2L11", "MTRNR2L13", "MTRNR2L6")
#
keepgenes <- c(1:nrow(liver_f0))[!(rownames(liver_f0)%in%as.character(filter_genes))]
liver_f1 <- subset(liver_f0,features=keepgenes);liver_f1
tmp1 <- NormalizeData(ULB_f1)
tmp1 <- FindVariableFeatures(tmp1, selection.method="vst", nfeatures=3000)
all.genes <- rownames(ULB_f1)
tmp1 <- ScaleData(tmp1, features = all.genes)
tmp1 <- RunPCA(tmp1, features = VariableFeatures(tmp1))
tmp1 <- RunUMAP(tmp1, dims = 1:20)
tmp1 <- RunTSNE(tmp1, dims = 1:20)
#plot
p1 <- FeaturePlot(tmp1, features=c("nCount_RNA", "nFeature_RNA"), reduction ="tsne")
ggsave("ULB-FeaturePlot-TSNE-Count-geneÈÈÍ¼.pdf", plot=p1, width = 10, height =5)
p2 <- FeaturePlot(tmp1, features=c("nCount_RNA", "nFeature_RNA"), reduction ="umap")
ggsave("ULB-FeaturePlot-UMAP-Count-GeneÈÈÍ¼.pdf", plot=p2, width = 5, height =5)

metadata <- liver@meta.data
p1 <- ggplot(metadata, aes(color=samInfo, x=nCount_RNA, fill=samInfo)) + geom_density(alpha=0.2)+
  theme_classic() + scale_x_log10()+ ylab("Kernel density of cells") + geom_vline(xintercept=1000)
ggsave("liver-原始counts-细胞分布曲线.pdf", plot=p1, width = 8, height =4)
#细胞数与gene数的分布曲线
p1 <- ggplot(metadata, aes(color=samInfo, x=nFeature_RNA, fill=samInfo)) + geom_density(alpha=0.2)+
  theme_classic() + scale_x_log10()+ ylab("Kernel density of cells") + geom_vline(xintercept=100)
ggsave("lung-原始gene_细胞分布曲线.pdf", plot=p1, width = 8, height =4)

#这里只保留基因数量在大于500个，counts数量大于2000个的cells
liver_f2 <- subset(liver_f1, subset = nFeature_RNA >500 & nCount_RNA >2000)
liver_f2
#过滤线粒体和核糖体基因占比高的cells
liver_f3 <- subset(liver_f2, subset = percent.mt<=20 & percent.rb<=30)
liver_f3
最终我们获得了 7311个高质量的liver 细胞数据。

#导出质控后的cells数量信息
write.csv(liver_f3@meta.data, file = "liver-cells-metadata-质控后最终.csv")
#
p1 <- VlnPlot(liver_f3, features = c("nCount_RNA", "nFeature_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) + NoLegend()
ggsave("liver-质控后-小提琴图.pdf", plot=p1, width = 20, height =8)
#保存R数据
saveRDS(liver_f3, file = "liver-clean-data.rds")
#liver_f3 <- readRDS(file="liver-clean-data.rds")

## https://satijalab.org/seurat/v3.2/cell_cycle_vignette.html
#先下载需要的文件: cell_cycle_vignette_files.zip 
#Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "F:/单细胞测序课题/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)
#
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
#We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes
g2m.genes
能够成功获取这2类基因的id list
#
# Cell-Cycle Scoring
?CellCycleScoring
liver_f3 <- CellCycleScoring(liver_f3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(liver_f3@meta.data, 3)
# Visualize the distribution of cell cycle markers across
p1 <- RidgePlot(liver_f3, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 4)
ggsave("liver-细胞周期marker-山岭图.pdf", plot=p1, width = 16, height =4)
#在对细胞进行细胞周期评分后，我想使用PCA确定细胞周期是否是我们数据集中变异的主要来源
#检查之前先对数据做个标准化使得不同测序深度的细胞可比
seurat_phase <- NormalizeData(liver_f3)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method="vst", nfeatures=3000)
seurat_phase <- ScaleData(seurat_phase, features=rownames(seurat_phase))
# Perform PCA
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase))
p1 <- DimPlot(seurat_phase, reduction = "pca", group.by= "Phase")
ggsave("liver-PCA-细胞周期影响.pdf", plot=p1, width = 8, height =8) 
#然后测试通过回归能否很好去除细胞周期影响（可忽略，不做这一步，后面SCT后再作图查看效果）
#很慢，多线程运行 
library(future)
plan('multiprocess', workers = 4)
seurat_phase <- ScaleData(seurat_phase, vars.to.regress = c("S.Score", "G2M.Score"), #features=rownames(seurat_phase))
                          plan('sequential');gc(reset = TRUE)
                          # Perform PCA
                          seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase))
                          p1 <- DimPlot(seurat_phase, reduction = "pca", group.by= "Phase")
                          ggsave("liver-PCA-细胞周期回归后.pdf", plot=p1, width = 8, height =8)
                          #有效果！
                          seurat_phase <- RunPCA(seurat_phase, features = c(s.genes, g2m.genes))
                          p1 <- DimPlot(seurat_phase, reduction = "pca", group.by= "Phase")
                          ggsave("liver-PCA-细胞周期回归后-markers.pdf", plot=p1, width = 16, height =8)
                          
#对数据进行负二项分布和中心化，并且在这一步我们回归掉细胞周期和线粒体基因的对数据影响，然后选取指定数量的variable.features（默认是n = 3000）到后续的分析
? SCTransform
#
liver_f3 <- SCTransform(liver_f3 , assay = "RNA", vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))
liver_f3
#导出
write.csv(liver_f3@meta.data, file = "liver-cells-metadata-SCT后的数量信息.csv")
p1 <- VlnPlot(liver_f3, features=c("nCount_SCT","nFeature_SCT","percent.mt"), ncol=3, pt.size=0.1, group.by="samInfo") + NoLegend()
ggsave("liver-SCT后counts_gene_线粒体占比-小提琴图.pdf", plot=p1, width = 8, height =8)
#导出SCT后的counts表（这很大可忽略）
#counts <- as.data.frame(as.matrix(liver_f3@assays$SCT@counts))
#write.csv(counts, file="liver-SCT后的gene_counts.csv")
#
#导出var.features counts
var.counts <- liver_f3@assays$SCT@counts[ rownames(counts) %in% liver_f3@assays$SCT@var.features, ]
write.csv(var.counts, file="liver-SCT后的var.gene_counts.csv")
#
#保存R数据
saveRDS(liver_f3, file = "lung-clean-data-SCT后.rds")

liver_f3 <- RunPCA(liver_f3)
#作碎石图选择dim值
p1 <- ElbowPlot(liver_f3, ndims=50, reduction="pca")  
ggsave("liver-PCA-ElbowPlot.pdf", plot=p1, width = 5, height =5) 
##DimHeatmap(liver_f3, dims = 1, cells = 2000, balanced = TRUE)
#
liver_f3 <- FindNeighbors(liver_f3, dims = 1:20)
liver_f3 <- RunUMAP(liver_f3, dims = 1:20)
liver_f3 <- RunTSNE(liver_f3, dims = 1:20)
#作图查看通过SCT后能否很好去除细胞周期影响(新增）
p1 <- DimPlot(liver_f3, reduction = "pca", group.by= "Phase")
ggsave("liver-PCA-细胞周期-SCT后-3个周期一起.pdf", plot=p1, width = 5, height =5)
p1 <- DimPlot(liver_f3, reduction="pca", group.by="Phase", split.by="Phase")
ggsave("liver-PCA-细胞周期-SCT后-3个周期分开.pdf", plot=p1, width = 15, height =5)

liver_f3 <- FindClusters(liver_f3, resolution = 0.5)
#
p1 <- DimPlot(liver_f3, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("liver-r0.5-dims20-UMAP.pdf", plot=p1, width = 6, height =5)
p2 <- DimPlot(liver_f3, reduction = "tsne", label = TRUE, repel = TRUE)
ggsave("liver-r0.5-dims20-UMAP.pdf", plot=p2, width = 6, height =5)

#作图查看每个cluster的基因数量和counts数量，研究哪些cluster是整体质量较差的类群
p1 <- VlnPlot(liver_f3, features="nCount_SCT", group.by="seurat_clusters", same.y.lims=TRUE, fill.by="seurat_clusters", pt.size = 0.1) + NoLegend()
ggsave("liver-cluster后-nCount_SCT-VlnPlot.pdf", plot=p1, width = 10, height =5)
p1 <- VlnPlot(liver_f3, features="nFeature_SCT", group.by="seurat_clusters", same.y.lims=TRUE, fill.by="seurat_clusters", pt.size = 0.1) + NoLegend()
ggsave("liver-cluster后-nFeature_SCT-VlnPlot.pdf", plot=p1, width = 10, height =5)

#UMAP-heatmap 从聚类图上展示细胞的基因和counts数量（另一种观察视野）
?FeaturePlot
p1 <- FeaturePlot(liver_f3, features="nCount_SCT", reduction ="umap", label=TRUE )
ggsave("liver-cluster后-nCount_SCT-FeaturePlot.pdf", plot=p1, width = 5, height =5)
p2 <- FeaturePlot(liver_f3, features="nFeature_SCT", reduction ="umap", label=TRUE )
ggsave("liver-cluster后-nFeature_SCT-FeaturePlot.pdf", plot=p2, width = 5, height =5)

#保存每一个cluster有的细胞数量
cluster <- table(muscle_f3@active.ident)
write.csv(cluster, file="spleen-cluster-r0.5-cell-number.csv")
#保存R数据
saveRDS(liver_f3, file = "liver-clean-data-cluster后.rds")
#读取
#liver_f3 <- readRDS(file = "liver-clean-data-cluster后.rds")
#liver_f3


#使用SingleR （预测）
Aran D, Looney AP, Liu L et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. _Nat. Immunology_ 20, 163–172.
#参考网站
http://www.bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
#安装一系列的依赖软件包
#install.packages(c("outliers", "pbmcapply", "doFuture"))
#BiocManager::install(c("GSVA","singscore"))
#BiocManager::install("SingleR") 
#
install.packages("scRNAseq")
BiocManager::install("scRNAseq")
BiocManager::install("scater")
BiocManager::install("SingleCellExperiment")
#
#加载
library(scater)
library(SingleR)
library(scRNAseq)
library(dplyr)

#加载自带的数据库
#we obtain reference data from the Human Primary Cell Atlas
#我们这里采用两个人类的参考集去做细胞注释，有两种：
#1   HPCA，人类原代细胞图谱(基于微阵列) 
#2  蓝图\_encode，一种结合蓝图、表基因组学和编码数据集(RNAseq)。
#
# The ENCODE Project Consortium (2012). An integrated encyclopedia of DNA elements in the human genome. Nature 489, pages 57–74.
# Martens JHA and Stunnenberg HG (2013). BLUEPRINT: mapping human blood cell epigenomes. Haematologica 98, 1487–1489.
# 数据格式方法学统一按照这篇文章的 Aran D, Looney AP, Liu L et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat. Immunol. 20, 163–172.
#
# Mabbott NA et al. (2013). An expression atlas of human primary cells: inference of gene function from coexpression 
networks. BMC Genomics 14, Article 632.
#
?HumanPrimaryCellAtlasData
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
?BlueprintEncodeData
bpe.se <- BlueprintEncodeData()
bpe.se
#小鼠参考数据库
?MouseRNAseqData
# Benayoun B et al. (2019). Remodeling of epigenome and transcriptome landscapes with aging in mice reveals widespread induction of inflammatory responses. Genome Res. 29, 697-709.
#
mbe.se <- MouseRNAseqData()
mbe.se
?ImmGenData
# Heng TS, Painter MW, Immunological Genome Project Consortium (2008). The Immunological Genome Project: networks of gene expression in immune cells. Nat. Immunol. 9, 1091-1094.
mig.se <- ImmGenData()
mig.se
#由于这些参考数据库都是即时下载生成的，我们考虑保存减少下一次重复下载时间?
saveRDS(hpca.se, file = "Human_PrimaryCellAtlas_SC_ref.rds")
saveRDS(bpe.se, file = "Human_BlueprintEncodeData_SC_ref.rds")
saveRDS(mbe.se, file = "Mouse_SC_ref.rds")
saveRDS(mig.se, file = "Mouse_immune_SC_ref.rds")
#
hpca.se <- readRDS(file = "F:/单细胞测序课题/cell_ref_database/Human_PrimaryCellAtlas_SC_ref.rds")
bpe.se <- readRDS(file = "F:/单细胞测序课题/cell_ref_database/Human_BlueprintEncodeData_SC_ref.rds")
mbe.se <- readRDS(file = "F:/单细胞测序课题/cell_ref_database/Mouse_SC_ref.rds")
mig.se <- readRDS(file = "F:/单细胞测序课题/cell_ref_database/Mouse_immune_SC_ref.rds")

#读取前面保存的seurat object文件
liver_f3 <- readRDS(file = "liver-clean-data-cluster后.rds")
liver_f3
#
#格式转换成SingleCellExperiment object
liver_f3.sce <- as.SingleCellExperiment(liver_f3)
#
#counts值log,因为ref数据库是用的logcounts所以需要匹配
liver_f3.sce <- logNormCounts(liver_f3.sce)
liver_f3.sce
#细胞注释，同时使用多个ref数据库
?SingleR
#注意参数较多需要学习
Anno <- SingleR(test=spleen_f3.sce, ref=list(HP=hpca.se, BP=bpe.se, MB=mbe.se, MI=mig.se), labels=list(hpca.se$label.main, bpe.se$label.main, mbe.se$label.main, mig.se$label.main), assay.type.test="logcounts", assay.type.ref="logcounts", method="cluster", cluster=spleen_f3.sce$seurat_clusters)
Anno
#前面使用method = "cluster", cluster = spleen_f3.sce$seurat_clusters参数使得对我们的18个cluster细胞组进行注释，#如果不用该参数的话，程序会对每个细胞单独进行注释输出结果 （不推荐，非常耗时！！）
#Anno2 <- SingleR(test=spleen_f3.sce, ref=list(HP=hpca.se, BP=bpe.se), labels=list(hpca.se$label.main, bpe.se$label.main), assay.type.test="logcounts", assay.type.ref="logcounts")
#Anno2
#导出
saveRDS(Anno, file = "lung-singleR-细胞cluster注释后.rds")
## saveRDS(Anno2, file = "liver-singleR-单个细胞注释后.rds")
## Anno2 <- readRDS(file = "liver-singleR注释后.rds")
#
#细胞类型注释结果导出成csv
celltype = data.frame(ClusterID=rownames(Anno), celltype= Anno$labels, ref=Anno$reference, score=Anno$scores, stringsAsFactors = F)
write.csv(celltype,"spleen-celltype_singleR-细胞cluster鉴定结果.csv",row.names = F)


#SingleR提供了强大的可视化工具。 plotScoreHeatmap()显示所有参考标签上的分数，这使用户可以检查整个数据集中预测标签的置信度。每个细胞的实际分配标签显示在顶部。理想情况下，每个cell（即热图的一列）应具有一个明显大于其余得分的分数，表明已将其明确分配给标签。
#
install.packages("pheatmap")
library(pheatmap)
?plotScoreHeatmap
plotScoreHeatmap(Anno, show_colnames=TRUE, show.labels=FALSE, fontsize_row=10, fountsize_col=10)
#手动另存为PDF
#将细胞注释信息重新添加到Seurat对象中去
#按照cluster分组注释的信息也添加进去 if `method="cluster"` was used:
head(liver_f3@meta.data, n=5)
#
liver_f3[["celltype.clusters"]] <- Anno$labels[match(liver_f3[[]][["seurat_clusters"]], rownames(Anno))]
head(liver_f3@meta.data, n=5)
#
#按照每个单个细胞注释结果逐个添加
liver_f3[["celltype.single"]] <- Anno2$labels[match(rownames(liver_f3@meta.data), rownames(Anno2))]
head(liver_f3@meta.data, n=5)
#其他把细胞注释信息添加回seurat的2种办法(备用）
#每个cluster手动灵活替换，重命名每个cluster为对应的细胞类型
##immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",  `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",  `10` = "Mk", `11` = "pDC", `12` = "Eryth")

#1 按照cluster分组进行注释的结果作图
p1 <- DimPlot(liver_f3, reduction = "umap", label = TRUE, group.by="celltype.clusters")
ggsave("liver-r0.5-dims20-UMAP-细胞按照cluster分组注释结果.pdf", plot=p1, width = 6, height =5)
#2 按照单个细胞进行注释后的结果作图
p1 <- DimPlot(liver_f3, reduction="umap", label=TRUE, repel=TRUE, group.by="celltype.single")
ggsave("liver-r0.5-dims20-UMAP-按照单个细胞注释结果.pdf", plot=p1, width = 10, height =6)

Seura提供多种差异分析的方法，默认wilcox方法，MASK方法和DESeq2方法也可以留意一下。MASK是专门针对单细胞数据差异分析设计的，DESeq2是传统bulkRNA数据差异分析的经典方法。还有area under the curve score
#
liver_f3@active.ident   #查看目前的cluster分组idents
?FindAllMarkers
#计算所有cluster的marker基因，非常耗时所以我们采用多核并行计算
library(future)
plan('multiprocess', workers = 3)
difgene1 <- FindAllMarkers(object=liver_f3, assay ="SCT", slot= "counts", min.pct=0.25, logfc.threshold=0.25, test.use="DESeq2", only.pos=TRUE)
plan('sequential');gc(reset = TRUE)
#从结果中过滤掉不显著的基因
difgene1 <- subset(difgene1, p_val_adj < 0.05)
head(difgene1, 5)
#导出
write.csv(difgene1, file="liver-r0.5-16cluster-marker基因.csv")
saveRDS(difgene1, file ="liver-marker-genes.rds")
Seura提供多种差异分析的方法，默认wilcox方法，MASK方法和DESeq2方法也可以留意一下。MASK是专门针对单细胞数据差异分析设计的，DESeq2是传统bulkRNA数据差异分析的经典方法。还有area under the curve score
#
liver_f3@active.ident   #查看目前的cluster分组idents
?FindAllMarkers
#计算所有cluster的marker基因，非常耗时所以我们采用多核并行计算
library(future)
plan('multiprocess', workers = 3)
difgene1 <- FindAllMarkers(object=liver_f3, assay ="SCT", slot= "counts", min.pct=0.25, logfc.threshold=0.25, test.use="DESeq2", only.pos=TRUE)
plan('sequential');gc(reset = TRUE)
#从结果中过滤掉不显著的基因
difgene1 <- subset(difgene1, p_val_adj < 0.05)
head(difgene1, 5)
#导出
write.csv(difgene1, file="liver-r0.5-16cluster-marker基因.csv")
saveRDS(difgene1, file ="liver-marker-genes.rds")

#DEGs分析
#为了进行个性化分组比较，我们将数据按照分组进行拆分
#由于findallmarker函数是根据ident来做检验的，所以需要设置ident符合我们指定的分组
levels(bulk2)  #查看目前的分组idents
bulk2@active.ident  #查看目前的分组idents
#将差异比较的ident设置为metadata中指定的分组
#head(bulk2@meta.data)
Idents(bulk2) <- bulk2@meta.data$group1
bulk2@active.ident
#
#确保ident正确后，将object按照metadata分组进行拆分,方便我们使用FindAllMarkers来个性化找DEGs
#比如现在我们按照SMT部位进行拆分，然后分别检验单个SMT部位内5个不同物种间的DEGs
tissue.split <- SplitObject(bulk2, split.by = "tissue")
tissue.split$SEM #查看
#Idents(tissue.split$SEM)  #确保ident正确
#保存
saveRDS(tissue.split, file = "tissue.split.rds")
#tissue.split <- readRDS(file="tissue.split.rds"); tissue.split
#同样按照物种拆分
species.split <- SplitObject(bulk2, split.by = "species")
saveRDS(species.split, file = "species.split.rds")
#species.split <- readRDS(file="species.split.rds"); species.split
#
#?FindAllMarkers
SEM <- FindAllMarkers(tissue.split$SEM, assay="RNA",logfc.threshold=1, min.pct=0.6, only.pos=TRUE)
#
#
write.csv(SEM, file="FindAllMarkers_SEM.csv")
#后台运行R脚本进行差异基因统计检验防止掉线
which Rscript
#
#将下面的命令写成R文件保存
#! /usr/bin/Rscript
library(future)
library(Seurat)
#
m7sct <- readRDS(file="pig-cleandata-m7sct-umap.rds");m7sct
m7sct <- FindClusters(m7sct, resolution=0.4)
#
?FindAllMarkers
plan("multiprocess", workers =10)
markers <- FindAllMarkers(object=m7sct, assay="SCT", slot= "counts", logfc.threshold=0.25, return.thresh =0.05, only.pos=TRUE)
#
write.csv(markers, file="FindAllMarkers-r0.4-31cluster-raw.csv")

#shell 执行
chmod +x test-findmarker.R
nohup Rscript ./test-findmarker.R &
  
  #回到R继续
  markers <- read.csv("FindAllMarkers-r0.4-31cluster-raw.csv", header=TRUE)
head(markers)
markers <- subset(markers, p_val_adj < 0.05)
markers <- markers[order(markers$avg_log2FC, decreasing=T),] #sort
markers <- markers[order(markers$cluster,decreasing=T),] #sort
#
#add ensembl gene ID into the marker gene results
#Load the homologous gene list between human and pigs
ortho <- read.csv("pig_to_human_orhtologous_add_MYH24.csv", header=TRUE)
head(ortho)
#load pig gene list
piggene <- read.csv("pig_ensembl102_gene_add_MYH24.csv", header=TRUE)
head(piggene)
#
markers$human_gene_id <- ortho$human_gene_id[match(markers$gene, ortho$pig_gene_name)]
markers$human_gene_name <- ortho$human_gene_name[match(markers$gene, ortho$pig_gene_name)]
markers$pig_gene_id <- piggene$pig_gene_id[match(markers$gene, piggene$pig_gene_name)]
head(markers)

#导出
write.csv(markers, file="markers-r0.4-31cluster.csv")



#单基因扫图
cd /media/disk3/zengbo/pig-SC-7sample
source ~/.bashrc
source activate /root/miniconda3/envs/Seurat
R

library(Seurat)
library(ggplot2)
m7sct <- readRDS(file="pig-cleandata-m7sct-umap.rds");m7sct
#
p1 <- FeaturePlot(m7sct, features="ADAMTS12", reduction ="umap")
ggsave("ADAMTS12.pdf", plot=p1, width = 6, height =5)
#


#单独找到cluster0所对应的marker（positive+negative）
?FindMarkers
#cluster0_markers <- FindMarkers(object= liver_f3, ident.1= 0, min.pct=0.25, logfc.threshold=0.25, #test.use="DESeq2", only.pos=TRUE)
#head(x = cluster0_markers, n = 5))
#write.csv(cluster0_markers, file="liver-r0.5-cluster0-差异基因.csv")
#找到区分cluster5和cluster0,3的marker
#cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#print(x = head(x = cluster5.markers, n = 5))

#针对给定的细胞和基因，用DoHeatmap函数画出表达热图（方便直观研究细胞类型对应的高表达marker基因）
#这里选定的是每个cluster中表达量前5的基因
difgene1
top5 <- difgene1 %>% group_by(cluster) %>% top_n(5, avg_logFC)
#
#分别用counts 和 scale.data作热图
p1 <- DoHeatmap(liver_f3, features=top5$gene, assay="SCT", slot="data")
p2 <- DoHeatmap(liver_f3, features=top5$gene, assay="SCT", slot="scale.data")

#先提取sct.variable基因（注意，这里不研究非variable的基因）
var.counts <- liver_f3@assays$SCT@counts[ rownames(liver_f3@assays$SCT@counts) %in% liver_f3@assays$SCT@var.features, ]
var.counts <- t(var.counts)  #转置一下
var.counts <- as.data.frame(var.counts) 
head(liver_f3@meta.data, n=5)
#
#将cluster分组列添加到counts表中
var.counts$seurat_clusters <- liver_f3@meta.data$seurat_clusters[ rownames(liver_f3@meta.data) %in% rownames(var.counts) ] 
head(var.counts, n=5)
#
#根据cluster 分组计算每个基因的平均表达量
library(dplyr)
cluster.mean <- var.counts %>% group_by(seurat_clusters) %>% summarize_each(funs(mean))
cluster.mean <- as.data.frame(cluster.mean)
#
cluster.mean2 <- t(cluster.mean)
#查看表达量高的前5个基因
cluster.mean2[ order(cluster.mean2[,1], decreasing=T) [1:5],]  #筛选表达量高的前5个基因
#
#导出平均表达量数据表
cluster.mean2 <- t(cluster.mean)   #head(cluster.mean, 10)
write.csv(cluster.mean2, file="liver-variable基因-cluster分组平均counts.csv")
saveRDS(cluster.mean2, file ="liver-variable基因-cluster分组平均counts.rds")

?top_n

#可视化个别gene的表达
VlnPlot(object = pbmc, features.plot = c("MS4A1", "CD79A"))
#可以选择用raw UMI counts
VlnPlot(object = pbmc, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)
#
#用PCA或者tSNE图可视化基因的表达
FeaturePlot(object = pbmc, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")



#基因表达气泡图，指定的各种基因在每种细胞类群中的表达量（待研究）
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("pDC", "Eryth", "Mk", "DC",     "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")
DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +    RotatedAxis()
# 山岭图 Visualize the distribution of cell cycle markers across （待研究）
RidgePlot(marrow, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

