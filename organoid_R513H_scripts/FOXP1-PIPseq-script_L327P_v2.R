#####FOXP1-R514H PIPseq analysis#####
#updating 20251202 with L327P data

library(dplyr)
library(Seurat)
library(patchwork)
#library(EnhancedVolcano)
library(RColorBrewer)
library(harmony)

setwd("/media/chang/HDD-11/kelsey/FOXP1/L327P/pipseq")

#load data
wt1.data <- Read10X(data.dir = "pipseeker/WT1/filtered_matrix/sensitivity_4")
wt2.data <- Read10X(data.dir = "pipseeker/WT2/filtered_matrix/sensitivity_4")
het15.data <- Read10X(data.dir = "pipseeker/HET15/filtered_matrix/sensitivity_5")
het18.data <- Read10X(data.dir = "pipseeker/HET18/filtered_matrix/sensitivity_4")

#initialize seurat object with non-normalized data
wt1 <- CreateSeuratObject(counts = wt1.data, project = "wt1", min.cells = 3, min.features = 200)
wt2 <- CreateSeuratObject(counts = wt2.data, project = "wt2", min.cells = 3, min.features = 200)
het15 <- CreateSeuratObject(counts = het15.data, project = "het15", min.cells = 3, min.features = 200)
het18 <- CreateSeuratObject(counts = het18.data, project = "het18", min.cells = 3, min.features = 200)

#######new code to try and manage differences in datasets that is causing problems later
#create l327p object
l327p <- merge(het15, y = c(wt1, wt2, het18), 
               add.cell.ids = c("het15", "wt1", "wt2", "het18"), project = "L327P")

#remove MT and low quality
l327p[["percent.mt"]] <- PercentageFeatureSet(l327p, pattern = "^MT-")
l327p <- subset(l327p, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 & nCount_RNA >= 1000)

#SCTransform and PCA
#following https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette
l327p <- SCTransform(l327p, vst.flavor = "v1", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident", assay.use = "SCT", plot_convergence = T) %>%
  RunUMAP(dims = 1:10, reduction = "harmony") %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5)
DefaultAssay(l327p) <- "RNA"
l327p <- NormalizeData(l327p, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(l327p)
l327p <- ScaleData(l327p, features = all.genes)


Idents(l327p) <- "seurat_clusters"
l327p <- subset(l327p, idents = "1", invert = TRUE) ## remove low quality

l327p <- SCTransform(l327p, vst.flavor = "v1", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>% 
  RunHarmony(group.by.vars = "orig.ident", assay.use = "SCT", plot_convergence = T) %>%
  RunUMAP(dims = 1:10, reduction = "harmony") %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5)
l327p <- NormalizeData(l327p, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(l327p)
l327p <- ScaleData(l327p, features = all.genes)

#some L327P analysis before integrating since it's being challenging
p1 <- DimPlot(l327p, reduction = "umap")
p2 <- DimPlot(l327p, reduction = "umap", group.by = "orig.ident")
p1 +p2
saveRDS(data, file = "L327P-0.5res.rds")

#find cluster markers
#PrepSCTFindMarkers(l327p, assay = "SCT")
data.markers <- FindAllMarkers(l327p, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "SCT")
data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(data.markers,"L327P_allmarkers-0.5res.csv")

Idents(l327p) <- "cell_type"
genes <- c("NES", "MKI67", "DLX2", "NEUROD2", "CUX2", "SATB2", "BCL11B", "TBR1", "TLE4", "ERBB4")
DotPlot(l327p, features = genes)

##DEG##
#genotype
l327p@meta.data$genotype[l327p@meta.data$orig.ident == "het15"] <- "L327P/WT"
l327p@meta.data$genotype[l327p@meta.data$orig.ident == "het18"] <- "L327P/WT"
l327p@meta.data$genotype[l327p@meta.data$orig.ident == "wt1"] <- "WT/WT"
l327p@meta.data$genotype[l327p@meta.data$orig.ident == "wt2"] <- "WT/WT"

Idents(l327p) <- "cell_type"
DefaultAssay(l327p) <- "SCT"
deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "DIV")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_DIV.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "EN-1")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-1.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "EN-2")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-2.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "EN-3")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-3.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "IN-1")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-1.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "IN-2")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-2.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "IN-3")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-3.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "IN-4")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-4.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "RG-1")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-1.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "RG-2")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-2.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "RG-3")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-3.csv")

deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "IPC")
write.csv(deg, "DEG/WT_vs_het/L327P-WT_vs_het_DEG_IPC.csv")

#make one file with all padj < 0.05 for all cell types for each comparison
df1 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_DIV.csv")
df2 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-1.csv")
df3 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-2.csv")
df4 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-3.csv")
df5 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-1.csv")
df6 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-2.csv")
df7 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-3.csv")
df8 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-4.csv")
df9 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_IPC.csv")
df10 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-1.csv")
df11 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-2.csv")
df12 <- read.csv("DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-3.csv")

df1 <- df1[df1$p_val_adj < 0.05,]
df1$cell_type <- "DIV"
df2 <- df2[df2$p_val_adj < 0.05,]
df2$cell_type <- "EN-1"
df3 <- df3[df3$p_val_adj < 0.05,]
df3$cell_type <- "EN-2"
df4 <- df4[df4$p_val_adj < 0.05,]
df4$cell_type <- "EN-3"
df5 <- df5[df5$p_val_adj < 0.05,]
df5$cell_type <- "IN-1"
df6 <- df6[df6$p_val_adj < 0.05,]
df6$cell_type <- "IN-2"
df7 <- df7[df7$p_val_adj < 0.05,]
df7$cell_type <- "IN-3"
df8 <- df8[df8$p_val_adj < 0.05,]
df8$cell_type <- "IN-4"
df9 <- df9[df9$p_val_adj < 0.05,]
df9$cell_type <- "IPC"
df10 <- df10[df10$p_val_adj < 0.05,]
df10$cell_type <- "RG-1"
df11 <- df11[df11$p_val_adj < 0.05,]
df11$cell_type <- "RG-2"
df12 <- df12[df12$p_val_adj < 0.05,]
df12$cell_type <- "RG-3"

df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12)
colnames(df)[1] <- "gene"

table(df$cell_type)
#DIV EN-1 EN-2 EN-3 IN-1 IN-2 IN-3 IN-4  IPC RG-1 RG-2 RG-3 
#322  282  282   32   55  132   27    3    7   68  258  122  

write.csv(df, "without_homozygous/DEG/WT_vs_het/WT_vs_het_allDEG_padj_0.05.csv")


##########transfer cell type labels from L327P onto R513H object###########
r513h <- readRDS("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/FOXP1-R514H_seuratobj.rds")
anchors <- FindTransferAnchors(reference = r513h, query = l327p, dims = 1:30,
                               reference.reduction = "harmony")
predictions <- TransferData(anchorset = anchors, refdata = r513h$cell_type, dims = 1:30)
l327p <- AddMetaData(l327p, metadata = predictions)

DimPlot(l327p, reduction = "umap", group.by = "predicted.id")


############full DEG analysis comparing those DEG in both L327P and R513H#############
#use predicted id from R513H mapped to L327P, annotate DEGs that are R513H biased or WT biased
#plot altogether on jjVolcano, highlight EN clusters of interest
#pearson correlation or FET?

table(l327p$cell_type, l327p$predicted.id)
#Dividing EN-1 EN-2 IN-1 IN-2   RG Subplate
#DIV       692    0    0    0    0   38        0
#EN-1       45   97  186    0    0   66     1013
#EN-2      361   39  111    1    3  134      715
#EN-3      348    0    3    0    0   17       43
#IN-1        1   21    0 1024  119    8        0
#IN-2      746    6    0   69  133   93        0
#IN-3      392    5    0   12  358    6        0
#IN-4       17  120   70    0    1    4        2
#IPC       245    0    0    0    0   22        0
#RG-1       66    0    0    0    0 1051        0
#RG-2       36    0    0    0    0  700        0
#RG-3      186    0    0    0    0  317        0

#DIV = DIV
#EN-1 = EN-3
#EN-2 = EN-3
#EN-3 = DIV
#IN-1 = IN-1
#IN-2 = DIV
#IN-3 = DIV, IN-2
#IN-4 = EN-1
#IPC = DIV
#RG-1 = RG
#RG-2 = RG
#RG-3 = RG

#one cell type at a time
df1 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_DIV.csv")
df2 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-1.csv")
df3 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-2.csv")
df4 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_EN-3.csv")
df5 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-1.csv")
df6 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-2.csv")
df7 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-3.csv")
df8 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_IN-4.csv")
df9 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_IPC.csv")
df10 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-1.csv")
df11 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-2.csv")
df12 <- read.csv("without_homozygous/DEG/WT_vs_het/L327P-WT_vs_het_DEG_RG-3.csv")

df1 <- df1[df1$p_val_adj < 0.05,]
df1$cell_type <- "DIV"
df2 <- df2[df2$p_val_adj < 0.05,]
df2$cell_type <- "EN-1"
df3 <- df3[df3$p_val_adj < 0.05,]
df3$cell_type <- "EN-2"
df4 <- df4[df4$p_val_adj < 0.05,]
df4$cell_type <- "EN-3"
df5 <- df5[df5$p_val_adj < 0.05,]
df5$cell_type <- "IN-1"
df6 <- df6[df6$p_val_adj < 0.05,]
df6$cell_type <- "IN-2"
df7 <- df7[df7$p_val_adj < 0.05,]
df7$cell_type <- "IN-3"
df8 <- df8[df8$p_val_adj < 0.05,]
df8$cell_type <- "IN-4"
df9 <- df9[df9$p_val_adj < 0.05,]
df9$cell_type <- "IPC"
df10 <- df10[df10$p_val_adj < 0.05,]
df10$cell_type <- "RG-1"
df11 <- df11[df11$p_val_adj < 0.05,]
df11$cell_type <- "RG-2"
df12 <- df12[df12$p_val_adj < 0.05,]
df12$cell_type <- "RG-3"

table(r513h$seurat_clusters, r513h$cell_type)
#0=DIV
#1=IN-1
#2=IN-2
#3=EN-1
#4=EN-2
#5=RG
#6=EN-3

#DIV
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust0-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df1$R513H_DEG[df1$X %in% r513h_mut$X] <- "R513H_enriched"
df1$R513H_DEG[df1$X %in% r513h_wt$X] <- "WT_enriched"

#EN-1
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust6-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df2$R513H_DEG[df2$X %in% r513h_mut$X] <- "R513H_enriched"
df2$R513H_DEG[df2$X %in% r513h_wt$X] <- "WT_enriched"

#EN-2
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust6-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df3$R513H_DEG[df3$X %in% r513h_mut$X] <- "R513H_enriched"
df3$R513H_DEG[df3$X %in% r513h_wt$X] <- "WT_enriched"

#EN-3
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust0-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df4$R513H_DEG[df4$X %in% r513h_mut$X] <- "R513H_enriched"
df4$R513H_DEG[df4$X %in% r513h_wt$X] <- "WT_enriched"

#IN-1
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust1-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df5$R513H_DEG[df5$X %in% r513h_mut$X] <- "R513H_enriched"
df5$R513H_DEG[df5$X %in% r513h_wt$X] <- "WT_enriched"

#IN-2
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust0-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df6$R513H_DEG[df6$X %in% r513h_mut$X] <- "R513H_enriched"
df6$R513H_DEG[df6$X %in% r513h_wt$X] <- "WT_enriched"

#IN-3
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust0-deg.csv")
deg2 <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust2-deg.csv")
r513h <- rbind(r513h_deg, deg2)
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df7$R513H_DEG[df7$X %in% r513h_mut$X] <- "R513H_enriched"
df7$R513H_DEG[df7$X %in% r513h_wt$X] <- "WT_enriched"

#IN-4
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust3-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df8$R513H_DEG[df8$X %in% r513h_mut$X] <- "R513H_enriched"
df8$R513H_DEG[df8$X %in% r513h_wt$X] <- "WT_enriched"

#IPC
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust0-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df9$R513H_DEG[df9$X %in% r513h_mut$X] <- "R513H_enriched"
df9$R513H_DEG[df9$X %in% r513h_wt$X] <- "WT_enriched"

#RG-1
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust5-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df10$R513H_DEG[df10$X %in% r513h_mut$X] <- "R513H_enriched"
df10$R513H_DEG[df10$X %in% r513h_wt$X] <- "WT_enriched"

#RG-2
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust5-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df11$R513H_DEG[df11$X %in% r513h_mut$X] <- "R513H_enriched"
df11$R513H_DEG[df11$X %in% r513h_wt$X] <- "WT_enriched"

#RG-3
r513h_deg <- read.csv("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/clust5-deg.csv")
r513h_mut <- subset(r513h_deg, r513h_deg$avg_log2FC > 0)
r513h_wt <- subset(r513h_deg, r513h_deg$avg_log2FC < 0)
df12$R513H_DEG[df12$X %in% r513h_mut$X] <- "R513H_enriched"
df12$R513H_DEG[df12$X %in% r513h_wt$X] <- "WT_enriched"

df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12)
colnames(df)[1] <- "gene"
write.csv(df, "without_homozygous/DEG/WT_vs_het/WT_vs_het_allDEG_padj_0.05_R513HDEG_anno.csv")


############change cell type anno for L327P based on R513H############
#DIV = DIV
#EN-1 = EN-3
#EN-2 = EN-3
#EN-3 = DIV
#IN-1 = IN-1
#IN-2 = DIV
#IN-3 = DIV, IN-2
#IN-4 = EN-1
#IPC = DIV
#RG-1 = RG
#RG-2 = RG
#RG-3 = RG


l327p@meta.data$cell_type_2[l327p$cell_type == "DIV"] <- "Dividing"
l327p@meta.data$cell_type_2[l327p$cell_type == "EN-1"] <- "EN-3"
l327p@meta.data$cell_type_2[l327p$cell_type == "EN-2"] <- "EN-2"
l327p@meta.data$cell_type_2[l327p$cell_type == "EN-3"] <- NA #using NAs for any cluster labeled as DIV that is not actually DIV, and unique to L327P
l327p@meta.data$cell_type_2[l327p$cell_type == "IN-1"] <- "IN-1"
l327p@meta.data$cell_type_2[l327p$cell_type == "IN-2"] <- NA
l327p@meta.data$cell_type_2[l327p$cell_type == "IN-3"] <- "IN-2"
l327p@meta.data$cell_type_2[l327p$cell_type == "IN-4"] <- "EN-1"
l327p@meta.data$cell_type_2[l327p$cell_type == "IPC"] <- NA
l327p@meta.data$cell_type_2[l327p$cell_type == "RG-1"] <- "Progenitors"
l327p@meta.data$cell_type_2[l327p$cell_type == "RG-2"] <- "Progenitors"
l327p@meta.data$cell_type_2[l327p$cell_type == "RG-3"] <- "Progenitors"

#for R513H
r513h <- readRDS("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/FOXP1-R514H_seuratobj.rds")
r513h$cell_type[r513h$cell_type == "Subplate"] <- "EN-3"
r513h$cell_type[r513h$cell_type == "RG"] <- "Progenitors"
r513h$condition[r513h$condition == "R514H"] <- "R513H"
saveRDS(r513h, "/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/FOXP1_R513H_object.rds")



############DEG using R513H label projection##############
Idents(l327p) <- "predicted.id"
deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "EN-3")
write.csv(deg, "DEG/R513H_labels/L327P-WT_vs_het_DEG_EN-3.csv")

l327p_df <- subset(deg, l327p_deg$p_val_adj < 0.05)
l327p_df$gene <- rownames(l327p_df)
l327p_df <- l327p_df[,c(2, 6)]
colnames(l327p_df)[1] <- "L327P_avg_log2FC"

deg_df <- merge(r513h_df, l327p_df, by.x = "gene", by.y = "gene", all = TRUE)

p <- ggscatter(deg_df, x = "R513H_avg_log2FC", y = "L327P_avg_log2FC",
               #color = "condition", fill = "condition",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "R513H_avg_log2FC", ylab = "L327P_avg_log2FC") + 
  ggtitle("EN-3 R513H L327P DEG correlation") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")

Idents(l327p) <- "cell_type_2"
deg <- FindMarkers(l327p, ident.1 = "L327P/WT", ident.2 = "WT/WT",
                   slot = 'data', group.by = 'genotype', subset.ident = "EN-3")
write.csv(deg, "DEG/R513H_labels/L327P-WT_vs_het_DEG_EN-3_cell_type_2.csv")

l327p_df <- subset(deg, l327p_deg$p_val_adj < 0.05)
l327p_df$gene <- rownames(l327p_df)
l327p_df <- l327p_df[,c(2, 6)]
colnames(l327p_df)[1] <- "L327P_avg_log2FC"

deg_df <- merge(r513h_df, l327p_df, by.x = "gene", by.y = "gene", all = TRUE)

p <- ggscatter(deg_df, x = "R513H_avg_log2FC", y = "L327P_avg_log2FC",
               #color = "condition", fill = "condition",
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "R513H_avg_log2FC", ylab = "L327P_avg_log2FC") + 
  ggtitle("EN-3 R513H L327P DEG correlation") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
  geom_vline(xintercept = 0, color = "black", linetype = "dashed")


########rename cell types in DEG list with R513H annotations#############
deg <- read.csv("/media/chang/HDD-11/kelsey/FOXP1/L327P/pipseq//DEG/WT_vs_het/WT_vs_het_allDEG_padj_0.05_R513HDEG_anno.csv")
deg$cell_type[deg$cell_type == "DIV"] <- "Dividing"
deg$cell_type[deg$cell_type == "EN-1"] <- "EN-3"
deg$cell_type[deg$cell_type == "EN-2"] <- "EN-2"
deg$cell_type[deg$cell_type == "EN-3"] <- NA
deg$cell_type[deg$cell_type == "IN-1"] <- "IN-1"
deg$cell_type[deg$cell_type == "IN-2"] <- NA
deg$cell_type[deg$cell_type == "IN-3"] <- "IN-2"
deg$cell_type[deg$cell_type == "IN-4"] <- "EN-1"
deg$cell_type[deg$cell_type == "IPC"] <- NA
deg$cell_type[deg$cell_type == "RG-1"] <- "Progenitors"
deg$cell_type[deg$cell_type == "RG-2"] <- "Progenitors"
deg$cell_type[deg$cell_type == "RG-3"] <- "Progenitors"
write.csv(deg, "/media/chang/HDD-11/kelsey/FOXP1/L327P/pipseq/DEG/WT_vs_het/WT_vs_het_allDEG_padj_0.05_R513HDEG_anno_R513H_celltype.csv")


