#####FOXP1-R514H PIPseq analysis#####

library(dplyr)
library(Seurat)
library(patchwork)
#library(EnhancedVolcano)
library(RColorBrewer)

setwd("/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2")

#load data
het3_1.data <- Read10X(data.dir = "/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/PIPseeker/het3-1_outs/filtered_matrix/sensitivity_4")
het3_2.data <- Read10X(data.dir = "/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/PIPseeker/het3-2_outs/filtered_matrix/sensitivity_4")
wt7_1.data <- Read10X(data.dir = "/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/PIPseeker/wt7-1_outs/filtered_matrix/sensitivity_5")
wt7_2.data <- Read10X(data.dir = "/media/chang/HDD-3/kelsey/FOXP1-R514H-pipseq/run2/PIPseeker/wt7-2_outs/filtered_matrix/sensitivity_4")

#initialize seurat object with non-normalized data
het3_1 <- CreateSeuratObject(counts = het3_1.data, project = "het3-1", min.cells = 3, min.features = 200)
het3_2 <- CreateSeuratObject(counts = het3_2.data, project = "het3-2", min.cells = 3, min.features = 200)
wt7_1 <- CreateSeuratObject(counts = wt7_1.data, project = "wt7-1", min.cells = 3, min.features = 200)
wt7_2 <- CreateSeuratObject(counts = wt7_2.data, project = "wt7-2", min.cells = 3, min.features = 200)

#combine objects#
data <- merge(het3_1, y = c(het3_2, wt7_1, wt7_2), add.cell.ids = c("het3-1", "het3-2", "wt7-1", "wt7-2"), project = "FOXP1-R514H")
data

#QC#
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
data <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 10 & nCount_RNA >= 1000)

#normalize data#
#data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
#find variable features#
#data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(data), 10)
top10

#scale data#
all.genes <- rownames(data)
#data <- ScaleData(data, features = all.genes)
#SCT replaces normalize, scale, and find variable features
data <- SCTransform(data, vars.to.regress = "percent.mt", verbose = FALSE)
#PCA#
data <- RunPCA(data, features = VariableFeatures(object = data))


#harmony#
data <- RunHarmony(data, group.by.vars = "orig.ident", assay.use = "SCT", plot_convergence = T)

#cluster cells#
data <- FindNeighbors(data, dims = 1:10)
data <- FindClusters(data, resolution = 0.2)
#umap#
data <- RunUMAP(data, dims = 1:10, reduction = "harmony")
p1 <- DimPlot(data, reduction = "umap")
p2 <- DimPlot(data, reduction = "umap", group.by = "orig.ident")
p1 +p2
saveRDS(data, file = "FOXP1-0.2res.rds")

data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(data.markers,"allmarkers-3-7-0.2res.csv")

#0.2 resolution clusters#
#clsuter 0 - IN
#cluster 1 - low quality
#cluster 2 - dividing
#cluster 3 - EN-1
#cluster 4 - EN-2
#cluster 5 - RG
#cluster 6 - EN-3

data2 <- subset(data, idents = c("0", "2", "3", "4", "5", "6"))
data2 <- SCTransform(data2, vars.to.regress = "percent.mt", verbose = FALSE)
data2 <- RunPCA(data2, features = VariableFeatures(object = data2))
data2 <- RunHarmony(data2, group.by.vars = "orig.ident", assay.use = "SCT",
                   plot_convergence = T)
data2 <- FindNeighbors(data2, dims = 1:10)
data2 <- FindClusters(data2, resolution = 0.2)
#umap#
data2 <- RunUMAP(data2, dims = 1:10, reduction='harmony')
p3 <- DimPlot(data2, reduction = "umap")
p4 <- DimPlot(data2, reduction = "umap", group.by = "orig.ident")
p3 +p4
saveRDS(data2, file = "FOXP1-3-7-0.2res-subset.rds")

data.markers <- FindAllMarkers(data2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(data.markers,"allmarkers-subset-3-7-0.2res.csv")

data.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
heatmap <- DoHeatmap(data2, features = top10$gene, label = F) + NoLegend()
heatmap

#dotplot#
Idents(data2) <- "cell_type"
genes <- c("NES", "MKI67", "DLX2", "CUX2", "BCL11B", "TBR1")
DotPlot(data2, features = genes)

##DEG##
data2@meta.data$condition[data2@meta.data$orig.ident == "het3-1"] = "R514H"
data2@meta.data$condition[data2@meta.data$orig.ident == "het3-2"] = "R514H"
data2@meta.data$condition[data2@meta.data$orig.ident == "wt7-1"] = "control"
data2@meta.data$condition[data2@meta.data$orig.ident == "wt7-2"] = "control"

het.deg <- FindMarkers(data2, ident.1 = "R514H", ident.2 = "control",
                        slot = 'data', group.by = 'condition')
write.csv(het.deg, "het_deg-3-7.csv")

Idents(data2) <- "orig.ident"
prop.table(table(Idents(data2)))
pt <- table(Idents(data2), data2$seurat_clusters)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

#####DEG by cluster#####
Idents(data2) <- "seurat_clusters"
clust0 <- subset(data2, idents = "0")
clust1 <- subset(data2, idents = "1")
clust2 <- subset(data2, idents = "2")
clust3 <- subset(data2, idents = "3")
clust4 <- subset(data2, idents = "4")
clust5 <- subset(data2, idents = "5")
clust6 <- subset(data2, idents = "6")

clust0_deg <- FindMarkers(clust0, ident.1 = "R514H", ident.2 = c("control"),
                          slot = 'data', group.by = 'condition')
clust1_deg <- FindMarkers(clust1, ident.1 = "R514H", ident.2 = c("control"),
                          slot = 'data', group.by = 'condition')
clust2_deg <- FindMarkers(clust2, ident.1 = "R514H", ident.2 = c("control"),
                          slot = 'data', group.by = 'condition')
clust3_deg <- FindMarkers(clust3, ident.1 = "R514H", ident.2 = c("control"),
                          slot = 'data', group.by = 'condition')
clust4_deg <- FindMarkers(clust4, ident.1 = "R514H", ident.2 = c("control"),
                          slot = 'data', group.by = 'condition')
clust5_deg <- FindMarkers(clust5, ident.1 = "R514H", ident.2 = c("control"),
                          slot = 'data', group.by = 'condition')
clust6_deg <- FindMarkers(clust6, ident.1 = "R514H", ident.2 = c("control"),
                          slot = 'data', group.by = 'condition')

write.csv(clust0_deg, file = "clust0-deg.csv")
write.csv(clust1_deg, file = "clust1-deg.csv")
write.csv(clust2_deg, file = "clust2-deg.csv")
write.csv(clust3_deg, file = "clust3-deg.csv")
write.csv(clust4_deg, file = "clust4-deg.csv")
write.csv(clust5_deg, file = "clust5-deg.csv")
write.csv(clust6_deg, file = "clust6-deg.csv")

#cluster 0 - dividing
#cluster 1 - IN
#cluster 2 - IN
#cluster 3 - EN/FOXP2
#cluster 4 - EN/FOXP1
#cluster 5 - RG
#cluster 6 - EN/subplate

p5 <- DimPlot(data2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
p5 + p4

data2@meta.data$cell_type[data2@meta.data$seurat_clusters %in% as.factor(c(0))] <- "Dividing"
data2@meta.data$cell_type[data2@meta.data$seurat_clusters %in% as.factor(c(1))] <- "IN-1"
data2@meta.data$cell_type[data2@meta.data$seurat_clusters %in% as.factor(c(2))] <- "IN-2"
data2@meta.data$cell_type[data2@meta.data$seurat_clusters %in% as.factor(c(3))] <- "EN-1"
data2@meta.data$cell_type[data2@meta.data$seurat_clusters %in% as.factor(c(4))] <- "EN-2"
data2@meta.data$cell_type[data2@meta.data$seurat_clusters %in% as.factor(c(5))] <- "Progenitors"
data2@meta.data$cell_type[data2@meta.data$seurat_clusters %in% as.factor(c(6))] <- "EN-3"

Idents(data2) <- "orig.ident"
prop.table(table(Idents(data2)))
pt <- table(Idents(data2), data2$cell_type)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

Idents(data2) <- "condition"
prop.table(table(Idents(data2)))
pt <- table(Idents(data2), data2$cell_type)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

saveRDS(data2, file = "FOXP1-R514H_seuratobj.rds")

EnhancedVolcano(clust6_deg, lab = clust6_deg$X, x = 'avg_log2FC', y = 'p_val_adj',
                pCutoff = 10e-6,
                FCcutoff = 0.2,
                title = 'EN-3: FOXP1 R514H vs. WT',
                pointSize = 3.0,
                legendLabSize = 15,
                legendIconSize = 5.0,
                legendPosition = 'right',
                maxoverlapsConnectors = Inf,
                drawConnectors = TRUE,
                xlim = c(-2, 2),
                selectLab = c("BCL11B", "NEUROD2", "SOX5", "TENM2"),
                labSize = 6.0)

het_fc <- FindMarkers(data2, ident.1 = "R514H", ident.2 = "control",
                    slot = 'data', group.by = 'condition', min.cells.group = 1, 
                    min.cells.feature = 1,
                    min.pct = 0,
                    logfc.threshold = 0,
                    only.pos = FALSE)
write.csv(het_fc, file = "allgenes_fc.csv")

clust6_fc <- FindMarkers(clust6, ident.1 = "R514H", ident.2 = "control",
                      slot = 'data', group.by = 'condition', min.cells.group = 1, 
                      min.cells.feature = 1,
                      min.pct = 0,
                      logfc.threshold = 0,
                      only.pos = FALSE)
write.csv(clust6_fc, file = "clust6_fc.csv")
clust4_fc <- FindMarkers(clust4, ident.1 = "R514H", ident.2 = "control",
                         slot = 'data', group.by = 'condition', min.cells.group = 1, 
                         min.cells.feature = 1,
                         min.pct = 0,
                         logfc.threshold = 0,
                         only.pos = FALSE)
write.csv(clust4_fc, file = "clust4_fc.csv")

#####analyze coexpression of FOXP1, FOXP2, FOXP4#####
a<-WhichCells(data2, expression = FOXP1 > 0)
data2$FOXP1 <- ifelse(colnames(data2) %in% a, "FOXP1_pos", "NA")
#FOXP1_pos = 8268 cells, 33.5%

x<-WhichCells(data2, expression = FOXP1 > 0 & FOXP2 > 0)
data2$FOXP2_coexp <- ifelse(colnames(data2) %in% x, "FOXP1_FOXP2_pos", "NA")
px <- DimPlot(data2, reduction = "umap", group.by = "FOXP2_coexp") + theme(legend.position = 'bottom')
px
count_x <- as.data.frame(x = table(data2$FOXP2_coexp))
#FOXP1_FOXP2_pos = 2095 cells, 8.5%
y <- WhichCells(data2, expression = FOXP1 > 0 & FOXP4 > 0)
data2$FOXP4_coexp <- ifelse(colnames(data2) %in% y, "FOXP1_FOXP4_pos", "NA")
py <- DimPlot(data2, reduction = "umap", group.by = "FOXP4_coexp") + theme(legend.position = 'bottom')
py
count_y <- as.data.frame(x = table(data2$FOXP4_coexp))
#FOXP1_FOXP4_pos = 510 cells, 2.1%
z <- WhichCells(data2, expression = FOXP1 > 0 & FOXP2 > 0 & FOXP4 > 0)
data2$FOXP2_FOXP4_coexp <- ifelse(colnames(data2) %in% z, "FOXP1_FOXP2_FOXP4_pos", "NA")
pz <- DimPlot(data2, reduction = "umap",  group.by = "FOXP2_FOXP4_coexp") + theme(legend.position = 'bottom')
pz
count_z <- as.data.frame(x = table(data2$FOXP2_FOXP4_coexp))
#FOXP1_FOXP2_FOXP4_pos = 305 cells, 1.2%

px + py + pz

#####create a data frame for bar graph#####
df <- data.frame(coexpression = c("FOXP1_FOXP2_pos", "FOXP1_FOXP4_pos", "FOXP1_FOXP2_FOXP4_pos"),
                 percent_positive_total = c(8.5, 2.1, 1.2))
df$coexpression <- factor(df$coexpression, levels=c("FOXP1_FOXP2_pos", "FOXP1_FOXP4_pos", "FOXP1_FOXP2_FOXP4_pos"))
ggplot(df, aes(x=coexpression, y=percent_positive_total, fill = coexpression)) + geom_bar(stat = "identity") + 
  theme(legend.position = "none") +  theme(text = element_text(size = 24)) + 
  ggtitle("FOXP1 FOXP2 FOXP4 Coexpression") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + xlab("")

#comparing coexpression against total FOXP1+ cells, instead of all cells
df2 <- data.frame(coexpression = c("FOXP1_FOXP2_pos", "FOXP1_FOXP4_pos", "FOXP1_FOXP2_FOXP4_pos"),
                  percent_positive_FOXP1 = c(25.3, 6.2, 3.7))
df2$coexpression <- factor(df2$coexpression, levels=c("FOXP1_FOXP2_pos", "FOXP1_FOXP4_pos", "FOXP1_FOXP2_FOXP4_pos"))
ggplot(df2, aes(x=coexpression, y=percent_positive_FOXP1, fill = coexpression)) + geom_bar(stat = "identity") + 
  theme(legend.position = "none") +  theme(text = element_text(size = 24)) + 
  ggtitle("FOXP1 FOXP2 FOXP4 Coexpression") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + xlab("")

#get coexpression data based on condition#
#total control cells = 11930, R514H = 12733
table(data2$FOXP1, data2$condition)
#FOXP1_pos control cells = 4163, R514H = 4105
table(data2$FOXP2_coexp, data2$condition)
#FOXP1_FOXP2_coexp control = 1075 (1075/4163 = 25.8%), R514H = 1020 (1020/4105 = 24.8%)
#NA control = 10855, R514H= 11713
table(data2$FOXP4_coexp, data2$condition)
#FOXP1_FOXP4_coexp control = 272 (272/4163 = 6.5%), R514H = 238 (238/4105 = 5.8%)
#NA control - 11658, R514H 12495
table(data2$FOXP2_FOXP4_coexp, data2$condition)
#FOXP1_FOXP2_FOXP4_coexp control = 171 (171/4163 = 4.1%), R514H = 134 (134/4105 = 3.3%)
#NA control = 11759, 12599

df3 <- data.frame(coexpression = c("FOXP1_FOXP2_pos", "FOXP1_FOXP4_pos", "FOXP1_FOXP2_FOXP4_pos", "FOXP1_FOXP2_pos", "FOXP1_FOXP4_pos", "FOXP1_FOXP2_FOXP4_pos"),
                  condition = c("control", "control", "control", "R514H", "R514H", "R514H"),
                  percent_positive_FOXP1 = c(25.8, 6.5, 4.1, 24.8, 5.8, 3.3))
df3$coexpression <- factor(df3$coexpression, levels=c("FOXP1_FOXP2_pos", "FOXP1_FOXP4_pos", "FOXP1_FOXP2_FOXP4_pos"))
ggplot(df3, aes(x=coexpression, y=percent_positive_FOXP1, fill = condition)) + geom_bar(stat = "identity", position = "dodge") + 
  theme(legend.position = "right") +  theme(text = element_text(size = 24)) + 
  ggtitle("FOXP1 FOXP2 FOXP4 Coexpression") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + xlab("")

#####analysis of cluster6 subset, recluster, find DEGs#####
clust6 <- SCTransform(clust6, vars.to.regress = "percent.mt", verbose = FALSE)
clust6 <- RunPCA(clust6, features = VariableFeatures(object = clust6))
#clust_6 <- RunHarmony(data2, group.by.vars = "orig.ident", assay.use = "SCT",
 #                   plot_convergence = T)
clust6 <- FindNeighbors(clust6, dims = 1:10)
clust6 <- FindClusters(clust6, resolution = 0.3)
#umap#
clust6 <- RunUMAP(clust6, dims = 1:10, reduction='harmony')
p8 <- DimPlot(clust6, reduction = "umap")
p9 <- DimPlot(clust6, reduction = "umap", group.by = "orig.ident")
p10 <- DimPlot(clust6, reduction = "umap", group.by = "condition")
p8 + p9
p8 + p10

data.markers <- FindAllMarkers(clust6, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(data.markers,"no-wt4/clust6_subset/clust6_subset_allmarkers.csv")

Idents(clust6) <- "condition"
clust6_copy <- clust6
Idents(clust6_copy) <- factor(Idents(clust6_copy), levels= c('control', 'R514H'))
VlnPlot(clust6_copy, features = c("FOXP1", "FOXP2", "FOXP4"), pt.size = 0)

Idents(clust6) <- "condition"
prop.table(table(Idents(clust6)))
pt <- table(Idents(clust6), clust6$seurat_clusters)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

#rename clusters as broad cluster 0 and 1#
clust6@meta.data$broad_clusters[clust6@meta.data$seurat_clusters %in% as.factor(c(0))] <- "0"
clust6@meta.data$broad_clusters[clust6@meta.data$seurat_clusters %in% as.factor(c(1))] <- "1"
clust6@meta.data$broad_clusters[clust6@meta.data$seurat_clusters %in% as.factor(c(2))] <- "1"
clust6@meta.data$broad_clusters[clust6@meta.data$seurat_clusters %in% as.factor(c(3))] <- "0"

p11 <- DimPlot(clust6, reduction = "umap", group.by = "broad_clusters")
p11 + p10

Idents(clust6) <- "condition"
prop.table(table(Idents(clust6)))
pt <- table(Idents(clust6), clust6$broad_clusters)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme(legend.title = element_blank())

#DEG between clusters 0 and 1#
subcluster_deg <- FindMarkers(clust6, ident.1 = "0", ident.2 = c("1"),
                          slot = 'data', group.by = 'broad_clusters')
write.csv(subcluster_deg, file = "no-wt4/clust6_subset/clust0_1_DEG.csv")

saveRDS(clust6, file = "no-wt4/clust6_subset/FOXP1-R514H_clust6subset.rds")


#####scatterplot of diff peaks and DEGs from cluster 6#####
scatter <- read.csv("CUTTag_DEG_overlap.csv")
rownames(scatter) <- scatter$hgnc_symbol
p <- ggplot(scatter, aes(x=log2fc_binding, y=log2fc_expression)) + geom_point() + geom_smooth(method = lm) + 
  geom_label(label=rownames(scatter), aes(size=14), vjust = 1.5) + theme(axis.title=element_text(size=20), 
                                                                         plot.title = element_text(size=26), legend.position = "none") + ggtitle("Correlation of FOXP1 differential peaks\nand gene expression")
p
