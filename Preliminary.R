setwd("/media/klp/Seagate Expansion Drive/elmallah/")

library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(sctransform)
library(ggplot2)
library(SoupX)
library(here)
library(pheatmap)
library(EnhancedVolcano)
library(slingshot)
library(RColorBrewer)
library(scales)

options(future.globals.maxSize= 8912896000000)
col_p = c("#F3E5C2","#EFB89A","#EF8F50", "#F5731D","#FA4602","#FE5100")


control.data3 <- Read10X("Control/outs/filtered_feature_bc_matrix/")


control.data2 <- load10X("Control/outs/")
control.data2 = autoEstCont(control.data2)
control.out = adjustCounts(control.data2)

DropletUtils:::write10xCounts("./Controlsoupx", control.out)

control.data <- Read10X("./Controlsoupx")

control <- CreateSeuratObject(counts = control.data,min.cells = 1)
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(control, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(control, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



control <- subset(control, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & nCount_RNA < 20000 & percent.mt <10)


control <- NormalizeData(control, normalization.method = "LogNormalize", scale.factor = 10000)
control <- FindVariableFeatures(control, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(control)
control <- ScaleData(control, features = all.genes)
control <- RunPCA(object = control, verbose = FALSE)
ElbowPlot(control, ndims = 30)
control <- FindNeighbors(object = control, dims = 1:20)
control <- FindClusters(object = control, resolution = 0.5)
control <- RunUMAP(object = control, dims = 1:20,min.dist = 1, n.neighbors = 40, spread = 5, local.connectivity = 20)

DimPlot(control, reduction = "umap", label = TRUE)

FeaturePlot(control, features = c("Sftpc"),cols=col_p, min.cutoff = "q10", order = T)

control_m <- FindAllMarkers(control, only.pos = T)
write.table(control_m,"control_m.txt",quote = F,sep='\t')



gaa.data3 <- Read10X("GAA_KO/outs/filtered_feature_bc_matrix/")


gaa.data2 <- load10X("GAA_KO/outs/")
gaa.data2 = autoEstCont(gaa.data2)
gaa.out = adjustCounts(gaa.data2)

DropletUtils:::write10xCounts("./GAA_KDsoupx", gaa.out)

gaa.data <- Read10X("./GAA_KDsoupx")


gaa <- CreateSeuratObject(counts = gaa.data,min.cells = 1)
gaa[["percent.mt"]] <- PercentageFeatureSet(gaa, pattern = "^mt-")
VlnPlot(gaa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(gaa, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(gaa, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(gaa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



gaa <- subset(gaa, subset = nFeature_RNA > 1200 & nFeature_RNA < 7500 & nCount_RNA < 25000 & percent.mt <10)


gaa <- NormalizeData(gaa, normalization.method = "LogNormalize", scale.factor = 10000)
gaa <- FindVariableFeatures(gaa, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gaa)
gaa <- ScaleData(gaa, features = all.genes)
gaa <- RunPCA(object = gaa, verbose = FALSE)
ElbowPlot(gaa, ndims = 30)
gaa <- FindNeighbors(object = gaa, dims = 1:22)
gaa <- FindClusters(object = gaa, resolution = 0.5)
gaa <- RunUMAP(object = gaa, dims = 1:22,min.dist = 1, n.neighbors = 40, spread = 5, local.connectivity = 20)

DimPlot(gaa, reduction = "umap", label = TRUE)

FeaturePlot(gaa, features = c("Sftpc"),cols=col_p, min.cutoff = "q10", order = T)

gaa_m <- FindAllMarkers(gaa, only.pos = T)
write.table(gaa_m,"gaa_m.txt",quote = F,sep='\t')


control$condition <- "Control"
gaa$condition <- "GAA KO"
gaa_merge <- merge(x=control,y=c(gaa))


gaa_merge <- NormalizeData(gaa_merge, normalization.method = "LogNormalize", scale.factor = 10000)
gaa_merge <- FindVariableFeatures(gaa_merge, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gaa_merge)
gaa_merge <- ScaleData(gaa_merge, features = all.genes)
gaa_merge <- RunPCA(object = gaa_merge, verbose = FALSE)
ElbowPlot(gaa_merge, ndims = 30)
gaa_merge <- FindNeighbors(object = gaa_merge, dims = 1:20)
gaa_merge <- FindClusters(object = gaa_merge, resolution = 1)
gaa_merge <- RunUMAP(object = gaa_merge, dims = 1:20,min.dist = 1, n.neighbors = 40, spread = 5, local.connectivity = 20)

DimPlot(gaa_merge, reduction = "umap", label = TRUE)
DimPlot(gaa_merge, reduction = "umap", label = F,group.by = "condition")


FeaturePlot(gaa_merge, features = c("Sftpc","Ager"),cols=col_p, min.cutoff = "q10", order = T)
FeaturePlot(gaa_merge, features = c("Pdgfra","Pecam1"),cols=col_p, min.cutoff = "q10", order = T)
FeaturePlot(gaa_merge, features = c("Col13a1","Col14a1"),cols=col_p, min.cutoff = "q10", order = T)


at2 <- subset(gaa_merge,idents = 9)
Idents(at2) <- at2$condition
at2_m <- FindMarkers(at2,ident.1 = "Control",ident.2 = "GAA KO",min.pct = 0.0,logfc.threshold = 0.01)

at1 <- subset(gaa_merge,idents = 20)
Idents(at1) <- at1$condition
at1_m <- FindMarkers(at1,ident.1 = "Control",ident.2 = "GAA KO",min.pct = 0.0,logfc.threshold = 0.01)

col13 <- subset(gaa_merge,idents = 0)
Idents(col13) <- col13$condition
col13_m <- FindMarkers(col13,ident.1 = "Control",ident.2 = "GAA KO",min.pct = 0.0,logfc.threshold = 0.01)

col14 <- subset(gaa_merge,idents = 3)
Idents(col14) <- col14$condition
col14_m <- FindMarkers(col14,ident.1 = "Control",ident.2 = "GAA KO",min.pct = 0.0,logfc.threshold = 0.01)

EnhancedVolcano(at2_m,
                title = "Control vs GAA KD",
                #  selectLab = m$V1,
                lab = rownames(at2_m),
                x = 'avg_logFC',
                y = 'p_val_adj',
                # ylim = c(-1,500),
                #xlim = c(-1.5,2),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-2,
                FCcutoff = 0.3,
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                legend=c('NS','log2 fold-change','P value',
                         'P value & log2 fold-change'),
                legendPosition = 'left',
                legendLabSize = 10,
                legendIconSize = 4.0,
                drawConnectors = T,widthConnectors = 0.75
)+NoLegend()


EnhancedVolcano(at1_m,
                title = "Control vs GAA KD",
                #  selectLab = m$V1,
                lab = rownames(at1_m),
                x = 'avg_logFC',
                y = 'p_val_adj',
                ylim = c(0,10),
                xlim = c(-1,1),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-2,
                FCcutoff = 0.35,
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                legend=c('NS','log2 fold-change','P value',
                         'P value & log2 fold-change'),
                legendPosition = 'left',
                legendLabSize = 10,
                legendIconSize = 4.0,
                drawConnectors = F,widthConnectors = 0.75
)+NoLegend()


EnhancedVolcano(col13_m,
                title = "Control vs GAA KD",
                #  selectLab = m$V1,
                lab = rownames(col13_m),
                x = 'avg_logFC',
                y = 'p_val_adj',
                #ylim = c(0,10),
                #xlim = c(-1,1),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-2,
                FCcutoff = 0.35,
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                legend=c('NS','log2 fold-change','P value',
                         'P value & log2 fold-change'),
                legendPosition = 'left',
                legendLabSize = 10,
                legendIconSize = 4.0,
                drawConnectors = F,widthConnectors = 0.75
)+NoLegend()

EnhancedVolcano(col14_m,
                title = "Control vs GAA KD",
                #  selectLab = m$V1,
                lab = rownames(col14_m),
                x = 'avg_logFC',
                y = 'p_val_adj',
                #ylim = c(0,10),
                #xlim = c(-1,1),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-2,
                FCcutoff = 0.35,
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                legend=c('NS','log2 fold-change','P value',
                         'P value & log2 fold-change'),
                legendPosition = 'left',
                legendLabSize = 10,
                legendIconSize = 4.0,
                drawConnectors = F,widthConnectors = 0.75
)+NoLegend()

