
###R scripts used for Single-cell transcriptome dissection of the toxic impact of Di (2-ethylhexyl) phthalate on primordial follicle assembly.


#################=======DoubletFinder=====#############
###################################################################

library(Seurat)
library(cowplot)
library(patchwork)
library('DoubletFinder')

####PD0读入数据
DEHP <- Read10X(data.dir = "C:\\Users\\Desktop\\DEHP\\filtered_feature_bc_matrix")
colnames(x = DEHP) <- paste('DEHP', colnames(x = DEHP), sep = '_')

DEHP <- CreateSeuratObject(counts = DEHP, project = "DEHP", min.cells = 3, min.features = 200)
DEHP[["percent.mt"]] <- PercentageFeatureSet(DEHP, pattern = "^mt-")
VlnPlot(DEHP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)  
VlnPlot(DEHP, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)  

DEHP <- subset(DEHP, subset =  nFeature_RNA > 1500 & nFeature_RNA < 7500 )
DEHP <- NormalizeData(DEHP, normalization.method = "LogNormalize", scale.factor = 10000)
DEHP <- FindVariableFeatures(DEHP, selection.method = "vst", nfeatures = 2000) 
all.genes <- rownames(DEHP)
DEHP <- ScaleData(DEHP, features = all.genes)
DEHP@meta.data$tech <- "DEHP"
DEHP <- RunPCA(DEHP, npcs = 30, verbose = FALSE)
DEHP <- FindNeighbors(DEHP, reduction = "pca", dims = 1:15)
DEHP <- FindClusters(DEHP, resolution = 0.5)
DEHP <- RunTSNE(object = DEHP, dims.use = 1:15, do.fast = TRUE)
DEHP <- RunUMAP(DEHP, reduction = "pca", dims = 1:15)

##Remove Doublets
DEHPdouble <- paramSweep_v3(DEHP, PCs = 1:15, sct = FALSE)
DEHPdouble <- summarizeSweep(DEHPdouble, GT = FALSE)
DEHPdouble <- find.pK(DEHPdouble)
mpK1<-as.numeric(as.vector(DEHPdouble$pK[which.max(DEHPdouble$BCmetric)]))
annotations1 <- DEHP@meta.data$seurat_clusters
homotypic.prop1 <- modelHomotypic(annotations1)           
nExp_poi <- round(0.075*length(DEHP@active.ident))  

DEHP <- doubletFinder_v3(DEHP, PCs = 1:15, pN = 0.25, pK = 0.24, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

DimPlot(DEHP, reduction = "umap",group.by="DF.classifications_0.25_0.24_346",pt.size = 1)
DEHP_Filtered <- subset(DEHP, subset = DF.classifications_0.25_0.24_346 == "Singlet")



##################========Seurat version=3.1.5========############
### Detail information see online vignettes of Seurat
### https://satijalab.org/seurat/vignettes.html

library(dplyr)
library(Seurat)
library(patchwork)

###Load Data
Seurat_Object <- Read10X(data.dir = "D:\\...\\SampleName\\filtered_feature_bc_matrix\\")

Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = "Seurat_Object", min.cells = 3, min.features = 200)

Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Seurat_Object, pattern = "^mt-")

VlnPlot(Seurat_Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Seurat_Object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_Object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

###Value setting of filteration
Seurat_Object <- subset(Seurat_Object, subset = nFeature_RNA > value1 & nFeature_RNA < value2)

###======Perform integration
pancreas.anchors <- FindIntegrationAnchors(object.list = list(PD0, PD0-DEHP, PD3, PD3-DEHP), dims = 1:30)

DEHP.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(object = DEHP.integrated) <- "integrated"

DEHPintegrated <- ScaleData(DEHP.integrated, verbose = FALSE)
DEHP.integrated <- RunPCA(DEHP.integrated, npcs = 30, verbose = FALSE)
DEHP.integrated <- RunUMAP(DEHP.integrated, reduction = "pca", dims = 1:20)

#### t-SNE and Clustering

DEHP.integration <- FindNeighbors(DEHP.integrated, reduction = "pca", dims = 1:20)

DEHP.integration <- FindClusters(DEHP.integration, resolution = 0.3)
table(DEHP.integration@meta.data$seurat_clusters)

p1 <- DimPlot(DEHP.integration, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(DEHP.integration, reduction = "umap", group.by = "seurat_clusters")
plot_grid(p1, p2)

###find markers
DefaultAssay(DEHP.integration) <- "RNA"
DEHP.integration.All.Markers <- FindAllMarkers(DEHP.integration, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25





##subset single cell type

cell_type <- subset(PF.integration, idents = c("cell_type"))
saveRDS(cell_type, file = "D:\\...\\cell_type.rds")


########====== single cell type with seurat

Seurat_Germs <- readRDS(file = "D:\\...\\Germ_cells.rds")
Seurat_Germs

Seurat_Germs[["percent.mt"]] <- PercentageFeatureSet(Seurat_Germs, pattern = "^mt-")
VlnPlot(Seurat_Granulosa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Seurat_Germs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seurat_Germs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Germ_cells <- subset(Seurat_Germs, subset = nFeature_RNA > 200 & nFeature_RNA < ValueX )

######## Data processing

Germ_cells <- FindVariableFeatures(Germ_cells, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(Germ_cells), 10)

all.genes <- rownames(Germ_cells)
Germ_cells <- ScaleData(Germ_cells, features = all.genes)

###Perform linear dimensional reduction
Germ_cells <- RunPCA(Germ_cells, features = VariableFeatures(object = Germ_cells))

print(Germ_cells[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Germ_cells, dims = 1:2, reduction = "pca")

DimPlot(Germ_cells, reduction = "pca")
DimHeatmap(Germ_cells, dims = 1, cells = 500, balanced = TRUE)

###Determine the ‘dimensionality’ of the dataset
Germ_cells <- JackStraw(Germ_cells, num.replicate = 100)
Germ_cells <- ScoreJackStraw(Germ_cells, dims = 1:20)

JackStrawPlot(Germ_cells, dims = 1:15)
ElbowPlot(Germ_cells)

Germ_cells <- FindNeighbors(Germ_cells, dims = 1:10)
Germ_cells <- FindClusters(Germ_cells, resolution = 0.2)

#####Run non-linear dimensional reduction (UMAP/tSNE)

Germ_cells <- RunUMAP(Germ_cells, dims = 1:10)

##Marker
Germ_cluster.markers <- FindAllMarkers(Germ_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)




#########=======Monocle version=2.10.1========#############
### Detail information see online vignettes of Monocle
### http://cole-trapnell-lab.github.io/monocle-release/docs/

library("Seurat")
library('monocle')
Germs <- readRDS(file = "D:\\...\\Germs_seurat.rds")
Germs

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(Germs@assays$RNA@data), 'sparseMatrix')
dim(data)
head(data)[1:5,1:5]
pd <- new('AnnotatedDataFrame', data = Germs@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
Germ_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


dim(Germ_cds)
head(fData(Germ_cds))
head(pData(Germ_cds))

dim(exprs(Germ_cds))

###Estimate size factors and dispersions

Germ_cds <- estimateSizeFactors(Germ_cds)
Germ_cds <- estimateDispersions(Germ_cds)


##Trajectory step 1: choose genes that define a cell's progress
## Select genes that differ between clusters/stages
seurat_var_genes = VariableFeatures(Germs)
head(seurat_var_genes)

Germ_cds <- estimateSizeFactors(Germ_cds)
Germ_cds <- estimateDispersions(Germ_cds)


Germ_cds <- setOrderingFilter(Germ_cds, seurat_var_genes)
plot_ordering_genes(Germ_cds)

###Trajectory step 2: reduce data dimensionality
Germ_cds <- reduceDimension(Germ_cds,max_components = 2,method = 'DDRTree')

#Trajectory step 3: order cells along the trajectory
Germ_cds <- orderCells(Germ_cds)
plot_cell_trajectory(Germ_cds, color_by = "State")

##########BEAM Function
Germ_cds <- orderCells(Germ_cds,root_state = 2)
BEAM_res <- BEAM(Germ_cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

library("colorRamps")
library("RColorBrewer")
hmcols<-colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(65)
plot_genes_branched_heatmap(subGerm_cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                                         branch_point = 1,
                                         num_clusters = 4,
                                         cores = 1,
                                         hmcols = hmcols,
                                         use_gene_short_name = F,
                                         branch_colors = c("#bebebe","#009E73", "indianred2"),
                                         show_rownames = F,return_heatmap = T)




################=================Gene differential expression analysis

#####load the dataset 
library("Seurat")

DEHP_seurat <- readRDS(file = "C:\\Users\\admin\\Desktop\\DEHP_seurat.rds")
DEHP_seurat


DEHP_PD0.markers <- FindMarkers(DEHP_seurat, ident.1 = "PD0-DEHP", ident.2 = "PD0")
DEHP_PD3.markers <- FindMarkers(DEHP_seurat, ident.1 = "PD3-DEHP", ident.2 = "PD3")


PD0_DEHP_DEG <- filter(DEHP_PD0.markers, avg_logFC > 0.25 | avg_logFC < -0.25 , p_val_adj < 1.00e-02)
PD3_DEHP_DEG <- filter(DEHP_PD3.markers, avg_logFC > 0.25 | avg_logFC < -0.25 , p_val_adj < 1.00e-02)


library(VennDiagram)
PD0 <- rownames(PD0_DEHP_DEG)
PD3 <- rownames(PD3_DEHP_DEG)
DEHP_DEGs  <- intersect(PD0,PD3)
DEHP_DEGs <-data.frame(DEHP_DEGs)

input <- list(PD0,PD3)
names(input)<-c("PD0","PD3")
venn.diagram(input,"DEHP_DEGs.venn.tif",col = "transparent",cex = 1.5,
             fill = c("cornflowerblue", "violetred2"),alpha = 0.50,cat.cex = 2)


###火山图
DEHP_PD0.markers <- read.csv(file="DEHP_PD0.markers.csv")
DEHP_PD3.markers <- read.csv(file="DEHP_PD3.markers.csv")

data <- data.frame(gene = row.names(DEHP_PD0.markers),
                   pval = -log10(DEHP_PD0.markers$p_val_adj), 
                   lfc = DEHP_PD0.markers$avg_logFC)

# Remove any rows that have NA as an entry
data <- na.omit(data)

library(dplyr)
library(ggrepel)
# Color the points which are up or down


data <- mutate(data, color = case_when(data$lfc > 0.25 & data$pval > 2 ~ "Increased",
                                       data$lfc < 0.25 & data$pval > 2 ~ "Decreased",
                                       data$pval < 2 ~ "nonsignificant"))


library(ggplot2)
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# Add ggplot2 layers
vol + ggtitle(label = "Volcano Plot of DEHP") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) + 
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#F8766D", Decreased = "#00BFC4", nonsignificant = "darkgray")) +
  theme_bw(base_size = 16) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("PD0-DEHP" / "PD0"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 2, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p")  # Scale yaxis due to large p-values





#############=================Perform integration of Germ and Granulosa cell

library(Seurat)
library(cowplot)
library(patchwork)
Granulosa <- readRDS(file = "E:\\...\\Granulosa_seurat.rds")
Germs <- readRDS(file = "E:\\...\\Germs_seurat.rds")

pancreas.anchors <- FindIntegrationAnchors(object.list = list(Granulosa, Germs), dims = 1:30)
#Integration 
Germ_GCs.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
Germ_GCs.integrated

Germ_GCs.integrated <- ScaleData(Germ_GCs.integrated, verbose = FALSE)
Germ_GCs.integrated <- FindVariableFeatures(Germ_GCs.integrated, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(Germ_GCs.integrated), 10)
plot1 <- VariableFeaturePlot(Germ_GCs.integrated)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


Germ_GCs.integrated <- RunPCA(Germ_GCs.integrated, npcs = 30, verbose = FALSE)
Germ_GCs.integrated <- RunUMAP(Germ_GCs.integrated, reduction = "pca", dims = 1:20)


## t-SNE and Clustering

Germ_GCs.integrated <- FindNeighbors(Germ_GCs.integrated, reduction = "pca", dims = 1:20)

Germ_GCs.integrated <- FindClusters(Germ_GCs.integrated, resolution = 0.2)
table(Germ_GCs.integrated@meta.data$seurat_clusters)

p1 <- DimPlot(Germ_GCs.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(Germ_GCs.integrated, reduction = "umap", group.by = "seurat_clusters")
plot_grid(p1, p2)


p1 <- DimPlot(Germ_GCs.integrated, reduction = "umap",pt.size = 0.8, 
              group.by = "orig.ident", cols = c("mediumorchid2","#99c9fb","darkcyan"))
p2 <- DimPlot(Germ_GCs.integrated, reduction = "umap",pt.size = 0.8, 
              group.by = "seurat_clusters", label = TRUE)
plot_grid(p1, p2)

DimPlot(Germ_GCs.integrated, reduction = "umap", split.by = "orig.ident")

#### Rename

Germ_GCs.integrated_remane <- RenameIdents(Germ_GCs.integrated, 
                                     `1` = "Germ cells",
                                     `2` = "Germ cells", 
                                     `6` = "Germ cells",
                                     `0` = "Granulosa cells", 
                                     `3` = "Granulosa cells", 
                                     `4` = "Granulosa cells",
                                     `5` = "Granulosa cells")
DimPlot(Germ_GCs.integrated_remane, label = TRUE)

germ <- subset(Germ_GCs.integrated_remane, idents = c("Germ cells"))
gcs <- subset(Germ_GCs.integrated_remane, idents = c("Granulosa cells"))
DefaultAssay(object = germ) <- "RNA"
DefaultAssay(object = gcs) <- "RNA"
VlnPlot(object = germ, features = c("Bmp2","Bmp4","Bmp15","Gdf9", "Bmpr1a","Bmpr1b", "Bmpr2","Id1", "Id2","Id3","Id4",
                                         "Smad2","Smad3","Smad4","Smad5", "Smad7"), group.by="orig.ident",pt.size = 0.01)
VlnPlot(object = gcs, features = c("Bmp2","Bmp4","Bmp15","Gdf9", "Bmpr1a","Bmpr1b", "Bmpr2","Id1", "Id2","Id3","Id4",
                                             "Smad2","Smad3","Smad4","Smad5", "Smad7"), group.by="orig.ident",pt.size = 0.01)





