######### Complete analysis of WT T. brucei Chromium scRNA-seq #########

## Script for analysis of T. brucei BSF parasites from Briggs et al. 2021
## Ran with R version 3.6.3 in Feburary 2021

### Required packages ###

# install.packages("BiocManager")
#install.packages("devtools")
#install.packages('Seurat')
#BiocManager::install("scran")
#BiocManager::install("scater")
#BiocManager::install("MAST")
#BiocManager::install("slingshot")
#BiocManager::install("tradeSeq")
#BiocManager::install("clusterExperiment")
#BiocManager::install("ggpubr")
#BiocManager::install("clustree")
#BiocManager::install("mgcv")
# For phate, first install phateR python package - https://github.com/KrishnaswamyLab/phateR
#install.packages("phateR")

###Load required packages##
# For prepocessing, replicate integration and clustering analysis
library(Seurat)
library(scran)
library(scater)
library(ggplot2)
library(RColorBrewer)
library(MAST)
library(SingleCellExperiment)
library(ggpubr)
library(clustree)

# For trajectory inference and DE analysis
library(slingshot)
library(tradeSeq)
library(clusterExperiment)
library(dplyr)
library(mgcv)
library(ComplexHeatmap)
library(phateR)
# Set path for the phate python module
Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")
reticulate::py_discover_config(required_module="phate")

## gglot2 plotting themes used for UMAP plots

UMAP_theme <- theme(axis.line=element_blank(), axis.ticks = element_blank(),  panel.background = element_rect(size=0.5,linetype="solid",color="black"), plot.title = element_text(size = 10, face = "bold", hjust = 0.05, vjust = -8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y =  element_blank(), legend.title = element_blank())


### Pre-processing of samples ###

# Load in VSGs and rRNA gene lists. These were generated by searching for "variant surface glycoprotein" and "rRNA" in the TREU927 genome on TriTrypDB. These will need to match the format of the scRNA object feature names.

VSG_list <- read.table("PATH-TO-FILE/VSG_list", quote="\"", comment.char="")
ribosomal_RNA <- ribosomal_RNA <- read.table("PATH-TO-FILE/ribosomal_RNA.txt", quote="\"", comment.char="")

# Load in data sets and generat seurat objects

WT_01_data <- Read10X(data.dir = "PATH-TO-FILE/filtered_feature_bc_matrix/")
WT_01_raw_object <- CreateSeuratObject(counts = cbind(WT_01_data), project = c("WT_01"), min.cells = 5, min.features = 3)

WT_02_data <- Read10X(data.dir = "PATH-TO-FILE/filtered_feature_bc_matrix/")
WT_02_raw_object <- CreateSeuratObject(counts = cbind(WT_02_data), project = c("WT_02"), min.cells = 5, min.features = 3)

# Add GEM classifications to object and subset. GEM classificatiosn were generated with Cellranger count function and can be found in Supplementary Data 1 (Briggs et al 2021)

#Read in as dataframe, were row names are cell barcodes
GEM_WT_01 <- as.data.frame(gem_classification)
GEM_WT_01 <- subset(GEM_WT_01, select = c("call"))
WT_01_raw_object <- AddMetaData(WT_01_raw_object, GEM_WT_01)
WT_01_raw_object_Tbrucei_subset <- subset(WT_01_raw_object, subset = call == "Tbrucei")

GEM_WT_02 <- as.data.frame(gem_classification_02)
GEM_WT_02 <- subset(GEM_WT_02, select = c("call"))
WT_02_raw_object <- AddMetaData(WT_02_raw_object, GEM_WT_02)
WT_02_raw_object_Tbrucei_subset <- subset(WT_02_raw_object, subset = call == "Tbrucei")

# Remove any transcripts mapping to L. mexicana in case of contamination 

WT_01_raw_object <- WT_01_raw_object_Tbrucei_subset[grep("Tbrucei", rownames(WT_01_raw_object_Tbrucei_subset)), ]
WT_02_raw_object <- WT_02_raw_object_Tbrucei_subset[grep("Tbrucei", rownames(WT_02_raw_object_Tbrucei_subset)), ]


# Calculate % of transcripts from KDNA

mt_marker_genes <- c("Tbrucei---ND7", "Tbrucei---COIII", "Tbrucei---MURF4", "Tbrucei---MURF1", "Tbrucei---ND1", "Tbrucei---MURF2", "Tbrucei---COI", "Tbrucei---COII", "Tbrucei---ND4", "Tbrucei---RPS12", "Tbrucei---ND5", "Tbrucei---Cyb")
WT_01_raw_object[["percent.MT"]] <- PercentageFeatureSet(WT_01_raw_object, features = mt_marker_genes)
WT_02_raw_object[["percent.MT"]] <- PercentageFeatureSet(WT_02_raw_object, features = mt_marker_genes)

# Calculate % of transcripts which are rRNA contaminants (need list of rRNA genes loaded, see above)
ribosomal_RNA <- as.character(ribosomal_RNA$V1)
WT_01_genes <- WT_01_raw_object@assays[["RNA"]]@counts
WT_01_genes <- c(WT_01_genes@Dimnames[[1]])
WT_01_rRNA_genes <- subset(ribosomal_RNA, ribosomal_RNA %in% WT_01_genes)
WT_01_raw_object[["percent.rRNA"]] <- PercentageFeatureSet(WT_01_raw_object, features = WT_01_rRNA_genes)

WT_02_genes <- WT_02_raw_object@assays[["RNA"]]@counts
WT_02_genes <- c(WT_02_genes@Dimnames[[1]])
WT_02_rRNA_genes <- subset(ribosomal_RNA, ribosomal_RNA %in% WT_02_genes)
WT_02_raw_object[["percent.rRNA"]] <- PercentageFeatureSet(WT_02_raw_object, features = WT_02_rRNA_genes)

# Visualize as scatter plots

UMI_MT <- FeatureScatter(WT_02_raw_object, feature1 = "nCount_RNA", feature2 = "percent.MT", pt.size = 0.1, smooth = FALSE) + labs(title = "") +
  scale_y_continuous(name = "% kDNA genes") + xlab("UMI Count") + 
  theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 9), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) +
  geom_hline(aes(yintercept= 2), color = "red", linetype="longdash") +
  geom_vline(aes(xintercept=4000), color = "red", linetype="longdash") 

UMI_rRNA <- FeatureScatter(WT_02_raw_object, feature1 = "nCount_RNA", feature2 = "percent.rRNA", pt.size = 0.1, smooth = FALSE) + labs(title = "") +
  scale_y_continuous(name = "% rRNA genes") + xlab("UMI Count") + 
  theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 9), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) +
  geom_hline(aes(yintercept=8), color = "red", linetype="longdash") +
  geom_vline(aes(xintercept=4000), color = "red", linetype="longdash")

UMI_feature <- FeatureScatter(WT_02_raw_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1, smooth = FALSE) + labs(title = "") +
  ylab("Feature Count") + xlab("UMI Count") +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) +
  geom_vline(aes(xintercept=4000), color = "red", linetype="longdash")  +
  geom_hline(aes(yintercept=250), color = "red", linetype="longdash") +
  geom_hline(aes(yintercept=2500), color = "red", linetype="longdash")

#pdf(file = "WT_QC_scatter_plots.pdf", width = 6.5, height = 2.5)
ggarrange(UMI_feature, UMI_MT, UMI_rRNA,  ncol = 3, legend = "none")
#dev.off()

# Select a subset of cells

WT_01 <- subset(WT_01_raw_object, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & nCount_RNA > 1000 & nCount_RNA < 4000 & percent.MT < 2 & percent.rRNA < 8)
WT_02 <- subset(WT_02_raw_object, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & nCount_RNA > 1000 & nCount_RNA < 4000 & percent.MT < 2 & percent.rRNA < 8)

# Remove rRNA from analysis
WT_01 <- WT_01[! rownames(WT_01) %in% WT_01_rRNA_genes,]
WT_02 <- WT_02[! rownames(WT_02) %in% WT_02_rRNA_genes,]

## Normalise and log transform data with SCRAN ##

# Export to a SingleCellExperiment object 

sce_WT_01 <- as.SingleCellExperiment(WT_01)
sce_WT_02 <- as.SingleCellExperiment(WT_02)

# Pre-cluster cells. Factors are first generated within clusters then rescaled to normalize between clusters
qclust_01 <- scran::quickCluster(sce_WT_01, min.size = 30)
qclust_02 <- scran::quickCluster(sce_WT_02, min.size = 30)

# Compute size factors - removes low abundance genes
sce_WT_01 <- scran::computeSumFactors(sce_WT_01, clusters = qclust_01)
sce_WT_01 <- scater::logNormCounts(sce_WT_01)

sce_WT_02 <- scran::computeSumFactors(sce_WT_02, clusters = qclust_02)
sce_WT_02 <- scater::logNormCounts(sce_WT_02)

## Convert back to Seurat object 
WT_01 <- as.Seurat(sce_WT_01, counts = "counts", data = "logcounts")
WT_02 <- as.Seurat(sce_WT_02, counts = "counts", data = "logcounts")

## Detect variable genes and remove VSG for variable gene list

# find vairble genes with scran
dec_01 <- modelGeneVar(sce_WT_01)
dec_02 <- modelGeneVar(sce_WT_02)

top.hvgs2_01 <- getTopHVGs(dec_01, n=3000)
top.hvgs2_02 <- getTopHVGs(dec_02, n=3000)

## variable genes with seurat
WT_01.features <- FindVariableFeatures(WT_01, selection.method = "vst", nfeatures = 3000, assay = "RNA") 
WT_02.features <- FindVariableFeatures(WT_02, selection.method = "vst", nfeatures = 3000, assay = "RNA") 

top3000_01 <- head(VariableFeatures(WT_01.features), 3000)
top3000_02 <- head(VariableFeatures(WT_02.features), 3000)

## Find those in common
common_var_genes_01 <- intersect(top.hvgs2_01, top3000_01)
common_var_genes_02 <- intersect(top.hvgs2_02, top3000_02)

## remove the VSGs
VSGs <-as.character(VSG_list$V1)

var_genes_01 <- subset(common_var_genes_01, ! common_var_genes_01 %in% VSGs)
var_genes_02 <- subset(common_var_genes_02, ! common_var_genes_02 %in% VSGs)

# Add variable genes to objects
WT_01@assays[["RNA"]]@var.features <- var_genes_01
WT_02@assays[["RNA"]]@var.features <- var_genes_02

# Save variable genes

write.csv(var_genes_01, file = "Common_variable_genes_WT_01")
write.csv(var_genes_02, file = "Common_variable_genes_WT_02")

# Save individual normalised objects
save(WT_01, file = "WT_01_seurat_norm")
save(WT_02, file = "WT_02_seurat_norm")

## Seurat integration ##

list <- list(WT_01, WT_02)

# Find common features for integratuon 
WT_features <- SelectIntegrationFeatures(object.list = list)

# Find integration anchors - 8 dims was found to give the most robust results by performing several iteration, varying this parameter and asses expression of known marker genes
WT.tryp.anchors <- FindIntegrationAnchors(object.list = list,  dims = 1:8, anchor.features = WT_features)

# Get all genes in both replicate experiment objects
all.genes <- row.names(list[[1]])
for (i in 2:length(list)) {
  all.genes <- intersect(all.genes, row.names(list[[i]]))
}

# Integrate samples
WT.integrated <- IntegrateData(anchorset = WT.tryp.anchors,  dims = 1:8, features.to.integrate=all.genes)

# Scale the data and regress variable due to total RNA
WT.integrated <- ScaleData(WT.integrated, vars.to.regress = "nCount_RNA", features = all.genes)

## PCA analysis and selectoion ##

WT.integrated <- RunPCA(object = WT.integrated, verbose = FALSE, features = WT_features)

# Determine percent of variation associated with each PC
pct <- WT.integrated@reductions[["pca"]]@stdev / sum(WT.integrated@reductions[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
pdf(file = "dim_eblbow_plot.pdf", width = 4, height = 3)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + xlab("Cumulative % variance") + ylab("% of variance") + guides(color=guide_legend(title="> 0.1% variance")) +
  geom_vline(xintercept = 90, color = "black", linetype="dashed") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + 
  theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), 
        axis.text.y =  element_text(size = 10), legend.text = element_text(size = 10))
dev.off()

## Clustering analysis
DefaultAssay(WT.integrated) <- "integrated"
WT.integrated <- FindNeighbors(WT.integrated, dims = 1:8, k.param = 30, nn.method = "annoy", annoy.metric = "euclidean")

## Run cluster to see results of different clustering resolutions
WT_clustree <- FindClusters(WT.integrated, resolution = c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1))

pdf(file = "clustree.pdf", width = 4, height = 8)
clustree(WT_clustree, prefix = "integrated_snn_res.", exprs = c("data", "counts", "scale.data"), assay = NULL)
dev.off()

# Rerun clustering with selected resolution parameter
WT.integrated <- FindClusters(WT.integrated, resolution = 0.4)

## UMAP vizualisation

# Run UMAP with same number of dims as for FindNeighbors
WT.integrated <- RunUMAP(WT.integrated, dims = 1:8, min.dist = 0.1)

# Order clusters as needed 
levels(WT.integrated) <- c("3", "2", "1", "0")

# Plot cells by replicate experiment or seurat cluateer (changee group.by parameter)
p <- DimPlot(WT.integrated, pt.size = 0.5, group.by = "orig.ident") 
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
#pdf(file = "UMAP_WT_intergation_replicate.pdf", width = 2, height = 2.2)
p + NoLegend() + UMAP_theme
#dev.off()

# Plot coloured by total RNA per cell
#p <- FeaturePlot(WT.integrated, features = "nCount_RNA")
#p + geom_point(aes(p$UMAP_1, p$UMAP_2), alpha = 0.5)

## Save integrated object
#save(WT.integrated, file = "WT_integrated_seurat_8dim")

### Identify cell types ###

## Plot known marker raw expression

DefaultAssay(WT.integrated) <- "RNA"

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.6.4280", min.cutoff = 0, max.cutoff = 4)
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_GAPHD.pdf", width = 1.8, height = 2)
p + labs(title = "GAPDH", color = "Expression") + NoLegend() + UMAP_theme
dev.off()

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.10.14140", min.cutoff = 0, max.cutoff = 4)
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_PYK1.pdf", width = 1.8, height = 2)
p + labs(title = "PYK1", color = "Expression") + NoLegend() + UMAP_theme
dev.off()

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.7.5940", min.cutoff = 0, max.cutoff = 4)
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_PAD2.pdf", width = 1.8, height = 2)
p + labs(title = "PAD2", color = "Expression") + NoLegend() + UMAP_theme
dev.off()

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.10.10260", min.cutoff = 0, max.cutoff = 4) + labs(title = "EP1", color = "Expression")
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_EP1.pdf", width = 1.8, height = 2)
p + NoLegend() + UMAP_theme
dev.off()

### Find all cluster marker genes
tryp.markers <- FindAllMarkers(WT.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", assay = "RNA")
write.csv(tryp.markers, file = "res0-4_WT_integrated_markers")

## Plot Heatmap 
# Identify top genes for each cluster
marker_genes_Wt <- WT_ZC3H20.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# Plot , gene names were added to markers manually in excel. 
pdf(file = "Heatmap_0_4res_WT_int_top10.pdf", width = 6, height = 6)
DoHeatmap(WT.integrated, features = as.character(marker_genes_Wt$gene), label = FALSE, draw.lines = TRUE, assay = "integrated") + scale_y_discrete(labels = rev(as.character(marker_genes_Wt$name)))
dev.off()

## Name clusters

# Name clusters, order to match that set by levels above
new.cluster.ids <- c("LS A", "LS B", "SS A", "SS B")
names(new.cluster.ids) <- levels(WT.integrated)
WT.integrated <- RenameIdents(WT.integrated, new.cluster.ids)
# Plot to check cluster names
DimPlot(WT.integrated)

## Caluclate proportions of each cluster per replicate
cell_proportions <- as.data.frame(prop.table(table(Idents(WT.integrated), WT.integrated$orig.ident), margin = 2))
write.csv(cell_proportions, file = "cell_proportions.csv")

## Marker gene violin plots

DefaultAssay(WT.integrated) <- "RNA"
pdf(file = "marker_vlnPlots.pdf", height = 8, width = 11)
VlnPlot(WT.integrated, features = c("Tbrucei---Tb927.1.4310", "Tbrucei---Tb927.11.880","Tbrucei---Tb927.9.7470", "Tbrucei---COII",
                                    "Tbrucei---Tb927.10.10590", "Tbrucei---Tb927.10.8940","Tbrucei---Tb927.3.2230","Tbrucei---Tb927.5.810",
                                    "Tbrucei---Tb927.10.2190", "Tbrucei---Tb927.3.3270","Tbrucei---Tb927.10.7410", "Tbrucei---Tb927.10.10260"), pt.size = 0, ncol = 4) 
dev.off()

#### PHATE, trajectory inference and differential expression analysis ###

# Get variable features from the integrated object
features <- WT.integrated@assays[["integrated"]]@var.features

## Run PHATE with integrated data

## extract data, setting cell as rows and features as coloums
data <- WT.integrated@assays[["integrated"]]@data
data_subset <- t(data[features, ])

# Run phate
phate_output <- phate(data_subset)
# Add phate result to seurat object
phate_output <- as.matrix(phate_output)
colnames(x = phate_output) <- paste0("PHATE_", 1:ncol(x = phate_output))
phate.reduction <- CreateDimReducObject(
  embeddings = phate_output,
  key = "PHATE_",
  assay = "integrated")
WT.integrated@reductions$phate <- phate.reduction

# Plotting the PHATE result (a bug in seurat means axis limits need to be set manually)
# Set PHATE plotting theme
PHATE_theme <- theme(axis.line=element_blank(), axis.ticks = element_blank(),
                     panel.background = element_rect(size=0.5,linetype="solid",color="black"),
                     plot.title = element_text(size = 10, face = "bold", hjust = 0.95, vjust = -8),
                     axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                     axis.text.y =  element_blank(), legend.title = element_blank()) + ylim(-0.01, 0.01) + xlim(0.03, -0.023)

# Plot by replicate sample
p <- DimPlot(WT.integrated, group.by = "orig.ident", reduction = "phate")
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "WT_intergation_sample.pdf", width = 3, height = 3)
p + PHATE_theme + NoLegend()
dev.off()

# Plot by cluster
p <- DimPlot(WT.integrated, group.by = "ident", reduction = "phate") 
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "WT_intergation_ident_PHATE.pdf", width = 3, height = 3.2)
p  + PHATE_theme + ylim(-0.01, 0.01) + xlim(0.033, -0.023) + NoLegend()
dev.off()

# Save again to keep PHATE
save(WT.integrated, file = "WT_integrated_seurat")

### Trajectory inference with slingshot

# convert to sce
sce <- as.SingleCellExperiment(WT.integrated, assay = "integrated")

# Run slingshot, setting the starting cluster
sce <- slingshot(sce, reducedDim = 'PHATE', clusterLabels = sce@colData@listData[["ident"]], start.clus = "LS A")

## Plots

# Colour by cluster
mycolours <- c("#f8766d", "#7cae00", "#01bfc4", "#c77cff")
pdf(file = "WT_intergation_ident_PHATE.pdf", width = 4.5, height = 4.5)
plot(reducedDims(sce)$PHATE, col = mycolours[sce$ident], pch = 16, cex = 0.5, bty='l', xlim=rev(c(-0.022, 0.032)))
lines(SlingshotDataSet(sce), col = "black", lwd = 2)
dev.off()   

# colour by pseudotime value
colors <- colorRampPalette(brewer.pal(11,'YlOrRd')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
pdf(file = "WT_intergation_pseudo_PHATE_TI.pdf", width = 4.8, height = 5.3)
plot(reducedDims(sce)$PHATE, col = plotcol, alpha = 0.5,pch = 16, cex = 0.5, xlim=rev(c(-0.023, 0.033)), axes = FALSE, ann = FALSE, xaxt='n', yaxt='n')
lines(SlingshotDataSet(sce), col = "black", lwd = 2)
dev.off() 

## Colour cells by marker gene expression level
DefaultAssay(WT.integrated) <- "RNA"

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.6.4280", min.cutoff = 0, max.cutoff = 4, reduction = "phate") 
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "PHATE_WT_GAPHD.pdf", width = 1.8, height = 2)
p + labs(title = "GAPDH", color = "Expression") + NoLegend() + PHATE_theme + ylim(-0.01, 0.01) + xlim(0.033, -0.023)
dev.off()

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.10.14140", min.cutoff = 0, max.cutoff = 4, reduction = "phate") 
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "PHATE_WT_PYK1.pdf", width = 1.8, height = 2)
p + labs(title = "PYK1", color = "Expression") + NoLegend() + PHATE_theme + ylim(-0.01, 0.01) + xlim(0.033, -0.023)
dev.off()

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.7.5940", min.cutoff = 0, max.cutoff = 4, reduction = "phate") 
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "PHATE_WT_PAD2.pdf", width = 1.8, height = 2)
p + labs(title = "PAD2", color = "Expression") + NoLegend() + PHATE_theme + ylim(-0.01, 0.01) + xlim(0.033, -0.023)
dev.off()

p <- FeaturePlot(object = WT.integrated, features = "Tbrucei---Tb927.10.10260", min.cutoff = 0, max.cutoff = 4, reduction = "phate") 
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "PHATE_WT_EP1.pdf", width = 1.8, height = 2)
p + labs(title = "EP1", color = "Expression") + NoLegend() + PHATE_theme + ylim(-0.01, 0.01) + xlim(0.033, -0.023)
dev.off() 

save(sce, file = "sce_slingshot_integrated_WT")

### Differential exression with tradeSeq

# Extract clusters
clusters <- sce@colData@listData[["ident"]]
# Get the lineage, selecting the same starting cluster if any
lin <- getLineages(SlingshotDataSet(sce), clusterLabels = clusters, start.clus = "LS A")
# Get the curves
crv <- getCurves(lin)

# 2d diffusion map with the slingshot trajectory and cluster, just as a check
plotGeneCount(curve = crv, clusters = clusters)

# You can set this control to try and prevent overfitting the data
control <- gam.control()
control$maxit <- 1000

# Get the counts data for the genes you want to use 
counts <- as.matrix(WT.integrated@assays[["RNA"]]@counts)

# Determine the appropriate number of knots (read on trade-seq page how to select best number of knots (points on trajectory))
# Test 3 to 15 knots, select best based on plots
#icMat <- evaluateK(counts = counts, sds = crv, k=3:15, nGenes = 200, verbose=FALSE)

# fit the GAM to the data (selected 6 knots). This takes a bit of time. 
sce_GAM_all_genes <- fitGAM(counts = counts, sds = crv, nknots = 6, control = control)
# save the single cell with GAM 
save(sce_GAM_all_genes, file = "sce_GAM_WT_integrated")

# Test gene expression association with trajectory
assoRes <- associationTest(sce_GAM_all_genes)
write.csv(assoRes, file = "assoRes_WT_integrated.csv")

## Select genes for clustering into gene modules
## select p value < 0.05 and FC > 2
diff.genes <- subset(assoRes, pvalue < 0.05)
diff.genes <- subset(diff.genes, meanLogFC > 0.301)
#Get list of genes
diff.genes <- as.character(diff.genes$Gene)
write.csv(diff.genes, file = "diff_genes.csv")

# Set number of points over trajector to test
nPointsClus <- 100

# Cluster genes and merge into coexpression modules
clusPat <- clusterExpressionPatterns(sce_GAM_all_genes, nPoints = nPointsClus,
                                     genes = diff.genes, mergeMethod="adjP", mergeCutoff=0.95)
save(clusPat, file = "clusPat_WT_inetgrated_100")

# Can test different merging thresholds but just running RSEC. 
RSEC <- RSEC(clusPat$rsec, eraseOld = FALSE, rerunClusterMany = FALSE, mergeMethod="adjP", mergeCutoff=0.95)

save(RSEC, file = "RSEC_WT_integrated")

## Plot expression of genes as heatmap using scaled expression of differentially expressed genes


# Get gene modules
clusterLabels <- primaryCluster(RSEC)
clusters_uniq <- unique(clusterLabels)

# Get scaled expression data for plotting
geneId <- rownames(clusPat$yhatScaled)
data <- clusPat$yhatScaled
rownames(data) <- geneId

# Order genes based on gene module (manually set)
data <- data[, seq_len(ncol(data))]
data <- cbind(data,clusterLabels)
data <- data[order(factor(data[,ncol(data)], levels = c("-1", "5", "6", "9", "1", "4", "3", "2", "8", "7"))),]

# Get gene clusters for each gene for labelling the heatmap, then remove from plotting data
clusters <- data[, 101]                               
write.csv(data, file = "Gene_modules_WT_integrated.csv")

data <-data[,-101]

# plot heatmap
colors <- setNames(colorRampPalette(brewer.pal(11,'YlOrRd')[-6])(100), 1:100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
ha <- HeatmapAnnotation(pseudotime = 1:100, col = list(pseudotime = colors),  show_legend = FALSE)
row_ha <- rowAnnotation(cluster = clusters, col = list(cluster = c("-1" = "grey",
                                                                   "5"="#F8766D",
                                                                   "6"="#DE8C00",
                                                                   "9"="#B79F00",
                                                                   "1"="#00BA38",
                                                                   "4"="#00BFC4", 
                                                                   "3"="#619CFF",
                                                                   "2"="#C77CFF",
                                                                   "8"="#F564E3",
                                                                   "7"="#FF64B0")))
                                                              

pdf(file = "WT_module_heatmap.pdf", width = 4.5, height = 5)
Heatmap(data, cluster_columns = FALSE, show_column_names = FALSE, cluster_rows = FALSE, show_row_names = FALSE,
        show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "expression"),
        top_annotation = ha, left_annotation = row_ha)
dev.off()

## Plot individual gene expression
#Add cluster info to the sce_GAM object for plotting

cluster_info <- WT.integrated@active.ident
sce_GAM_all_genes$cluster <- cluster_info

# Get raw counts data
counts <- WT.integrated@assays[["RNA"]]@counts

# Plot

pdf(file = "Tb927.11.15100_WT_smooth.pdf", width = 1.8, height = 1.5)
plotSmoothers(sce_GAM_all_genes, counts, gene = "Tbrucei---Tb927.11.15100", lwd = 1, sample = 1, alpha = 0.1, pointCol = "cluster") +
  labs(title = "Tb927.11.15100") + scale_color_manual(values = c("#f8766d", "#7cae00", "#01bfc4", "#c77cff")) + guides(color = FALSE) 
  theme(plot.title = element_text(size = 9, face = "bold"), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y =  element_text(size = 8))
  dev.off()
  
  
### DE analysis of slender vs stumpy cells for each replicate experiment, results detailed in Figure S3

#Rename seurat object to avoid confusion

WT.integrated_02 <- WT.integrated
# Repeat clustering with lower resolution to cluster slender and stumpy cell types
WT.integrated_02 <- FindClusters(WT.integrated_02, resolution = 0.1)
DimPlot(WT.integrated_02)    

# Subset cells from each replicate experiment to analysis seperately 
WT_01 <- subset(WT.integrated_02, orig.ident == "WT_01")
WT_02 <- subset(WT.integrated_02, orig.ident == "WT_02")

# Plot each replicate
p <- DimPlot(WT_02, pt.size = 0.5) 
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_02_clus.pdf", width = 2, height = 2.2)
p + NoLegend() + UMAP_theme
dev.off()

# Identify DE genes between slender and stumpy cells for each replicate

WT_01_markers <- FindAllMarkers(WT_01,test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", assay = "RNA")
write.csv(WT_01_markers, file = "WT_01_markers.csv")

WT_02_markers <- FindAllMarkers(WT_02, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST", assay = "RNA")
write.csv(WT_02_markers, file = "WT_02_markers.csv")

# Generate feature plots of individual genes
p <- FeaturePlot(object = WT_02, features = "Tbrucei---Tb927.9.14290", min.cutoff = 0, max.cutoff = 4) + labs(title = "CIF2", color = "Expression")
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_02_CIF2.pdf", width = 1.8, height = 2)
p + NoLegend() + UMAP_theme
dev.off()


genes <- tryp_subset@assays[["integrated"]]@var.features

slender_var_gene_pec <- PrctCellExpringGene(tryp_subset, genes = genes)

