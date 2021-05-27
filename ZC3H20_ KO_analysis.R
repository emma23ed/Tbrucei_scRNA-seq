### Complete analysis of WT T. brucei Chromium scRNA-seq ###


#install.packages('Seurat')
#BiocManager::install("scran")
#BiocManager::install("scater")
#BiocManager::install("MAST")
#BiocManager::install("slingshot")
#BiocManager::install("tradeSeq")
#BiocManager::install("clusterExperiment")
#BiocManager::install("ggpubr")
#BiocManager::install("clustree")
#install.packages("phateR")
#BiocManager::install("ComplexHeatmap")
#remotes::install_github("carmonalab/STACAS")
#install.packages("bigmemory")


#Load packages
library(Seurat)
library(MAST)
library(scater)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(slingshot)
library(tradeSeq)
library(clusterExperiment)
library(ComplexHeatmap)
library(ggpubr)
library(dplyr)
library(STACAS)
library(gridExtra)

## Load in WT replicated already normalised and var genes found

# Load in VSGs and rRNA

VSG_list <- read.table("~/Documents/Work/Data/scRNA-seq/Complete_data/Final_analysis_NatureComms2021/VSG_list", quote="\"", comment.char="")
ribosomal_RNA <- read.table("~/Documents/Work/Data/scRNA-seq/Complete_data/Final_analysis_NatureComms2021/ribosomal_RNA.txt", quote="\"", comment.char="")

## MUTANT analysis

ZC3H20_data <- Read10X(data.dir = "/Users/emmabriggs/Documents/Work/Data/scRNA-seq/Complete_data/Final_analysis_NatureComms2021/ZC3H20_KO/filtered_feature_bc_matrix/")
ZC3H20_raw_object <- CreateSeuratObject(counts = ZC3H20_data, project = "ZC3H20_KO", min.cells = 5, min.features = 3)

gem_classification_ZC3H20 <- read.csv("~/Documents/Work/Data/scRNA-seq/Complete_data/Final_analysis_NatureComms2021/ZC3H20_KO/gem_classification_ZC3H20.csv", row.names=1)
GEM_ZC3H20 <- as.data.frame(gem_classification_ZC3H20)
GEM_ZC3H20 <- subset(GEM_ZC3H20, select = c("call"))
ZC3H20_raw_object <- AddMetaData(ZC3H20_raw_object, GEM_ZC3H20)
ZC3H20_raw_object_Tbrucei_subset <- subset(ZC3H20_raw_object, subset = call == "Tbrucei")

ZC3H20_raw_object <- ZC3H20_raw_object_Tbrucei_subset[grep("Tbrucei", rownames(ZC3H20_raw_object_Tbrucei_subset)), ]


# Calculate % of transcripts from KDNA
mt_marker_genes <- c("Tbrucei---ND7", "Tbrucei---COIII", "Tbrucei---MURF4", "Tbrucei---MURF1", "Tbrucei---ND1", "Tbrucei---MURF2", "Tbrucei---COI", "Tbrucei---COII", "Tbrucei---ND4", "Tbrucei---RPS12", "Tbrucei---ND5", "Tbrucei---Cyb")
ZC3H20_raw_object[["percent.MT"]] <- PercentageFeatureSet(ZC3H20_raw_object, features = mt_marker_genes)




# Calculate % of transcripts which are rRNA contaminants (need list of rRNA genes loaded)
ribosomal_RNA <- as.character(ribosomal_RNA$V1)
ZC3H20_genes <- ZC3H20_raw_object@assays[["RNA"]]@counts
ZC3H20_genes <- c(ZC3H20_genes@Dimnames[[1]])
ZC3H20_rRNA_genes <- subset(ribosomal_RNA, ribosomal_RNA %in% ZC3H20_genes)
ZC3H20_raw_object[["percent.rRNA"]] <- PercentageFeatureSet(ZC3H20_raw_object, features = "Tbrucei---tmp.1.10")

# Visualize QC metrics as a violin plot
featureVln <- VlnPlot(ZC3H20_raw_object, features = "nFeature_RNA", pt.size = 0) + labs(title = "Features") + scale_y_continuous(name = "Count", breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000)) + xlab("") + theme(plot.title = element_text(size = 10, face = "plain"), axis.text.x = element_blank(), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) + geom_hline(aes(yintercept=250), color = "red", linetype="longdash") + geom_hline(aes(yintercept=2250), color = "red", linetype="longdash")
UMIVln <- VlnPlot(ZC3H20_raw_object, features = "nCount_RNA", pt.size = 0, cols = "grey") + labs(title = "UMIs") + scale_y_continuous(name = "Count", breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000)) + xlab("Cells") + theme(plot.title = element_text(size = 10, face = "plain"), axis.text.x = element_blank(), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) + geom_hline(aes(yintercept=4000), color = "red", linetype="longdash")
MTVln <- VlnPlot(ZC3H20_raw_object, features = "percent.MT", pt.size = 0, cols = "grey") + labs(title = "kDNA Genes") + scale_y_continuous(name = "% Features") + xlab("Cells") + theme(plot.title = element_text(size = 10, face = "plain"), axis.text.x = element_blank(), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) + geom_hline(aes(yintercept=1.6), color = "red", linetype="longdash")
rRNAVln <- VlnPlot(ZC3H20_raw_object, features = "percent.rRNA", pt.size = 0, cols = "grey") + labs(title = "rRNA Genes") + scale_y_continuous(name = "% Features") + xlab("Cells") + theme(plot.title = element_text(size = 10, face = "plain"), axis.text.x = element_blank(), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) + geom_hline(aes(yintercept=8), color = "red", linetype="longdash")
pdf(file = "QC_vln_plots.pdf", width = 6.5, height = 2.5)
ggarrange(featureVln, UMIVln, MTVln, rRNAVln, ncol = 4, legend = "none", label.x = 0)
dev.off()

# Visualize as scatter plots
# ZC3H20_01

UMI_MT <- FeatureScatter(ZC3H20_raw_object, feature1 = "nCount_RNA", feature2 = "percent.MT", pt.size = 0.1, smooth = FALSE) + labs(title = "") +
  scale_y_continuous(name = "% kDNA genes") + xlab("UMI Count") + 
  theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 9), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) +
  geom_hline(aes(yintercept= 1), color = "red", linetype="longdash") +
  geom_vline(aes(xintercept=4000), color = "red", linetype="longdash") 

UMI_rRNA <- FeatureScatter(ZC3H20_raw_object, feature1 = "nCount_RNA", feature2 = "percent.rRNA", pt.size = 0.1, smooth = FALSE) + labs(title = "") +
  scale_y_continuous(name = "% rRNA genes") + xlab("UMI Count") + 
  theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 9), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) +
  geom_hline(aes(yintercept=3), color = "red", linetype="longdash") +
  geom_vline(aes(xintercept=4000), color = "red", linetype="longdash")

UMI_feature <- FeatureScatter(ZC3H20_raw_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1, smooth = FALSE) + labs(title = "") +
  ylab("Feature Count") + xlab("UMI Count") +
  theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8)) +
  geom_vline(aes(xintercept=4000), color = "red", linetype="longdash")  +
  geom_hline(aes(yintercept=250), color = "red", linetype="longdash") +
  geom_hline(aes(yintercept=2500), color = "red", linetype="longdash")

pdf(file = "ZC3H20_QC_scatter_plots.pdf", width = 6.5, height = 2.5)
ggarrange(UMI_feature, UMI_MT, UMI_rRNA,  ncol = 3, legend = "none")
dev.off()

# Select a subset of cells

ZC3H20 <- subset(ZC3H20_raw_object, subset = nFeature_RNA > 250 & nFeature_RNA < 2500 & nCount_RNA > 1000 & nCount_RNA < 4000 & percent.MT < 2 & percent.rRNA < 8)

# Remove rRNA from analysis
ZC3H20 <- ZC3H20[! rownames(ZC3H20) %in% ZC3H20_rRNA_genes,]

## Normalise data with SCRAN

# Export to a SingleCellExperiment object 

sce_ZC3H20 <- as.SingleCellExperiment(ZC3H20)

# Pre-cluster cells. Factors are first generated within clusters then rescaled to normalize between clusters
qclust <- scran::quickCluster(sce_ZC3H20, min.size = 30)

# Compute size factors - removes low abundance genes
sce_ZC3H20 <- scran::computeSumFactors(sce_ZC3H20, clusters = qclust)
sce_ZC3H20 <- scater::logNormCounts(sce_ZC3H20)

## Convert back to Seurat object 
ZC3H20 <- as.Seurat(sce_ZC3H20, counts = "counts", data = "logcounts")

## Detect variable genes and remove VSG
VSGs <-as.character(VSG_list$V1)

# find variable genes with scran

dec <- modelGeneVar(sce_ZC3H20)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)

top.hvgs2 <- getTopHVGs(dec, n=3000)
write.csv(top.hvgs2, file = "scran_HVGs")

## variable genes with seurat

ZC3H20.features <- FindVariableFeatures(ZC3H20, selection.method = "vst", nfeatures = 3000, assay = "RNA") 

top3000 <- head(VariableFeatures(ZC3H20.features), 3000)
write.csv(top3000, file = "seurat_HVGs")

## find those in common

common_var_genes <- intersect(top.hvgs2, top3000)

## remove the VSGs
var_genes <- subset(common_var_genes, ! common_var_genes %in% VSGs)

# Add to object 

ZC3H20@assays[["RNA"]]@var.features <- var_genes
write.csv(var_genes, file = "Common_variable_genes_ZC3H20")
save(ZC3H20, file = "ZC3H20_seurat_norm")

## STACAS integration
DefaultAssay(WT.integrated) <- "integrated"
list <- list(WT.integrated, ZC3H20)

tryp.features <- SelectIntegrationFeatures(object.list = list)

ref.anchors <- FindAnchors.STACAS(list, dims=1:8, anchor.features=tryp.features)


plots <- PlotAnchors.STACAS(ref.anchors)
g.cols <- 2
g.rows <- as.integer((length(plots)+2)/g.cols)
g <- do.call("arrangeGrob", c(plots, ncol=g.cols, nrow=g.rows))

plot(g)

ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors)

all.genes <- row.names(list[[1]])
for (i in 2:length(list)) {
  all.genes <- intersect(all.genes, row.names(list[[i]]))
}

mySampleTree <- SampleTree.STACAS(ref.anchors.filtered)

ref.integrated <- IntegrateData(anchorset=ref.anchors.filtered, dims=1:8, features.to.integrate=all.genes, sample.tree = mySampleTree)


## Scale data and regress RNA count
ref.integrated <- ScaleData(ref.integrated, verbose = TRUE, vars.to.regress = "nCount_RNA", features = all.genes)

# Run PCA using varible genes 
ref.integrated <- RunPCA(ref.integrated, verbose = FALSE)

# Determine percent of variation associated with each PC
pct <- ref.integrated[["pca"]]@stdev / sum(ref.integrated[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.05%.
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
  geom_text() + xlab("Cumulative % variance") + ylab("% of variance") + guides(color=guide_legend(title="> 0.05% variance")) +
  geom_vline(xintercept = 90, color = "black", linetype="dashed") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + 
  theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), 
        axis.text.y =  element_text(size = 10), legend.text = element_text(size = 10))
dev.off()

DefaultAssay(ref.integrated) <- "integrated"

ref.integrated <- FindNeighbors(ref.integrated, dims = 1:8, k.param = 30, nn.method = "annoy", annoy.metric = "euclidean", assay = "integrated")

ref.integrated <- FindClusters(ref.integrated, resolution = 0.4)
# Run UMAP with same number of dims as for FindNeighbors
ref.integrated <- RunUMAP(ref.integrated, reduction = "pca", dims = 1:7, min.dist = 0.01)

p <- DimPlot(ref.integrated, label = FALSE) 
p[[1]]$layers[[1]]$aes_params$alpha = 0.3
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_KO_intergation_cluster.pdf", width = 3, height = 3)
p + NoLegend() + UMAP_theme
dev.off()

# Plot coloured by total RNA per cell
p <- FeaturePlot(ref.integrated, features = "nCount_RNA")
p + geom_point(aes(p$UMAP_1, p$UMAP_2), alpha = 0.5)

save(ref.integrated, file = "mutant_WT_integrated_seurat_scaleall")


## Add line information
sample <- ref.integrated@meta.data[["orig.ident"]]
sample[sample == "WT_01"] <-  "WT"
sample[sample == "WT_02"] <-  "WT"
ref.integrated <- AddMetaData(ref.integrated, metadata = sample, col.name = "line")

# Generate plots

p <- DimPlot(ref.integrated, group.by = "line", label = FALSE, label.size = 5)
p[[1]]$layers[[1]]$aes_params$alpha = 0.3
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_ZC3H20.integrated_line.pdf", width = 3, height = 3)
p + UMAP_theme + NoLegend()
dev.off()

p <- DimPlot(ref.integrated, group.by = "line", label = FALSE, label.size = 5)
p[[1]]$layers[[1]]$aes_params$alpha = 0.3
p[[1]]$layers[[1]]$aes_params$shape = 16

pdf(file = "UMAP_WT_ZC3H20.integrated_line.pdf", width = 3, height = 3)
p + UMAP_theme + NoLegend()
dev.off()

## Plot known markers

p <- FeaturePlot(object = ref.integrated, features = "Tbrucei---Tb927.6.4280", min.cutoff = 0, max.cutoff = 4)
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_KO_GAPHD.pdf", width = 1.8, height = 2)
p + labs(title = "GAPDH", color = "Expression") + NoLegend() + UMAP_theme
dev.off()

p <- FeaturePlot(object = ref.integrated, features = "Tbrucei---Tb927.10.14140", min.cutoff = 0, max.cutoff = 4)
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_KO_PYK1.pdf", width = 1.8, height = 2)
p + labs(title = "PYK1", color = "Expression") + NoLegend() + UMAP_theme
dev.off()

p <- FeaturePlot(object = ref.integrated, features = "Tbrucei---Tb927.7.5940", min.cutoff = 0, max.cutoff = 4)
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_KO_PAD2.pdf", width = 1.8, height = 2)
p + labs(title = "PAD2", color = "Expression") + NoLegend() + UMAP_theme
dev.off()

p <- FeaturePlot(object = ref.integrated, features = "Tbrucei---Tb927.10.10260", min.cutoff = 0, max.cutoff = 4) + labs(title = "EP1", color = "Expression")
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_WT_KO_EP1.pdf", width = 1.8, height = 2)
p + NoLegend() + UMAP_theme
dev.off()

# Name clusters
new.cluster.ids <- c("LS A.1", "LS A.2", "LS B.1", "LS B.2","SS B", "SS A")
new.cluster.ids <- c("SS B", "LS B.1", "LS B.2", "SS A", "LS A.2", "LS A.1")
names(new.cluster.ids) <- levels(ref.integrated)
ref.integrated <- RenameIdents(ref.integrated, new.cluster.ids)
levels(ref.integrated) <- c("LS A.1", "LS A.2", "LS B.1", "LS B.2","SS A", "SS B")

p <- DimPlot(ref.integrated, group.by = "ident")
p[[1]]$layers[[1]]$aes_params$alpha = 0.5
p[[1]]$layers[[1]]$aes_params$shape = 16
pdf(file = "UMAP_intergation_tryp.integrated_ZC3H20_cluster.pdf", width = 4.6, height = 3.5)
p + UMAP_theme
dev.off()

## Find cluster proportions

cell_proportions <- as.data.frame(prop.table(table(Idents(ref.integrated), ref.integrated$line), margin = 2))
write.csv(cell_proportions, file = "cell_proportions_mutant_integration.csv")
pdf(file = "WT_cell_proportions.pdf", width = 2.5, height = 2.5)
ggplot(data=cell_proportions, aes(x=cell_proportions$Var2, y=cell_proportions$Freq, fill=cell_proportions$Var1)) + geom_bar(stat="identity", color="black") + labs(x="sample", y="Proportion of Cells", fill="Cluster")
dev.off()

## Markers between clusters
DefaultAssay(ref.integrated) <- "RNA"
WT_ZC3H20.markers <- FindAllMarkers(ref.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(WT_ZC3H20.markers, file = "WT_ZC3H20_integrated_markers.csv")

top10 <- WT_ZC3H20.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DefaultAssay(ref.integrated) <- "integrated"
pdf(file = "Heatmap_WT_ZC3H20_top_10_markers_RNA.pdf", width = 6, height = 6)
DoHeatmap(ref.integrated, features = top10$gene, label = FALSE)
dev.off()

## violin plots

DefaultAssay(ref.integrated) <- "RNA"
pdf(file = "LSB2_marker_vlnPlots.pdf", height = 3, width = 12)
VlnPlot(ref.integrated, features = c("Tbrucei---Tb927.7.6860",
                                     "Tbrucei---Tb927.11.11560",
                                     "Tbrucei---Tb2.NT.8",
                                     "Tbrucei---Tb927.10.1040",
                                     "Tbrucei---Tb927.4.710"), pt.size = 0, ncol = 5) 
dev.off()


## PHATE and trajectory analysis
features <- ref.integrated@assays[["integrated"]]@var.features

data <- ref.integrated@assays[["integrated"]]@data
data_subset <- t(data[features, ])
Sys.setenv(RETICULATE_PYTHON = "/usr/local/bin/python3")
phate_output <- phate(data_subset)


phate_output <- as.matrix(phate_output)
colnames(x = phate_output) <- paste0("PHATE_", 1:ncol(x = phate_output))

phate.reduction <- CreateDimReducObject(
  embeddings = phate_output,
  key = "PHATE_",
  assay = "integrated")

ref.integrated@reductions$phate <- phate.reduction

pdf(file = "PHATE_ref_intergation_line.pdf", width = 3, height = 3)
DimPlot(ref.integrated, group.by = "line", reduction = "phate") + NoLegend() + PHATE_theme + ylim(0.015, -0.015) + xlim(-0.025, 0.033)
dev.off()

pdf(file = "PHATE_ref_intergation_ident.pdf", width = 3, height = 3)
DimPlot(ref.integrated, group.by = "ident", reduction = "phate") + NoLegend() + PHATE_theme + ylim(0.015, -0.015) + xlim(-0.025, 0.033)
dev.off()

save(ref.integrated, file = "ref_integrated_seurat")

## convert to and infer trajectory

library(slingshot)

sce <- as.SingleCellExperiment(ref.integrated, assay = "integrated")

sce <- slingshot(sce, reducedDim = 'PHATE', clusterLabels = sce@colData@listData[["ident"]], start.clus = "LS A.1")
save(sce, file = "sce_int_mutant_WT_slingshot")

mycolours <- c("#d39200","#f8766d", "#7cae00", "#01bfc4", "#619cff", "#c77cff")

# Plot by ident
mycolours <- c("#f8766d", "#01bfc4")

sce$line <- as.factor(sce$line)
sce$phase <- ref.integrated$Phase
pdf(file = "ref_intergation_ident_PHATE.pdf", width = 4.5, height = 4.5)
plot(reducedDims(sce)$PHATE, col = mycolours[sce$line], pch = 16, cex = 0.5, bty='l', ylim = rev(c(-0.015, 0.015)),  axes = FALSE, ann = FALSE, xaxt='n', yaxt='n')
lines(SlingshotDataSet(sce), col = "black", lwd = 2)
dev.off()                                             

# Plot by pseudotime

colors <- colorRampPalette(brewer.pal(11,'YlOrRd')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
pdf(file = "ref_intergation_pseudo1_PHATE_Traj.pdf", width = 6, height = 6)
plot(reducedDims(sce)$PHATE, col = plotcol, pch = 16, cex = 0.5, bty='l', ylim = rev(c(-0.015, 0.015)), axes = FALSE, ann = FALSE, xaxt='n', yaxt='n')
lines(SlingshotDataSet(sce), col = "black", lwd = 2)
dev.off() 

colors <- colorRampPalette(brewer.pal(11,'YlGnBu')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_2, breaks=100)]
pdf(file = "ref_intergation_pseudo2_PHATE_Traj.pdf", width = 6, height = 6)
plot(reducedDims(sce)$PHATE, col = plotcol, pch = 16, cex = 0.5, bty='l', ylim = rev(c(-0.015, 0.015)), axes = FALSE, ann = FALSE, xaxt='n', yaxt='n')
lines(SlingshotDataSet(sce), col = "black", lwd = 2)
dev.off()

## Colour cells by gene
colors <- colorRampPalette(brewer.pal(11,'Blues')[-6])(100)
gene_logcounts <- sce@assays@data@listData[["logcounts"]]
plotcol <- colors[cut(gene_logcounts["Tb927.11.6870", ], breaks=100)]
pdf(file = "ref_intergation_PHATE_PYK1.pdf", width = 4.5, height = 4.5)
plot(reducedDims(sce)$PHATE, col = plotcol, pch = 16, cex = 0.5, bty='l', ylim = rev(c(-0.015, 0.015)))
dev.off() 

## Differential expression with tradeseq

library(tradeSeq)
library(mgcv)
library(clusterExperiment)
library(cowplot)

# You can set this control to try and prevent overfitting the data
control <- gam.control()
control$maxit <- 1000

# Get the counts data for the genes you want to use (here, top 2000)
# Get all counts
counts <- as.matrix(ref.integrated@assays[["RNA"]]@counts)

#all_genes <- rownames(ref.integrated)

#smoothers_all <- predictSmooth(sce_GAM_line, gene = all_genes, tidy = TRUE)

## Fit a generlised additive model using the linage curves identified by slingshot
# Get cluster ids
clusters <- sce@colData@listData[["ident"]]
# Get the lineage, selecting the same starting cluster if any
lin <- getLineages(SlingshotDataSet(sce), clusterLabels = clusters, start.clus = "LS A.1")
# Get the curves
crv <- getCurves(lin)

# Determine the appropriate number of knots (read on trade-seq page how to select best number of knots (points on trajectory))

# Test 3 to 15 knots
icMat <- evaluateK(counts = counts, sds = crv, k=3:15, nGenes = 200,
                   verbose=FALSE)

# fit the GAM to the data (I selected 9 knots). This takes a bit of time. 
sce_GAM_line <- fitGAM(counts = counts, sds = crv, nknots = 9, control = control, conditions =as.factor(sce$line))
  2# save the single cell with GAM 
save(sce_GAM_line, file = "sce_GAM_WT_ZC3H20_perline_all_genes2")

# Test gene expression association with trajectories, setting lineages to true will test each lineage seperatlely
assoRes_line <- associationTest(sce_GAM_line, lineage = TRUE)


head(assoRes)
write.csv(assoRes_line, file = "assoRes_line_WT_ZC3H20")

## Plot individual gene expression (full of bugs to had to make a work around)

smooth <- predictSmooth(sce_GAM_line, gene = "Tbrucei---Tb927.7.970")
smooth_data_WT <- subset(smooth, condition == "WT")
smooth_data_WT <- subset(smooth_data_WT, lineage == "1")
smooth_data_WT <- pivot_wider(smooth_data_WT, names_from = gene, values_from = yhat)

smooth_data_ZC3H20 <- subset(smooth, condition == "ZC3H20_KO")
smooth_data_ZC3H20 <- subset(smooth_data_ZC3H20, lineage == "2")
smooth_data_ZC3H20 <- pivot_wider(smooth_data_ZC3H20, names_from = gene, values_from = yhat)

pCol <- as.character("#f8766d")

mycolours <- c("#f8766d", "#00bfc4")
pdf(file = "NMD3_mutant_smooth.pdf", width = 1.5, height = 1.5)
plotSmoothers(sce_GAM_line, counts, gene = "Tbrucei---Tb927.7.970", lwd = 0.3, size = 1/10, pointCol = "line", plotLineages = FALSE) +
  labs(title = "NMD3") + scale_color_manual(values =c("#f8766d", "#f8766d", "#00bfc4")) + NoLegend() + geom_line(data = smooth_data_WT, aes(x = smooth_data_WT$time, y = smooth_data_WT$`Tbrucei---Tb927.7.970`)) + geom_line(data = smooth_data_ZC3H20, aes(x = smooth_data_ZC3H20$time, y = smooth_data_ZC3H20$`Tbrucei---Tb927.7.970`), color="#00bfc4") +
  theme(plot.title = element_text(size = 9, face = "bold"), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y =  element_text(size = 8)) 
dev.off()

## Plot relative expression of differentiation associated genes at heatmap

library(clusterExperiment)
library(ComplexHeatmap)

# Get differentiation assoaicted genes

diff.genes.mutantvWT <- row.names(assoRes_line)
differentiation.genes <- subset(diff.genes.mutantvWT, diff.genes.mutantvWT %in% diff_genes)
yhat_ZC3H20 <- predictSmooth(sce_GAM_ZC3H20, gene = diff.genes)


all_genes <- rownames(ref.integrated)
smoothers_all <- predictSmooth(sce_GAM_line, gene = all_genes, tidy = TRUE)
# Get line information and subset out only ZC3H20 KO cells
sce_GAM_line$line <- ref.integrated$line
sce_GAM_ZC3H20 <- subset(sce_GAM_line, , line == "ZC3H20_KO")
sce_GAM_ZC3H20 <- subset(sce_GAM_line_2, , lineage == "ZC3H20_KO")
# Set plots to analyse
nPointsClus <- 100

diff.genes <- as.character(Differentiation$V1)
# Run clustering and calculte relative expression over trajectory for KO cells

smooth <- predictSmooth(sce_GAM_line, gene = diff.genes, tidy = TRUE)


yhatPatScaled <- t(scale(t(yhatPat)))

clusPat_diff_genes_ZC3H20 <- clusterExpressionPatterns(sce_GAM_ZC3H20, nPoints = nPointsClus,
                                                       genes = order)

save(clusPat_diff_genes_ZC3H20, file = "clusPat_diff_genes_ZC3H20")

# Plot heatmap

smooth <- predictSmooth(sce_GAM_line, gene = order, tidy = TRUE)
smooth_subset_KO <- subset(smooth, smooth$condition == 'ZC3H20_KO')
smooth_subset_KO1 <- subset(smooth_subset_KO, smooth_subset$lineage == '1')
smooth_subset_KO2 <- subset(smooth_subset_KO, smooth_subset$lineage == '2')
smooth_subset_WT <- subset(smooth, smooth$condition == 'WT')
smooth_subset_WT1 <- subset(smooth_subset_WT, smooth_subset$lineage == '1')
smooth_subset_WT2 <- subset(smooth_subset_WT, smooth_subset$lineage == '2')



KO1 <- pivot_wider(smooth_subset_KO1, names_from = gene, values_from = yhat)
KO2 <- pivot_wider(smooth_subset_KO2, names_from = gene, values_from = yhat)
WT1 <- pivot_wider(smooth_subset_WT1, names_from = gene, values_from = yhat)
WT2 <- pivot_wider(smooth_subset_WT2, names_from = gene, values_from = yhat)


rownames(KO1) <- KO1$time
KO1 <- KO1[ ,!(names(KO1) %in% c("lineage", "time", "condition"))]
KO1 <- t(KO1)
rownames(KO2) <- KO2$time
KO2 <- KO2[ ,!(names(KO2) %in% c("lineage", "time", "condition"))]
KO2 <- t(KO2)
rownames(WT1) <- WT1$time
WT1 <- WT1[ ,!(names(WT1) %in% c("lineage", "time", "condition"))]
WT1 <- t(WT1)
rownames(WT2) <- WT2$time
WT2 <- WT2[ ,!(names(WT2) %in% c("lineage", "time", "condition"))]
WT2 <- t(WT2)

data2 <- cbind(KO1, KO2)


data_scale <- t(scale(t(data2)))
data_scale_subset <- data_scale[ ,101:200]


order <- rownames(data)
data_scale_order <- data_scale_subset[order,]


write.csv(data, file = "Gene_modules_WT_integrated.csv")

colors <- setNames(colorRampPalette(brewer.pal(11,'YlGnBu')[-6])(100), 1:100)
plotcol <- colors[cut(sce$slingPseudotime_2, breaks=100)]
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


pdf(file = "Zc3H20_module_heatmap4.pdf", width = 4.5, height = 5)
Heatmap(data_scale_order, cluster_columns = FALSE, show_column_names = FALSE, cluster_rows = FALSE, show_row_names = FALSE,
        show_heatmap_legend = TRUE, heatmap_legend_param = list(title = "expression"), 
        top_annotation = ha)
dev.off()


## Early differentially expressed gene analysis) ##

## plot the knot positions 
mycolours <- c("#d39200","#f8766d", "#7cae00", "#01bfc4","#c77cff", "#619cff")
sce@colData@listData[["line"]] -> line
crv <- AddMetaData(crv, metadata = line, col.name = "line")
pdf(file = "PHATE_cluster_knots.pdf", width = 7, height = 7)
plotGeneCount(curve = crv, counts = counts,
              clusters = apply(slingClusterLabels(crv), 1, which.max),
              models = sce_GAM_line) + scale_y_reverse() + NoLegend() + scale_color_manual(values= mycolours)
dev.off()

# find early differentialy expressed genes 
earlyDERes <- earlyDETest(sce_GAM_line, knots = c(1, 3))
write.csv(earlyDERes, file = "early_res_test.csv")


### Cell cycle ##

genes <- ref.integrated@assays[["integrated"]]@data
genes <- genes@Dimnames[[1]]

s.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$S.phase %in% genes)
s.genes <- s.genes$S.phase
g2m.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$G2.M.phase %in% genes)
g2m.genes <- g2m.genes$G2.M.phase
early.g1.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$Early.G1 %in% genes)
early.g1.genes <- early.g1.genes$Early.G1
late.g1.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$Late.G1 %in% genes)
late.g1.genes <- late.g1.genes$Late.G1

ref.integrated <- MetaFeature(ref.integrated, features = s.genes, meta.name = "S.aggregate")
ref.integrated <- MetaFeature(ref.integrated, features = g2m.genes, meta.name = "G2M.aggregate")
ref.integrated <- MetaFeature(ref.integrated, features = early.g1.genes, meta.name = "Early.G1.aggregate")
ref.integrated <- MetaFeature(ref.integrated, features = late.g1.genes, meta.name = "Late.G1.aggregate")

df <- data.frame(ref.integrated@meta.data[["S.aggregate"]], ref.integrated@meta.data[["G2M.aggregate"]], ref.integrated@meta.data[["Early.G1.aggregate"]], ref.integrated@meta.data[["Late.G1.aggregate"]])
colnames(df) <- c("S", "G2M", "Early G1", "Late G1")
cell_cycle_stage <- colnames(df)[apply(df,1,which.max)]
cell_cycle_stage <- new_phase_mutant$V1
pident <- as.factor(cell_cycle_stage)
ref.integrated <- AddMetaData(ref.integrated, pident, col.name = "Phase")

save(ref.integrated, file = "ref_zc3h20_integrated_object")

pdf(file = "UMAP_WT_ZC3H20_Phase.pdf", width = 5, height = 3.5)
DimPlot(ref.integrated, group.by = "Phase", label = FALSE,  cols = c("Early G1" = "#f8766d", "Late G1" = "#7cae00", "S" = "#01bfc4", "G2M" = "#c77cff", "Unassigned" = "grey")) + UMAP_theme
dev.off()

DimPlot(ref.integrated, group.by = "Phase") + UMAP_theme


cell_proportions <- as.data.frame(prop.table(table(ref.integrated$Phase, ref.integrated@active.ident), margin = 2))
write.csv(cell_proportions, file = "cell_proportions_phase_mutant_integration_new.csv")
pdf(file = "WT_cell_proportions.pdf", width = 2.5, height = 2.5)
ggplot(data=cell_proportions, aes(x=cell_proportions$Var2, y=cell_proportions$Freq, fill=cell_proportions$Var1)) + geom_bar(stat="identity", color="black") + labs(x="sample", y="Proportion of Cells", fill="Cluster")
dev.off()





