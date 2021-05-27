## Cell cycle analysis 

## Circle trajectory inference
library(princurve)
library(tradeSeq)
library(cowplot)
library(clusterExperiment)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
# Plotting theme for graph asthetics
UMAP_theme <- theme(axis.line=element_blank(), axis.ticks = element_blank(),  panel.background = element_rect(size=0.5,linetype="solid",color="black"), plot.title = element_text(size = 10, face = "bold", hjust = 0.05, vjust = -8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y =  element_blank(), legend.title = element_blank())

####### PART ONE #######
## Assign cell cycle stage using marker genes
## Save the marker genes in the same format at my t. brucei genes, with the same headers.
## Call it "Cell_cycle_regulated_genes" and import it. 

## Cell cycle analysis with MetaFeature
tryp <- WT.integrated
# Get all the genes identified in your data set
genes <- tryp@assays[["RNA"]]@counts
genes <- genes@Dimnames[[1]]

# Get list of marker genes present in your data set
s.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$S.phase %in% genes)
s.genes <- s.genes$S.phase

g2m.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$G2.M.phase %in% genes)
g2m.genes <- g2m.genes$G2.M.phase

early.g1.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$Early.G1 %in% genes)
early.g1.genes <- early.g1.genes$Early.G1

late.g1.genes <- subset(Cell_cycle_regulated_genes, Cell_cycle_regulated_genes$Late.G1 %in% genes)
late.g1.genes <- late.g1.genes$Late.G1


# Calculate an expression score for each phase and save it to the seurat object
tryp <- MetaFeature(tryp, features = s.genes, meta.name = "S.aggregate")
tryp <- MetaFeature(tryp, features = g2m.genes, meta.name = "G2M.aggregate")
tryp <- MetaFeature(tryp, features = early.g1.genes, meta.name = "Early.G1.aggregate")
tryp <- MetaFeature(tryp, features = late.g1.genes, meta.name = "Late.G1.aggregate")

# Creat and dataframe with the expression score of each cell and each phase
df <- data.frame(tryp@meta.data[["S.aggregate"]], tryp@meta.data[["G2M.aggregate"]], tryp@meta.data[["Early.G1.aggregate"]], tryp@meta.data[["Late.G1.aggregate"]])
colnames(df) <- c("S", "G2M", "Early G1", "Late G1")
# Find the top scoring phase of each cell
cell_cycle_stage <- colnames(df)[apply(df,1,which.max)]
pident <- as.factor(cell_cycle_stage)
# Save the phase for each cell
tryp <- AddMetaData(tryp, pident, col.name = "Phase")

mycolours <- c("#f8766d", "#7cae00", "#01bfc4", "#c77cff", "grey")
# Plot the cells by phase on a UMAP
pdf(file = "PHATE_WT_phase.pdf", width = 3.7, height = 2.5)
DimPlot(object = tryp, group.by = "Phase", reduction = "umap") + UMAP_theme
dev.off()

cell_proportions <- as.data.frame(prop.table(table(tryp$Phase, tryp@active.ident), margin = 2))
write.csv(cell_proportions, file = "cell_proportions_phase_mutant_integration_testing_replicate.csv")
pdf(file = "WT_cell_proportions.pdf", width = 2.5, height = 2.5)
ggplot(data=cell_proportions, aes(x=cell_proportions$Var2, y=cell_proportions$Freq, fill=cell_proportions$Var1)) + geom_bar(stat="identity", color="black") + labs(x="sample", y="Proportion of Cells", fill="Cluster")
dev.off()



####### PART TWO #######

## Conduct peudotime analysis with the cycling cells, in my data the slender clusters

tryp_subset <- subset(tryp, subset = seurat_clusters == c("LS A", "LS B"))

# Find top variable genes for slenders only
tryp_subset <- FindVariableFeatures(tryp_subset, nfeatures = 2000, assay = "integrated")
slender_variable_genes <- tryp_subset@assays[["integrated"]]@var.features
write.csv(slender_variable_genes, file = "slender_variable_genes.csv")

# Repeat PCA analysis with these variable genes
tryp_subset  <- RunPCA(tryp_subset, verbose = FALSE)

## Recluster if you want to, but not nessary 
# Determine percent of variation associated with each PC
pct <- tryp_subset[["pca"]]@stdev / sum(tryp_subset[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
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
pdf(file = "dim_plot.pdf", width = 4, height = 3)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + xlab("Cumulative % variance") + ylab("% of variance") + guides(color=guide_legend(title="> 0.05% variance")) +
  geom_vline(xintercept = 90, color = "black", linetype="dashed") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + 
  theme(plot.title = element_text(size = 12, face = "bold"), axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11), axis.title.y = element_text(size = 11), 
        axis.text.y =  element_text(size = 10), legend.text = element_text(size = 10))
dev.off()
"integrated" -> DefaultAssay(tryp_subset)

tryp_subset <- FindNeighbors(tryp_subset, dims = 1:9, k.param = 30, nn.method = "annoy", annoy.metric = "euclidean")
# Can identify clusters
tryp_subset <- FindClusters(tryp_subset, resolution = 0.1)

tryp_subset <- RunUMAP(tryp_subset, dims = 1:9, reduction = "pca", min.dist = 0.1)

DimPlot(tryp_subset, reduction = "umap", label = FALSE,
        label.size = 4,
        pt.size = 0.05, group.by = "Phase") + UMAP_theme

pdf(file = "UMAP_slender_subset_cluster.pdf", width = 2.2, height = 2)
DimPlot(tryp_subset, reduction = "umap", group.by = "Phase") + 
  UMAP_theme + NoLegend()
dev.off()

save(tryp_subset, file = "slender_integrated_subset")

## Get UMAP dims
rd <- tryp_subset@reductions[["umap"]]@cell.embeddings
# check plot
mycolours <- c("#f8766d", "#7cae00", "#01bfc4", "#c77cff")
plot(rd, pch=16, col = mycolours[as.factor(tryp_subset$Phase)] ) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(),
        plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8))


# Draw a seperate curve for trajectory inference 
pcc <- principal_curve(rd, smoother="periodic_lowess")
#add line to plot
lines(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2], col="black", lwd=2)


colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pcc$lambda, breaks=100)]

data <- as.data.frame(rd)

data$pseudotime <- pcc$lambda

pdf(file = "UMAP_slender_pseudotime_princure.pdf", width = 3.1, height = 2)
ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) + geom_point(size = 1) + scale_color_gradientn(colours = colors) +
  geom_path(x=pcc$s[order(pcc$lambda),1], y=pcc$s[order(pcc$lambda),2], color = "black", lwd = 1) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(),
        plot.title = element_text(size = 10, face = "bold"), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9), axis.text.y =  element_text(size = 8))
dev.off()     

# fit smoothers on raw data
nPointsClus <- 50
counts <- as.matrix(tryp_subset@assays[["RNA"]]@counts)
cWeights <- rep(1,ncol(counts))
pseudoT <- matrix(pcc$lambda,nrow=ncol(counts),ncol=1)
gamList <- fitGAM(counts, pseudotime=pseudoT, cellWeights=cWeights, nknots=5)
save(gamList, file = "slenders_GAM")

# Test for association of expression with the trajectory
assocTestRes_cellCycle <- associationTest(gamList)
assoc.genes <- subset(assocTestRes_cellCycle, pvalue < 0.05)
assoc.genes <- subset(assoc.genes, meanLogFC > 0.301)
assoc.genes <- rownames(assoc.genes)

write.csv(assocTestRes_cellCycle, file = "assocTestRes_cellCycle.csv")

## Plot gene expression 

# Add phase info to object
gamList$phase <- tryp_subset$Phase

pdf(file = "PKA-R_cell_cycle_smooth.pdf", width = 1.8, height = 1.5)
plotSmoothers(gamList, counts, gene = "Tbrucei---Tb927.11.4610", lwd = 0.75, sample = 1, alpha = 0.5, pointCol = "phase") +
  labs(title = "PKA-R") + scale_color_manual(values = c("#f8766d", "#7cae00", "#01bfc4", "#c77cff")) + NoLegend() +
  theme(plot.title = element_text(size = 9, face = "bold"), axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8), axis.text.y =  element_text(size = 8))
dev.off()



clusPat_slender <- clusterExpressionPatterns(gamList, nPoints = nPointsClus,
                                     genes = assoc.genes)


# Get data
geneId <- rownames(clusPat_slender$yhatScaled)
data <- clusPat$yhatScaled
rownames(data) <- geneId

##Order genes based on order

data <- data[, seq_len(ncol(data))]
data <- cbind(data,clusterLabels)
data <- data[order(factor(data[,ncol(data)], levels = c("-1", "5", "6", "9", "1", "4", "3", "2", "8", "7"))),]

# Get gene clusters for each gene then remove from plotting data
clusters <- data[, 101]                               
write.csv(data, file = "Gene_modules_WT_integrated.csv")
data <-data[,-101]

# plot
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
        top_annotation = ha, left_annotation = row_ha, )
dev.off()



