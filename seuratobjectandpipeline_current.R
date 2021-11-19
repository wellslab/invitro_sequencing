#download (install.packages('')) and load (library ('') following packages
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
#load matrix files created through aggr in cell ranger without depth normalization with file output called filtered_feature_bc_matrix
matrix_dir = "/Users/nraja/Desktop/Single Cell processing via Seurat/matrix/unnormalized cellranger output matrix files/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2 #ensure feature names V refers to the right V column so that gene IDs are picked up instead of ensembl IDs (important for removing mitochondrial genes later on)
#create a seurat object. I have altered min.features to be 90 as sample 7 as a median gene number of 94 per cell.
seurat_object <- CreateSeuratObject(counts = mat, project = "singlecell", min.cells = 3, min.features = 90)
seurat_object #to see what this consists on
dense.size <- object.size(as.matrix(mat))
dense.size
sparse.size <- object.size(mat)
sparse.size
dense.size/sparse.size
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") #lower case mt if for mouse mitochondrial genes. capital for human.
head(seurat_object@meta.data, 20) #shows QC metrics for the first 20 cells
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#remove unwanted cells. below are default settings but you can modify these
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) 
#now you have removed unwanted cells, it is time to normalize the data. By default, Seurat employs a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
#you can alternatively input seurat_object <- NormalizeData(seurat_object) instead of above.


VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#we now identify highly variable features in order to determine a subset of features that exhibit high cell-to-cell variation in the dataset.
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#now we apply a linear transformation (scaling) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)
#we can visualise both cells and features that define the PCA
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
print(seurat_object[["pca"]], dims = 1:8, nfeatures = 2)
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
DimPlot(seurat_object, reduction = "pca")
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)
#in order to determine the dimensionality of the dataset, we can first visualise the distribution of pvalues for each PC with a uniform distribution. Significant PCs will show a strong enrichment of features with low pvalues.
seurat_object <- JackStraw(seurat_object, num.replicate = 100)
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
JackStrawPlot(seurat_object, dims = 1:15)
#this can be visualised using an Elbow Plot
ElbowPlot(seurat_object)
#to cluster cells, I have used the number of sig. PCs that I observed in the above plots. The findneighbors function is for constuction of a KNN graph based on euclidean distance in PCA space and refines edge weights between any two cells based on the shared overlap in their local neighborhoods (jaccard similarity). It uses the input of previously defined dimensionality of the dataset.
seurat_object <- FindNeighbors(seurat_object, dims = 1:8)
#now to actually cluster the cells, we apply modularity optimisation techniques (default is Louvain algorithm). The findclusters function contains a resolution parameter which sets the granularity of downstream clustering. Settings are recommended between 0.4-1.2, but may need to be increased for larger datasets.
seurat_object <- FindClusters(seurat_object, resolution =0.80)
head(Idents(seurat_object), 5) #to have a look at cluster IDs of first 5 cells.

#run non-linear dimensional reduction (UMAP/tSNE)
seurat_object <- RunUMAP(seurat_object, dims = 1:8, resolution =0.80)
DimPlot(seurat_object, reduction = "umap")

# find markers for every cluster compared to all remaining cells, reporting only the positive ones (number of genes shown can be changed)
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

FeaturePlot(seurat_object, features = c("IL10"))
                         
seurat_object.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10
DoHeatmap(seurat_object, features = top10$gene) + NoLegend()      

VlnPlot(seurat_object, features = c("HMGB1"))
VlnPlot(seurat_object, features = c("IL1B"), slot = "counts", log = TRUE)
data("seurat_object")
cd_genes <- c("IL1B", "S100A9", "CCL4", "CCL3","CCL5","CCL4L2")
DotPlot(object = seurat_object, features = cd_genes)
FeaturePlot(seurat_object, features = c("IL1B", "S100A9", "CCL3", "CCL5", "CCL4", "CCL4L2"))
FeaturePlot(seurat_object, features = c("HMGB1"))

cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#to seperate out barcode.names file to identify samples
library(tidyr)
seurat_object@meta.data$names <- row.names(seurat_object@meta.data) #makes another column that has the row names of the data frame
groups <- separate(seurat_object@meta.data, col = names, into = c("a", "b"), sep = "-")
seurat_object <- RunUMAP(seurat_object, dims = 1:8)
seurat_object@meta.data$group = groups$b #added group column to seurat_object meta.data
DimPlot(seurat_object, reduction = "umap", group.by = "group")

# Rename groups
sample_names = c("PSCM control", "PSCM acute", "PSCM 18hr", "PSCM re-stimulation", "MDM control", "MDM acute", "MDM 18hr", "MDM re-stimulation")
sample_labels = as.factor(seurat_object@meta.data$group)
levels(sample_labels) = sample_names
sample_labels

seurat_object@meta.data$group = sample_labels
VlnPlot(seurat_object, features = c("STAT2"), group.by = "group")
FeaturePlot(seurat_object, features = c("STAT2"))

#cellxgene file output
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

SaveH5Seurat(seurat_object, filename = "/Users/nraja/Desktop/Single Cell processing via Seurat/cellxgene.h5Seurat")
Convert("/Users/nraja/Desktop/Single Cell processing via Seurat/cellxgene.h5Seurat", dest = "h5ad")

#setting up a loop for aggregation for projection
#we want to take for each cluster within the curly braces we want to select out the count matrix that corresponds to that count cluster (everything you want to do is in the braces)
df <- data.frame(row.names = 0:33537)
for (i_cluster in unique(seurat_object@meta.data$seurat_clusters)) 
{
  message(sprintf("Current cluster:%s\n", i_cluster))
  seurat_object@meta.data[seurat_object@meta.data[,"seurat_clusters"] == i_cluster,]
  temp <- mat[,row.names(seurat_object@meta.data[seurat_object@meta.data[,"seurat_clusters"] == i_cluster,])]
  print(dim(temp))
  colSums(temp!=0)
  colnames(temp) <- sample.int(floor(ncol(temp)/150), ncol(temp), replace = TRUE)
  print(colnames(temp))
  for (i_subset in unique(colnames(temp))) 
  {temp_subcluster <- temp[,colnames(temp)==i_subset]
  i_name<-paste(i_cluster, i_subset, sep= "_")
  print(i_name)
  df[i_name] <- rowSums(temp_subcluster)
  }  
}
rownames(df)<-feature.names$V1
colSums(df)
#export file now for projection on atlas
write.table(df, file = "C:/Users/nraja/Desktop/Single Cell processing via Seurat/df.tsv", quote = FALSE, sep = "\t" )
#make sample description file
df_samples <- data.frame(row.names = colnames(df))
df_samples$names <- row.names(df_samples)
df_samples <- separate(df_samples, col = names, into = c("cluster", "sub_cluster"), sep = "_")
write.table(df_samples, file = "C:/Users/nraja/Desktop/Single Cell processing via Seurat/df_seurat_samples.tsv", quote = FALSE, sep = "\t" )

seurat_object <- RunPCA(object = seurat_object)
DimPlot(object = seurat_object, reduction = "pca")
VariableFeaturePlot(object = seurat_object)
FeaturePlot(object = seurat_object, features = c("IL1B", "MMP9"), blend = TRUE)

levels(seurat_object)
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
seurat_object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(seurat_object.markers, file = "C:/Users/nraja/Desktop/Single Cell processing via Seurat/seurat_object.markers.tsv", quote = FALSE, sep = "\t" )

cluster.markers <- FindMarkers(seurat_object, ident.1 = 6, ident.2 = c(10), min.pct = 0.25)
head(cluster.markers, n = 10)

#aggregation based on cluster name
df <- data.frame(row.names = 0:33537)
for (i_cluster in unique(seurat_object@meta.data$group)) 
{
  message(sprintf("Current cluster:%s\n", i_cluster))
  seurat_object@meta.data[seurat_object@meta.data[,"group"] == i_cluster,]
  temp <- mat[,row.names(seurat_object@meta.data[seurat_object@meta.data[,"group"] == i_cluster,])]
  print(dim(temp))
  colSums(temp!=0)
  colnames(temp) <- sample.int(floor(ncol(temp)/50), ncol(temp), replace = TRUE)
  print(colnames(temp))
  for (i_subset in unique(colnames(temp))) 
  {temp_subcluster <- temp[,colnames(temp)==i_subset]
  i_name<-paste(i_cluster, i_subset, sep= "_")
  print(i_name)
  df[i_name] <- rowSums(temp_subcluster)
  }  
}
rownames(df)<-feature.names$V1
colSums(df)
#export file now for projection on atlas
write.table(df, file = "C:/Users/nraja/Desktop/Single Cell processing via Seurat/df_samplenames_50.tsv", quote = FALSE, sep = "\t" )
#make sample description file
df_samples <- data.frame(row.names = colnames(df))
df_samples$names <- row.names(df_samples)
df_samples <- separate(df_samples, col = names, into = c("cluster", "sub_cluster"), sep = "_")
write.table(df_samples, file = "C:/Users/nraja/Desktop/Single Cell processing via Seurat/df_seurat_samplenames_50.tsv", quote = FALSE, sep = "\t" )

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#pseudobulk DE analysis
BiocManager::install("DESeq2")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(df, 
                              colData = df_samples, 
                              design = ~ cluster)

rld <- rlog(dds, blind=TRUE)
DESeq2::plotPCA(rld, intgroup = "cluster")

rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
install.packages("pheatmap")
library(pheatmap)
pheatmap(rld_cor, annotation = df_samples[, c("cluster"), drop=F])

dds <- DESeq(dds)
plotDispEsts(dds)
df_samples$cluster
contrast <- c("cluster", "MDM acute", "PSCM acute")

res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)


res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_tbl
write.csv(res_tbl,
          "C:/Users/nraja/Desktop/Single Cell processing via Seurat/de_MDM acute_vs_PSCM acute_genes.csv" ,
          quote = FALSE, 
          row.names = FALSE)

VlnPlot(seurat_object, features = c("VIM"))
FeaturePlot(seurat_object, features = c("HLA-DRB1"))