#download (install.packages('')) and load (library ('') following packages
library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

#load matrix files created through aggr in cell ranger without depth normalization with file output called filtered_feature_bc_matrix
matrix_dir = "/Users/rajabn/Desktop/SSAnalysis/Merged_Aggr_Nonorm_FilteredFeature/"
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

#create a seurat object.
seurat_object <- CreateSeuratObject(counts = mat, project = "singlecell", min.cells = 3, min.features = 200)
seurat_object #to see what this consists on
dense.size <- object.size(as.matrix(mat))
dense.size
sparse.size <- object.size(mat)
sparse.size
dense.size/sparse.size

#Look at some metrics for setting cut-offs.
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-") #lower case mt if for mouse mitochondrial genes. capital for human.
head(seurat_object@meta.data, 20) #shows QC metrics for the first 20 cells
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#remove unwanted cells. 
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

#in order to determine the dimensionality of the dataset, we can first visualise the distribution of pvalues for each PC with a uniform distribution. 
#This can be visualised with eiher a JackStraw plot which will show ignificant PCs will show a strong enrichment of features with low pvalues.
#If run takes too long, an elbow plot can also be done.
# For running a JackStraw plot: seurat_object <- JackStraw(seurat_object, num.replicate = 100)
#seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20)
#JackStrawPlot(seurat_object, dims = 1:15)
#For an ElbowPlot:
ElbowPlot(seurat_object)

#to cluster cells, I have used the number of sig. PCs that I observed in the above plots/or in the case of an ElbowPlot, where I see levelling out of components. 
#The findneighbors function is for constuction of a KNN graph based on euclidean distance in PCA space and refines edge weights between any two cells based on the shared overlap in their local neighborhoods (jaccard similarity). It uses the input of previously defined dimensionality of the dataset.
seurat_object <- FindNeighbors(seurat_object, dims = 1:8)

#now to actually cluster the cells, we apply modularity optimisation techniques (default is Louvain algorithm).
# The findclusters function contains a resolution parameter which sets the granularity of downstream clustering. 
#I played with the resolution number in order to capture features/subsets in populations for constrast analysis.
#Subsets coloured consistently at a range of resolution settings for PSCMs. I kept this at 0.9. 
# At higher resolutions, some little featues in MDM clusters were coloured. The biology was explored, but noting high resolution was required to capture those aspects.
seurat_object <- FindClusters(seurat_object, resolution =0.90)
head(Idents(seurat_object), 5) #to have a look at cluster IDs of first 5 cells.

#run non-linear dimensional reduction (UMAP/tSNE)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "umap")

#run non-linear dimensional reduction (UMAP/tSNE)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
DimPlot(seurat_object, reduction = "umap")

#To colour with different colours:
p<-DimPlot(seurat_object, reduction = "umap")
p<-p + scale_color_manual(values = c('#E6194B', '#3CB44B', '#FFE119', '#4363D8', '#F58231', '#911EB4', '#46F0F0', '#F032E6', '#BCF60C', '#FABEBE', '#008080', '#E6BEFF', '#9A6324', '#FFFAC8', '#800000', '#AAFFC3', '#808000', '#FFD8B1', '#000075', '#808080', '#FFFFFF', '#000000'))
p


#find all markers distinguishing cluster 11 from cluster 8
seurat_object.PSCMcontrolsubsets <- FindMarkers(seurat_object, ident.1= 11, ident.2 = 8, min.pct=0.25)
head(seurat_object.PSCMcontrolsubsets, n=5)
write.table(seurat_object.PSCMcontrolsubsets, file = "C:/Users/rajabn/Desktop/SSAnalysis/PSCM_controlsubsets.tsv", quote = FALSE, sep = "\t" )

#find all markers distinguishing cluster 0 from 10
seurat_object.PSCMcontrolsubsets <- FindMarkers(seurat_object, ident.1= 0, ident.2 = 10, min.pct=0.25)
head(seurat_object.PSCMcontrolsubsets, n=5)
write.table(seurat_object.PSCMcontrolsubsets, file = "C:/Users/rajabn/Desktop/SSAnalysis/PSCM_acutesubsets.tsv", quote = FALSE, sep = "\t" )


-------------------------------------------

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
seurat_object@meta.data$sample_labels = sample_labels
DimPlot(seurat_object, reduction = "umap", group.by = "sample_labels")
sample_labels

#aggregation based on cluster name
df <- data.frame(row.names = 0:33537)
for (i_cluster in unique(seurat_object@meta.data$sample_labels))
{
  message(sprintf("Current cluster:%s\n", i_cluster))
  seurat_object@meta.data[seurat_object@meta.data[,"sample_labels"] == i_cluster,]
  temp <- mat[,row.names(seurat_object@meta.data[seurat_object@meta.data[,"sample_labels"] == i_cluster,])]
  print(dim(temp))
  colSums(temp!=0)
  colnames(temp) <- sample.int(floor(ncol(temp)/100), ncol(temp), replace = TRUE)
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
write.table(df, file = "C:/Users/rajabn/Desktop/SSAnalysis/merged_samplenames_100_new.tsv", quote = FALSE, sep = "\t" )
#make sample description file
df_samples <- data.frame(row.names = colnames(df))
df_samples$names <- row.names(df_samples)
df_samples <- separate(df_samples, col = names, into = c("cluster", "sub_cluster"), sep = "_")
write.table(df_samples, file = "C:/Users/rajabn/Desktop/SSAnalysis/Merged_seurat_samplenames_100_new.tsv", quote = FALSE, sep = "\t" )


--------------------------------
#pseudobulk DE analysis
BiocManager::install("DESeq2")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(df,
                              colData = df_samples,
                              design = ~ cluster)

dds <- DESeq(dds)
plotDispEsts(dds)
df_samples$cluster
contrast <- c("cluster", "MDM acute", "PSCM acute")
res <- results(dds,
               contrast = contrast,
               alpha = 0.05)
library("tibble")
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_tbl
write.csv(res_tbl,
          "C:/Users/rajabn/Desktop/SSAnalysis/de_MDM acute_vs_PSCM acute_genes.csv" ,
          quote = FALSE,
          row.names = FALSE)


-----------------------
                         
#Functions for different plots
VlnPlot(seurat_object, features = c("HMGB1"))
VlnPlot(seurat_object, features = c("IL1B"), slot = "counts", log = TRUE)
data("seurat_object")
cd_genes <- c("IL1B", "S100A9", "CCL4", "CCL3","CCL5","CCL4L2")
DotPlot(object = seurat_object, features = cd_genes)
FeaturePlot(seurat_object, features = c("IL1B", "S100A9", "CCL3", "CCL5", "CCL4", "CCL4L2"))
FeaturePlot(seurat_object, features = c("HMGB1"))

seurat_object@meta.data$group = sample_labels
VlnPlot(seurat_object, features = c("STAT2"), group.by = "group")
FeaturePlot(seurat_object, features = c("STAT2"))

----------------------
#cellxgene file output
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

SaveH5Seurat(seurat_object, filename = "/Users/nraja/Desktop/Single Cell processing via Seurat/cellxgene.h5Seurat")
Convert("/Users/nraja/Desktop/Single Cell processing via Seurat/cellxgene.h5Seurat", dest = "h5ad")


