# install.packages('Seurat')
# install.packages('SeuratDisk')
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)

set.seed(123)

Convert("adataT.h5ad", dest = "h5seurat", overwrite = TRUE)
adata2seurat <- LoadH5Seurat("adataT.h5seurat")
adata2seurat

saveRDS(adata2seurat, file = "adataT2seurat.rds")

# readRDS(file, refhook = NULL)
adata2seurat2 <- readRDS("adataT2seurat.rds")

########################################################################
########################################################################

norm <- NormalizeData(adata2seurat2, normalization.method = "LogNormalize", scale.factor = 10000)

########################################################################
########################################################################

all.genes <- rownames(norm)
norm <- ScaleData(norm, features = all.genes)

########################################################################
########################################################################

norm <- RunPCA(norm, features = VariableFeatures(object = norm))

print(norm[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(norm, dims = 1:2, reduction = "pca")

DimPlot(norm, reduction = "pca")


norm <- JackStraw(norm, num.replicate = 100)
norm <- ScoreJackStraw(norm, dims = 1:50)
JackStrawPlot(norm, dims = 1:50)
ElbowPlot(norm)

########################################################################
########################################################################

norm <- FindNeighbors(norm, dims = 1:50)
norm <- FindClusters(norm, resolution = 0.5)

head(Idents(norm), 5)

norm <- RunUMAP(norm, dims = 1:50)
DimPlot(pbmc, reduction = "umap", label=TRUE)

saveRDS(pbmc, file = "reanalysis_output.rds")

########################################################################
########################################################################

norm.markers <- FindAllMarkers(norm, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

