source("Loadings.R")

# > Oligodendrocytes CCA Integration ####################################################################
# Read in the oligodendrocyte Seurat objects
oligo_1stBatch <- readRDS("Sobj_1stSort_OligoSubset.rRNA_nCountFiltered_3Feb2022.rds")
oligo_2ndBatch <- readRDS("Sobj_Sorted2Last_OligoSubset_Aged_750Var_SCTres0.5_24Sep2021.rds")

oligo.list <- list(oligo_2ndBatch, oligo_1stBatch)

com.genes <- SelectIntegrationFeatures(object.list = oligo.list,
                                       assay = c("SCT","SCT"),
                                       nfeatures = 750
)
oligo.list <- PrepSCTIntegration(object.list = oligo.list, 
                                 assay = "SCT",
                                 anchor.features = com.genes, 
                                 verbose = TRUE)

oligo.list <- lapply(X = oligo.list, FUN = function(x) {
  x <- RunPCA(x, features = com.genes, verbose = T)
})

oligo.anchors <- FindIntegrationAnchors(object.list = oligo.list,
                                        #reduction = "rpca", 
                                        dims = 1:20,
                                        #k.anchor = 3,
                                        normalization.method = "SCT", 
                                        anchor.features = com.genes, verbose = TRUE)
oligo.integrated <- IntegrateData(anchorset = oligo.anchors, 
                                  dims = 1:20,
                                  normalization.method = "SCT", 
                                  verbose = TRUE)

DefaultAssay(oligo.integrated) <- "RNA"
oligo.integrated %<>% NormalizeData(.)
DefaultAssay(oligo.integrated) <- "integrated"

oligo.integrated$Batch <- ifelse(oligo.integrated$orig.ident == "Sorted_14" | oligo.integrated$orig.ident == "Sorted_15",
                                 "1stRound",
                                 "2ndRound"
)
# PCA ====================================================================
oligo.integrated %<>% RunPCA(., dims = 1:30, assay = 'integrated')

# PC selection for downstream
ElbowPlot(oligo.integrated, ndims = 20)
DimHeatmap(oligo.integrated, dims = c(1, 8:12))
# UMAP ====================================================================
oligo.integrated %<>% RunUMAP(., dims = 1:10,
                              seed.use = 6,
                              n.neighbors = 50,
                              #local.connectivity = 100,
                              #n.epochs = 500,
                              #min.dist = 0.05,
                              #reduction.name = "umap_dim15",
                              assay = "integrated"
)

# Marker gene expression example
myPalette <- colorRampPalette(c(
  rep("midnightblue", 3), rep("blue2", 1),
  rep("grey", 1),
  rep("orange", 1), rep("red", 3)
))(256)
DefaultAssay(oligo.integrated) <- "RNA"
FeaturePlot(oligo.integrated,
            features = 'Ifi27l2a',
            slot = 'data',
            order = T,
            pt.size = 2.5
) + scale_color_gradientn(colours = myPalette) + coord_fixed(ratio = 0.75)

# Clustering ====================================================================
DefaultAssay(oligo.integrated) <- 'integrated'
oligo.integrated %<>% FindNeighbors(., 
                                    reduction = 'pca',
                                    assay = 'integrated',
                                    dims = 1:10
                                    #prune.SNN = 1/10,
                                    #k.param = 20
)

oligo.integrated %<>% Seurat::FindClusters(., 
                                           algorithm = 1,
                                           resolution = c(0.4, 0.5, 0.6)
)
library(clustree)
clustree(oligo.integrated, prefix = 'integrated_snn_res.')

DimPlot(oligo.integrated, reduction = 'umap', pt.size = 2,
        group.by = 'integrated_snn_res.0.5',
        #split.by = 'orig.ident',
        label = T, label.size = 10
) + coord_fixed(ratio = 0.75)

oligo.clusters <- list(
  "0" = "Oligo2b",
  "1" = "Oligo1a",
  "2" = "Oligo1b",
  "3" = "Oligo1c",
  "4" = "Oligo2a",
  "5" = "ARO",
  "6" = "IRO"
)

oligo.integrated$clusters <- unlist(oligo.clusters[oligo.integrated$integrated_snn_res.0.5]) %>% as.factor()
oligo.integrated$clusters <- factor(oligo.integrated$clusters, 
                                    levels = levels(oligo.integrated$clusters)[c(3,4,5,6,7,1,2)])

oligo.clusters.merged <- list(
  'Oligo1a' = "Oligo1",
  'Oligo1b' = "Oligo1",
  'Oligo1c' = "Oligo1",
  'Oligo2a' = "Oligo2",
  'Oligo2b' = "Oligo2",
  'ARO' = "ARO",
  'IRO' = "IRO"
)

oligo.integrated$clusters_merged <- unlist(oligo.clusters.merged[oligo.integrated$clusters]) %>% as.factor()
oligo.integrated$clusters_merged <- factor(oligo.integrated$clusters_merged, 
                                           levels = levels(oligo.integrated$clusters_merged)[c(3,4,1,2)])

DimPlot(oligo.integrated,
        reduction = 'umap',
        group.by = 'clusters', cols = palette.10x.clusters,
        pt.size = 2.5) +
  coord_fixed(ratio = 0.75)
