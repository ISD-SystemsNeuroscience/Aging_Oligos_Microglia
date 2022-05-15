source("Loadings.R")

# This is 
# > Latest sorted libraries (2020-08-14) ####################################################################
df16 <- Read10X_h5("20200915_10x_2lanes/library_01/cellranger/outs/filtered_feature_bc_matrix.h5")
Sobj16 <- CreateSeuratObject(df16, project = "Library_16.Sorted", min.cells = 3, min.features = 200)
Sobj16$ExperimentDate <- "20200814"
Sobj16$Region <- "GreyM"
Sobj16$Gender <- "Male"
Sobj16$Age <- "24mo"
Sobj16$Genotype <- "WT"
#
df17 <- Read10X_h5("20200915_10x_2lanes/library_02/cellranger/outs/filtered_feature_bc_matrix.h5")
Sobj17 <- CreateSeuratObject(df17, project = "Library_17.Sorted", min.cells = 3, min.features = 200)
Sobj17$ExperimentDate <- "20200814"
Sobj17$Region <- "WhiteM"
Sobj17$Gender <- "Male"
Sobj17$Age <- "24mo"
Sobj17$Genotype <- "WT"

dataset.sorted <- merge(Sobj16, Sobj17)

# QC and filtering ====================================================================
dataset.sorted$mito.percent <- PercentageFeatureSet(dataset.sorted, pattern = "^mt-")

VlnPlot(dataset.sortedlast2, 
        features = c("mito.percent","nCount_RNA","nFeature_RNA"), pt.size = 0.1,
        cols = c(isdblue,isdorange),
        group.by = "Region")

dataset.sorted %<>% subset(., subset = mito.percent <= 10)
dataset.sorted %<>% subset(., subset = nCount_RNA <= 50000 & nFeature_RNA <= 8000)
# sctransform ====================================================================
dataset.sorted %<>% SCTransform(.,
                                assay = "RNA",
                                variable.features.n = NULL,
                                variable.features.rv.th = 1.4,
                                return.only.var.genes = F,
                                vars.to.regress = 'mito.percent',
                                do.scale = FALSE,
                                seed.use = 3
)

p0 <- VariableFeaturePlot(dataset.sorted)
top30 <- head(VariableFeatures(dataset.sorted), 30)

LabelPoints(plot = p0, points = top30, repel = T)
# PCA ====================================================================
dataset.sorted %<>% RunPCA(., dims = 1:30)
DimPlot(dataset.sorted, reduction = "pca", 
        group.by = 'Region', pt.size = 2)

ElbowPlot(dataset.sorted, ndims = 30)
DimHeatmap(dataset.sorted, dims = c(1:2, 29:30))

# UMAP ====================================================================
dataset.sorted %<>% RunUMAP(., dims = 1:30,
                            seed.use = 6,
                            n.neighbors = 50,
                            local.connectivity = 100,
                            n.epochs = 200,
                            #min.dist = 0.05,
                            #reduction.name = "umap_dim15",
                            assay = "SCT"
)

DimPlot(dataset.sorted,
        reduction = 'umap',
        #group.by = "Region", cols = c(isdblue,isdorange),
        group.by = 'CellTypeAnnotation2', cols = brewer.pal(8,"Paired"),
        label = T, label.size = 7, repel = T,
        pt.size = 1.5
) +
  coord_fixed(ratio = 1.25)

# Clustering ====================================================================
dataset.sorted %<>% FindNeighbors(., #features = VariableFeatures(dataset.sorted)
                                  reduction = 'pca', 
                                  dims = 1:30,
                                  #prune.SNN = 1/10,
                                  k.param = 20
)

dataset.sorted %<>% FindClusters(., algorithm = 2,
                                 resolution = c(0.1, 0.2, 0.3)
)

library(clustree)
clustree(dataset.sorted, prefix = 'SCT_snn_res.')

DimPlot(dataset.sorted, reduction = 'umap', pt.size = 1,
        group.by = 'SCT_snn_res.0.2',
        #split.by = 'orig.ident',
        label = T
)

# .. Find markers based on unbiased clustering --------------------------------------------------------------------
Idents(dataset.sorted) <- dataset.sorted$SCT_snn_res.0.1
markers_sct_res0.1 <- FindAllMarkers(dataset.sorted,test.use = "MAST", 
                                     only.pos = T, 
                                     min.diff.pct = 0.6,
                                     logfc.threshold = 1) %>% 
  arrange(cluster, desc(avg_logFC))
#
Idents(dataset.sorted) <- dataset.sorted$SCT_snn_res.0.2
markers_sct_res0.2 <- FindAllMarkers(dataset.sorted,test.use = "wilcox", only.pos = T, 
                                     #min.pct = 0.3,
                                     logfc.threshold = 0.5)
markers_sct_res0.2 %<>% filter(p_val_adj <= 0.05)

Idents(dataset.sorted) <- "SCT_snn_res.0.1"
markers_0vs2 <- FindMarkers(dataset.sorted, ident.1 = '0', ident.2 = '2',
                            test.use = "MAST",
                            only.pos = F
)
# Cell type annotation ====================================================================
markers.supp2 <- readRDS("supp2markergenes.rds")

FeaturePlot(dataset.sortedlast2, 
            #features = 'Bcas1', pt.size = 1, slot = 'data',
            features = "Sgk1", pt.size = 2,
            order = T
            #label = T, 
            #cols = c("#3E8EDE", "#F39E61")
) + scale_color_gradientn(colours = myPalette(256)) + coord_fixed(ratio = 1.25)

VlnPlot(dataset.sorted, features = "Bcas1", slot = 'scale.data')

FeaturePlot(dataset.sorted, 
            features = c('Ctss', 'Mbp', 'Aldoc', 'Sox4', 'Ccl5', 'S100a9'), 
            pt.size = 0.5,
            order = F,
            coord.fixed = T,
            ncol = 3
            #label = T, 
            #cols = c("#3E8EDE", "#F39E61")
) + scale_color_gradientn(colours = myPalette)
#coord_fixed(ratio = 1)

DotPlot(dataset.sorted,
        features = (c('Slc1a2','Ctss','S100a9','Sox4','Mbp','Ccl5','Ttr','Rgs5')),
        #features = c('Plp1', 'Ctss', 'Ccl5', 'Aldoc', 'Sox4', 'S100a9'),
        group.by = "CellTypeAnnotation2",
        #cols = c("#3E8EDE", "#F39E61"),
        dot.scale = 20
) + scale_color_gradientn(colours = myPalette)

cluster.markers <- list(
  '0' = 'Oligodendrocytes',
  '1' = 'Oligodendrocytes',
  '2' = 'Microglia',
  '3' = 'Oligodendrocytes',
  '4' = 'T/NK Cells',
  '5' = 'Astrocytes',
  '6' = 'NSCs/NBs',
  '7' = 'Mural/Ependymal Cells',
  '8' = 'Mural/Ependymal Cells',
  '9' = 'Neutrophils'
)

dataset.sorted$CellTypeAnnotation2 <- unlist(cluster.markers[dataset.sorted$SCT_snn_res.0.2])

dataset.sorted$CellTypeAnnotation2 %<>% as.factor(.)
dataset.sorted$CellTypeAnnotation2 <- factor(dataset.sorted$CellTypeAnnotation2,
                                             levels = levels(dataset.sorted$CellTypeAnnotation2)[c(1,2,4,5,6,7,3)]
)

# _ ----
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
