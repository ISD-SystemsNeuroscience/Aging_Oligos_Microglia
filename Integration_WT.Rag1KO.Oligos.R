# > RPCA Integration Oligo RagKO + WT 24 mo ####################################################################
# Prepare the list of Seurat objects to integrate
oligo.RagKO_1 <- readRDS("Sobj_Rag1KO_Batch1_Oligo.rds")

oligo.RagKO_2 <- readRDS("Sobj_Rag1KO_Batch2_Oligo.rds")

##
oligo.WT_1 <- readRDS("Sobj_WT_Batch1_Oligo.rds")

oligo.WT_2 <- readRDS("Sobj_WT_Batch1_Oligo.rds")


# . SCT Integration --------------------------------------------------------------------
oligo.list <- list(oligo.WT_1, oligo.WT_2,oligo.RagKO_1, oligo.RagKO_2)

com.genes <- SelectIntegrationFeatures(object.list = oligo.list,
                                       assay = rep("SCT",4),
                                       nfeatures = 500
)

oligo.list <- PrepSCTIntegration(object.list = oligo.list,
                                 anchor.features = com.genes,
                                 verbose = TRUE)
oligo.list <- lapply(X = oligo.list, FUN = function(x) {
  x <- RunPCA(x, features = com.genes, 
              npcs = 20, 
              verbose = T)
})
oligo.anchors <- FindIntegrationAnchors(object.list = oligo.list,
                                        dims = 1:20,
                                        k.anchor = 3,
                                        reduction = "rpca",
                                        normalization.method = "SCT",
                                        anchor.features = com.genes, verbose = TRUE)
st <- matrix(c(-1,1,2,-2,-3,-4), ncol = 2)
oligo.merge <- IntegrateData(anchorset = oligo.anchors,
                             dims = 1:20,
                             normalization.method = "SCT",
                             sample.tree = st,
                             verbose = TRUE)

DefaultAssay(oligo.merge) <- "RNA"
oligo.merge %<>% NormalizeData(.)
DefaultAssay(oligo.merge) <- "integrated"

oligo.merge$Genotype_Tissue <- paste0(oligo.merge$Genotype,'_',oligo.merge$Region)
# QC ====================================================================
oligo.merge %<>% subset(., subset = nCount_RNA <= 25000)

VlnPlot(oligo.merge, 
        features = c("mito.percent", "nFeature_RNA", "nCount_RNA"),
        #cols = c(isdblue,isdorange),
        #cols = brewer.pal(8,"Paired")[c(1:2,5:6)],
        pt.size = 0.1,
        group.by = "orig.ident")
# PCA ====================================================================
oligo.merge %<>% RunPCA(., dims = 1:20, assay = 'integrated')

# Pick the number of PCs for downstream
ElbowPlot(oligo.merge, ndims = 20)
DimHeatmap(oligo.merge, dims = c(1:2, 12:15))

# UMAP ====================================================================
oligo.merge %<>% RunUMAP(., dims = 1:15,
                         seed.use = 6,
                         n.neighbors = 50,
                         #local.connectivity = 20,
                         assay = "integrated"
)

# Example feature plot
DefaultAssay(oligo.merge) <- "RNA"
myPalette <- colorRampPalette(c(
  rep("midnightblue", 3), rep("blue2", 1),
  rep("grey", 1),
  rep("orange", 1), rep("red", 2)
))(256)

FeaturePlot(oligo.merge,
            features = "Ifi27l2a",
            pt.size = 2,
            #split.by = "Genotype",
            #ncol = 2,
            order = T, 
            label = F
) & scale_color_gradientn(colours = myPalette) & coord_fixed(ratio = 0.8)

# Clustering ====================================================================
DefaultAssay(oligo.merge) <- "integrated"
oligo.merge %<>% FindNeighbors(., #features = VariableFeatures(oligo.merge)
                               reduction = 'pca',
                               dims = 1:15,
                               assay = 'integrated',
                               #prune.SNN = 1/10,
                               #k.param = 20
)

oligo.merge %<>% FindClusters(., 
                              resolution = c(0.4, 0.5, 0.6))

DimPlot(oligo.merge, reduction = 'umap', pt.size = 2,
        group.by = 'integrated_snn_res.0.5',
        label.size = 10,
) + coord_fixed(ratio = 0.8)

# Cluster Assignment ====================================================================
merge.list <- list(
  "0" = "Oligo1c",
  "1" = "Oligo1a",
  "2" = "Oligo2b",
  "3" = "ARO",
  "4" = "Oligo1b",
  "5" = "Oligo2a",
  "6" = "Oligo2c",
  "7" = "IRO"
)

oligo.merge$cluster.int <- unlist(merge.list[oligo.merge$integrated_snn_res.0.5]) %>% as.factor()
oligo.merge$cluster.int <- factor(oligo.merge$cluster.int,
                                  levels = levels(oligo.merge$cluster.int)[c(3:8,1,2)])

merge.list <- list(
  "0" = "Oligo1",
  "1" = "Oligo1",
  "2" = "Oligo2",
  "3" = "ARO",
  "4" = "Oligo1",
  "5" = "Oligo2",
  "6" = "Oligo2",
  "7" = "IRO"
)

oligo.merge$cluster.int.merged <- unlist(merge.list[oligo.merge$integrated_snn_res.0.5]) %>% as.factor()
oligo.merge$cluster.int.merged <- factor(oligo.merge$cluster.int.merged,
                                         levels = levels(oligo.merge$cluster.int.merged)[c(3,4,1,2)])

palette.10x.clusters.8 <- c("#a3e166", "#66cd00", "#51a400", "#f4de66", "#eec900", "#e0bb00", "#2d7dd2", "#e50000")
DimPlot(oligo.merge, reduction = 'umap', pt.size = 2,
        group.by = 'cluster.int', cols = palette.10x.clusters.8,
        label.size = 10,
) + coord_fixed(ratio = 0.8)