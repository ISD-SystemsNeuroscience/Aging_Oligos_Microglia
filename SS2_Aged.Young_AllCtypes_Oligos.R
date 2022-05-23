source("Loadings.R")

SS2_allctypes <- readRDS("Sobj_SS2_allctypes.rds")
SS2_allctypes$Age_Tissue <- paste0(SS2_allctypes$Age,'_',SS2_allctypes$Condition)


# QC and filtering ====================================================================
SS2_allctypes$mito.percent <- PercentageFeatureSet(SS2_allctypes, pattern = "^mt-")
SS2_allctypes$ReadCountsLog10 <- log10(SS2_allctypes$nCount_RNA)

VlnPlot(SS2_allctypes,
        pt.size = 0.1,
        features = c("nFeature_RNA"),
        cols = brewer.pal(8,"Paired"),
        group.by = "Animal")

SS2_allctypes %<>% subset(., subset = mito.percent <= 10)
SS2_allctypes %<>% subset(., subset = nCount_RNA <= 1*10^6 & 
                        nCount_RNA >= 1000
)

# SCTransform ====================================================================  
SS2_allctypes %<>% SCTransform(., assay = "RNA",
                           variable.features.n = NULL,
                           variable.features.rv.th = 1.8,
                           return.only.var.genes = F,
                           seed.use = 2,
                           do.scale = F
)
# PCA ====================================================================
SS2_allctypes %<>% RunPCA(., dims = 1:30, assay = "SCT")
DimPlot(SS2_allctypes, reduction = "pca", group.by = 'Region', pt.size = 1)

ElbowPlot(SS2_allctypes, ndims = 20)
DimHeatmap(SS2_allctypes, dims = c(1:2,7:10))

# UMAP ====================================================================
SS2_allctypes %<>% RunUMAP(.,
                       dims = 1:8, # 8 PCs for Lognorm and SCT
                       seed.use = 6,
                       n.neighbors = 50,
                       #local.connectivity = 2,
                       #n.epochs = 500,
                       min.dist = 0.6, 
                       #spread = 0.5,
                       assay = "SCT"
)

myPalette <- colorRampPalette(c(
  rep("midnightblue", 3), rep("blue2", 1),
  rep("grey", 1),
  rep("orange", 1), rep("red", 2)
))(256)

# Check expression of cell type markers
FeaturePlot(SS2_allctypes,
            features = c("Gimap4", "Ccnd2", "Mbp", "Mog", "Pdgfra", "Aldoc"), 
            pt.size = 0.75
) + coord_fixed(ratio = 1) + 
  scale_colour_gradientn(colours = myPalette)

# Clustering ====================================================================
SS2_allctypes %<>% FindNeighbors(., #features = VariableFeatures(SS2_allctypes)
                             reduction = 'pca', 
                             dims = 1:8
                             #prune.SNN = 1/10,
                             #k.param = 20
)

SS2_allctypes %<>% FindClusters(., algorithm = 2,
                            resolution = c(0.1, 0.2, 0.3, 0.4, 0.5,
                                           0.6, 0.8)
)

DimPlot(SS2_allctypes, reduction = 'umap', pt.size = 2,
        group.by = 'SCT_snn_res.0.6',
        #split.by = 'Genotype',
        label = T,
        label.size = 10
) + coord_fixed(ratio = 0.6)

# _ ----
# > Oligo Subset ####################################################################
# Select oligos based on canonical marker expression and clustering on the all cell types umap
c1 <- CellSelector(plot = DimPlot(SS2_allctypes))
oligo.subset <- SS2_allctypes[,c1]

oligo.subset$Age_Tissue <- paste0(oligo.subset$Age,'_',oligo.subset$Condition)

oligo.subset %<>% subset(., Mbp >= 6) # SCT, filter out non-oligos

# sctransform ====================================================================
oligo.subset %<>% SCTransform(.,
                              variable.features.n = 2000,
                              #variable.features.rv.th = 1.4,
                              return.only.var.genes = F,
                              #vars.to.regress = 'mito.percent',
                              do.scale = T,
                              seed.use = 3
)

p0 <- VariableFeaturePlot(oligo.subset)
top30 <- head(VariableFeatures(oligo.subset), 30)

LabelPoints(plot = p0, points = top30, repel = T)
# PCA ====================================================================
oligo.subset %<>% RunPCA(., dims = 1:20, assay = "SCT")

ElbowPlot(oligo.subset, ndims = 20)
DimHeatmap(oligo.subset, dims = c(1:2,3:6))

# UMAP ====================================================================
oligo.subset %<>% RunUMAP(.,
                          dims = 1:10, # 10 for sct, 7 for LogNorm (2K var) or 5 (1K var)
                          seed.use = 6,
                          n.neighbors = 50,
                          assay = "SCT"
)

myPalette <- colorRampPalette(c(
  rep("midnightblue", 3), rep("blue2", 1),
  rep("grey", 1),
  rep("orange", 1), rep("red", 2)
))

FeaturePlot(oligo.subset,
            features = 'Ifi27l2a',
            order = T,
            pt.size = 2.5
) + coord_fixed(ratio = 1) + 
  scale_colour_gradientn(colours = myPalette(256))