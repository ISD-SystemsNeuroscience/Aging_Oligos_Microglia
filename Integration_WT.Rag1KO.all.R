source("Loadings.R")

# > Integration all ctypes WT (FACS) + WT (Myelin clean up) + Rag1KO ####################################################################
dataset.RagKO <- readRDS("../RagKO_20210120/Sobj_allRagKOLibraries_allctypes_sctransform_7Apr2022.rds")
data.FACS <- readRDS("Sobj_SortedIntegrated_allctypes_1stand2ndRound.rds")
data.mye.clean <- readRDS("Sobj.Training_Batch3_Batch1GM_all.Integration_2020Nov27_updated.rds")
data.mye.clean$Age <- ifelse(is.na(data.mye.clean$Age), '24mo', 
                                  data.mye.clean$Age)

DefaultAssay(dataset.RagKO) <- "RNA"
DefaultAssay(data.FACS) <- "RNA"
DefaultAssay(data.mye.clean) <- "RNA"

rag.list <- SplitObject(dataset.RagKO, split.by = 'orig.ident')
wt.facs <- SplitObject(data.FACS, split.by = 'orig.ident')
wt.myelin <- SplitObject(data.mye.clean, split.by = 'orig.ident')

integrated.all.list <- list(wt.myelin[[1]],wt.myelin[[2]],wt.myelin[[3]],
                            wt.myelin[[4]], wt.facs[[1]], wt.facs[[2]],
                            wt.facs[[3]], wt.facs[[4]], rag.list[[1]],
                            rag.list[[2]], rag.list[[3]], rag.list[[4]]
)
###
integrated.all.list <- lapply(X = integrated.all.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = integrated.all.list)
integrated.all.list <- lapply(X = integrated.all.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

int.anchors <- FindIntegrationAnchors(object.list = integrated.all.list, 
                                      anchor.features = features, 
                                      reduction = "rpca")


wt.rag.all <- IntegrateData(anchorset = int.anchors)

# QC ====================================================================
VlnPlot(wt.rag.all, 
        features = c("nFeature_RNA", "mito.percent"),
        #features = 'nCount_RNA',
        group.by = 'Batch', ncol = 1,
        #cols = c(isdblue,isdorange),
        cols = brewer.pal(12,"Paired"),
        pt.size = 0
)

wt.rag.all$Genotype <- ifelse(wt.rag.all$orig.ident == "Sorted_14",
                              "WT", wt.rag.all$Genotype)
wt.rag.all$Genotype <- ifelse(wt.rag.all$orig.ident == "Sorted_15",
                              "WT", wt.rag.all$Genotype)
wt.rag.all$Genotype <- ifelse(wt.rag.all$orig.ident == "WhiteM",
                              "WT", wt.rag.all$Genotype)
wt.rag.all$Genotype <- ifelse(wt.rag.all$orig.ident == "GreyM",
                              "WT", wt.rag.all$Genotype)

wt.rag.all$Genotype_Tissue <- paste0(wt.rag.all$Genotype,'_',wt.rag.all$Region)

# Tidy up method and batch metadata
method.list <- list(
  "L1" = "FACS",
  "L2" = "FACS",
  "Library_05" = "MyelinClean",
  "Library_06" = "MyelinClean",
  "Library_16.Sorted" = "FACS",
  "Library_17.Sorted" = "FACS",
  "R1WM_D1" = "FACS",
  "R2WM_D2" = "FACS",
  "Sorted_14" = "FACS",
  "Sorted_15" = "FACS",
  "WhiteM" = "MyelinClean",
  "GreyM" = "MyelinClean"
)
wt.rag.all$Method <- unlist(method.list[wt.rag.all$orig.ident])

batch.list <- list(
  "L1" = "Rag1KO_1st_GM",
  "L2" = "Rag1KO_1st_WM",
  "Library_05" = "MyeCl_2nd_GM",
  "Library_06" = "MyeCl_2nd_WM",
  "Library_16.Sorted" = "FACS_2nd_GM",
  "Library_17.Sorted" = "FACS_2nd_WM",
  "R1WM_D1" = "Rag1KO_2nd_WM_1",
  "R2WM_D2" = "Rag1KO_2nd_WM_2",
  "Sorted_14" = "FACS_1st_GM",
  "Sorted_15" = "FACS_1st_WM",
  "WhiteM" = "MyeCl_1st_WM",
  "GreyM" = "MyeCl_1st_GM"
)

wt.rag.all$Batch <- unlist(batch.list[wt.rag.all$orig.ident])
# PCA ====================================================================
wt.rag.all %<>% ScaleData(., verbose = T,
                          vars.to.regress = 'mito.percent'
)
wt.rag.all %<>% RunPCA(., dims = 1:30, assay = 'integrated')
# UMAP ====================================================================
wt.rag.all %<>% RunUMAP(.,
                        #features = top_n_markers$gene,
                        dims = 1:30,
                        seed.use = 6,
                        n.neighbors = 50,
                        #local.connectivity = 100,
                        #n.epochs = 500,
                        #spread = 0.5,
                        #min.dist = 0.15,
                        #reduction.name = "umap_topgenes",
                        assay = "integrated"
)

# Check marker gene expression
myPalette <- colorRampPalette(c(
  rep("midnightblue", 3), rep("blue2", 1),
  rep("grey", 1),
  rep("orange", 1), rep("red", 3)
))(256)
DefaultAssay(wt.rag.all) <- "RNA"
FeaturePlot(wt.rag.all, 
            #features = 'Nkg7', 
            features = c('Ctss','Plp1','Slc1a2','Nkg7'),
            #slot = 'data',
            order = F,
            pt.size = 0.4
            #label = T, 
) & scale_color_gradientn(colors = myPalette) & coord_fixed(ratio = 0.65)
DefaultAssay(wt.rag.all) <- "integrated"
# Clustering ====================================================================
wt.rag.all %<>% FindNeighbors(., #features = VariableFeatures(wt.rag.all)
                              reduction = 'pca', 
                              dims = 1:30,
                              assay = "integrated"
                              #slot = 'scale.data'
                              #prune.SNN = 1/10,
                              #k.param = 20
)

wt.rag.all %<>% FindClusters(., algorithm = 2,
                             resolution = c(0.4, 0.5, 0.6))


clustree(wt.rag.all, prefix = 'integrated_snn_res.')

DimPlot(wt.rag.all, reduction = 'umap',
        pt.size = 0.6,
        group.by = "integrated_snn_res.0.5",
        label = T, label.size = 5
) + coord_fixed(ratio = 0.8)


# Find Markers ====================================================================
Idents(wt.rag.all) <- "integrated_snn_res.0.5"
markers_wt.rag.all_res0.5 <- FindAllMarkers(wt.rag.all,
                                            assay = "RNA",
                                            test.use = "MAST",
                                            #min.pct = 0.25,
                                            only.pos = T,
                                            logfc.threshold = 1.5
) %>% filter(p_val_adj <= 0.05) %>% arrange(cluster, desc(avg_log2FC))

# Cell Type Annotation ====================================================================
ctypeannot <- list(
  "0" = "Microglia",
  "1" = "Oligodendrocytes",
  "2" = "Oligodendrocytes",
  "3" = "Microglia",
  "4" = "Microglia",
  "5" = "Oligodendrocytes",
  "6" = "Astrocytes",
  "7" = "Microglia",
  "8" = "T/NK Cells",
  "9" = "Oligodendrocytes",
  "10" = "Ependymal Secretory",
  "11" = "Astrocytes",
  "12" = "Astrocytes",
  "13" = "ECs",
  "14" = "Neurons/NSCs",
  "15" = "Ependymal Cells",
  "16" = "Macrophages",
  "17" = "Ependymal Cilia",
  "18" = "Erythrocyte-like",
  "19" = "19",
  "20" = "Neuroblasts",
  "21" = "Pericytes",
  "22" = "Microglia",
  "23" = "23",
  "24" = "OPCs",
  "25" = "Neutrophils",
  "26" = "Fibroblasts",
  "27" = "Oligodendrocytes",
  "28" = "Neutrophils",
  "29" = "29"
)
wt.rag.all$integrated_snn_res.0.5 %<>% as.character()
wt.rag.all$CellTypeAnnot <- unlist(ctypeannot[wt.rag.all$integrated_snn_res.0.5])

wt.rag.all.facs <- subset(wt.rag.all, Method == "FACS")

myPalette2 <- colorRampPalette(c(
  rep("midnightblue", 4), rep("blue2", 1),
  rep("grey", 1),
  rep("orange", 1), rep("red", 2)
))(256)

DefaultAssay(wt.rag.all) <- "RNA"
DotPlot(wt.rag.all,
        features = (c('Slc1a2','Nkg7','Hexb','Vtn','Sox4','Mbp','Lyz2','S100a9')),
        group.by = "CellTypeAnnot",
        dot.scale = 20
) + scale_color_gradientn(colours = myPalette2)
