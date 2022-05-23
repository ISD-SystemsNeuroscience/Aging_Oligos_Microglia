source("Loadings.R")
# > Integration Microglia RagKO + WT (with new libraries) ####################################################################
mg.WT_2 <- readRDS("../20200907_10x/dataset_3_20200929.rds")
mg.WT_2$Region <- ifelse(mg.WT_2$Tissue == "GM", "GreyM", "WhiteM")
mg.WT_2$Batch <- "Batch3"
mg.WT_2$Library <- paste0(mg.WT_2$Batch,'_',mg.WT_2$Region)

mg.WT_1 <- readRDS("../10X_WMGM_1st/Training_Sobj.microglia.rds")
mg.WT_1$Region <- mg.WT_1$orig.ident
mg.WT_1$Batch <- "Training"
mg.WT_1$Library <- paste0(mg.WT_1$Batch,'_',mg.WT_1$Region)

mg.RagKO_1 <- readRDS("Sobj_mgfiltered_RagKO_Mar252021.rds")
mg.RagKO_1$Batch <- "RagKO_2020"
mg.RagKO_1$Library <- paste0(mg.RagKO_1$Batch,'_',mg.RagKO_1$Region)

mg.RagKO_2$Batch <- "RagKO_2022"
mg.RagKO_2$Library <- paste0(mg.RagKO_2$orig.ident,'_',
                                      mg.RagKO_2$Region)

mg.list <- list(mg.WT_1, mg.WT_2, mg.RagKO_1, mg.RagKO_2)

com.genes <- SelectIntegrationFeatures(object.list = mg.list,
                                       assay = rep("SCT",4),
                                       nfeatures = 500
)

mg.list <- PrepSCTIntegration(object.list = mg.list,
                              anchor.features = com.genes,
                              verbose = TRUE)
mg.list <- lapply(X = mg.list, FUN = function(x) {
  x <- RunPCA(x, features = com.genes, 
              npcs = 30, 
              verbose = T)
})
mg.anchors <- FindIntegrationAnchors(object.list = mg.list,
                                     dims = 1:15,
                                     k.anchor = 4,
                                     reduction = "rpca",
                                     max.features = 100,
                                     normalization.method = "SCT",
                                     anchor.features = com.genes, verbose = TRUE)

#st <- matrix(c(-1,1,2,-2,-3,-4), ncol = 2)
mg.merge <- IntegrateData(anchorset = mg.anchors,
                          dims = 1:15,
                          normalization.method = "SCT",
                          #sample.tree = st,
                          verbose = TRUE)

DefaultAssay(mg.merge) <- "RNA"
mg.merge %<>% NormalizeData(.)
DefaultAssay(mg.merge) <- "integrated"
mg.merge$Genotype <- ifelse(is.na(mg.merge$Genotype), "WT", mg.merge$Genotype)

mg.merge$Genotype_Region <- paste0(mg.merge$Genotype,"_",mg.merge$Region)