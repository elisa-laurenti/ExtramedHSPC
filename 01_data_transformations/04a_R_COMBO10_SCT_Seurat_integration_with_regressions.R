library(Seurat)
library(future)

plan(strategy = "multicore", workers = 10)

options(future.globals.maxSize = +Inf) # (1024 * 32) * 1024^2 )



Sys.time()


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


data.path = '~/databoard/users/2021/BloodPaper/rds/'

data <- readRDS(file = paste0(data.path, '20210909_COMBO10_XTRA_filtered_counts_Seurat3_obj.rds') )


data.list <- SplitObject(object = data, split.by = "donor")

for (i in 1:length(x = data.list)) {

    # normalize data with SCTransform()
    # ----------------------------------------------------------------------------------
    data.list[[i]] <- SCTransform(object = data.list[[i]],
                      assay = 'RNA',
                      new.assay.name = 'SCT',
                      vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent_mitoc"),
                      verbose = TRUE)

    # perform cell cycle analysis (make sure to specify the "assay" parameter
    # ----------------------------------------------------------------------------------
    data.list[[i]] <- CellCycleScoring( data.list[[i]],
      s.features = s.genes,
      g2m.features = g2m.genes,
      assay = 'SCT',
      set.ident = TRUE,
      verbose = TRUE
    )

    # normalise again but this time including also the cell cycle scores
    # ----------------------------------------------------------------------------------
    sample <- SCTransform( data.list[[i]],
      assay = 'RNA',
      new.assay.name = 'SCT',
      vars.to.regress = c('percent_mitoc', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score'),
      verbose = TRUE
    )

}


features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)

data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)


reference.list <- data.list[ c("DOD1", "DOD2", "DOD3", "DOD4", "TQ198",
                               "BP62j", "BP37d", "BP74", "BP1c","BP59h") ]


data.anchors <- FindIntegrationAnchors(object.list = reference.list,
                                       normalization.method = "SCT",
                                       anchor.features = features,
                                       dims = 1:15, scale = FALSE)


saveRDS(file = paste0(data.path, '20210910_COMBO10_NO_SPL3_Seurat3_SCT_and_REGRESS_anchors.rds'), data.anchors)


data.integrated <- IntegrateData(anchorset = data.anchors, , normalization.method = "SCT", dims = 1:15)


data.integrated <- RunPCA(object = data.integrated, npcs = 30, verbose = FALSE)
data.integrated <- RunUMAP(object = data.integrated, reduction = "pca", dims = 1:30)


saveRDS(file = paste0(data.path, '20210910_COMBO10_NO_SPL3_Seurat3_integrated_SCT_and_REGRESS.rds'), data.integrated)


Sys.time()
 
