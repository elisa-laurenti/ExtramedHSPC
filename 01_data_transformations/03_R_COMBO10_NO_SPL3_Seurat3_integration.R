 
library(Seurat)


Sys.time()


data.path = '~/datafloor/users/2020/SLX19841/analysis/matrices/'

data <- readRDS(file = paste0(data.path, 'COMBO10_NO_SPL3_filtered_counts_Seurat3_obj.rds') )


data.list <- SplitObject(object = data, split.by = "donor")

for (i in 1:length(x = data.list)) {
    data.list[[i]] <- NormalizeData(object = data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(object = data.list[[i]], 
        selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}


reference.list <- data.list[ c("DOD1", "DOD2", "DOD3", "DOD4", "TQ198", 
                               "BP62j", "BP37d", "BP74", "BP1c","BP59h") ]
                               
                               
Sys.time()
data.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:15)
Sys.time()                            
                               
                               
saveRDS(file = paste0(data.path, 'COMBO10_NO_SPL3_Seurat3_anchors.rds'), data.anchors)                            
                               

Sys.time()
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:15)
Sys.time()


# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = data.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
data.integrated <- ScaleData(object = data.integrated, verbose = FALSE)
data.integrated <- RunPCA(object = data.integrated, npcs = 30, verbose = FALSE)
data.integrated <- RunUMAP(object = data.integrated, reduction = "pca", 
    dims = 1:30)
    
    
saveRDS(file = paste0(data.path, 'COMBO10_NO_SPL3_Seurat3_integrated.rds'), data.integrated)


Sys.time()
