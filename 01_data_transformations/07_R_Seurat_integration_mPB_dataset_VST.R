library(Seurat)
library(future)
library(anndata)
library(batchelor)
library(scater)


plan(strategy = "multicore", workers = 12)
options(future.globals.maxSize = +Inf)


Sys.time()

data.path = 'source/'

infile <- paste0(data.path, 'h5ad/20211110_COMBO_PB10PLUS_clean.h5ad')
data = read_h5ad(infile)

obj <- SingleCellExperiment( assays = List(counts = as(t(data$X), "CsparseMatrix")),
                             colData=data$obs,
                             rowData=data$var
                           )

obj = as.Seurat(obj, data = 'counts')                           
                           
saveRDS(file = paste0(data.rds, 'rds/20211110_COMBO_PB10PLUS_filtered_counts_Seurat3_input_obj.rds'), obj)


data <- readRDS(file = paste0(data.path, 'rds/20211110_COMBO_PB10PLUS_filtered_counts_Seurat3_input_obj.rds') ) 


data.list <- SplitObject(object = data, split.by = "donor")


for (i in 1:length(x = data.list)) {
    data.list[[i]] <- NormalizeData(object = data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(object = data.list[[i]], 
        selection.method = "vst", nfeatures = 3000, verbose = FALSE)
}

reference.list <- data.list[ c("DOD1", "DOD2", "DOD3", "DOD4", "TQ198", "BP62j",
                               "BP37d", "BP74", "BP1c", "BP59h", "a_0", "a_1", "b_0", "b_1") ]

# find anchors                               
data.anchors <- readRDS(file = paste0(data.path, 'rds/20211110_COMBO_PB10PLUS_Seurat3_VST_classic_anchors.rds') )                            


# integrate features
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:15)

saveRDS(file = paste0(data.path, 'rds/20211110_COMBO_PB10PLUS_Seurat3_VST_classic_integrated.rds'), data.integrated)

# Export/save anndata
library(anndata)

#Making sure cell/gene orders are the same
data.integrated = data.integrated[row.names(data), colnames(data)]

#Convertin to AnnData and saving
x = data.integrated@assays$integrated@data
x = AnnData(X = t(x),
            obs = data.integrated@meta.data,
            var = data@assays$originalexp@meta.features)

            
outfile <-  paste0(data.path, 'h5ad/20211110_COMBO_PB10PLUS_Seurat3_VST_classic_integrated.h5ad')     


x$write_h5ad(outfile, compression = 'lzf')

Sys.time()
                            

 
