library(Seurat)
library(scDiffCom)
library(dplyr)
library(future)

options(future.globals.maxSize = 10 * 1024 * 1024^2)

pbmc.all <- readRDS(file.path('data','seurat','pbmc_all.rds'))
pbmc.all$celltype <- ifelse(pbmc.all$predicted.celltype.l1 %in% c('B','CD4 T','CD8 T'),
                            pbmc.all$predicted.celltype.l2,
                            pbmc.all$predicted.celltype.l1)

plan(multicore, workers = 4)

scdiffcom_object <- run_interaction_analysis(
  seurat_object = pbmc.all,
  LRI_species = "human",
  seurat_celltype_id = "celltype",
  iterations=10000,
  seurat_condition_id = list(
    column_name = "donor_type",
    cond1_name = "donor",
    cond2_name = "recipient"
  )
)

saveRDS(scdiffcom_object, file.path('data','scDiffCom','scdiffcom_object.rds'))
