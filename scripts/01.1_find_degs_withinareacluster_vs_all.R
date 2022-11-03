library(Seurat)
library(tidyverse)
library(here)

#load metadta

meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))
regions <- meta %>% distinct(region) %>% pull()
region_oi <- regions[8]

meta_sub <- meta %>%
  filter(region == region_oi) %>%
  group_by(within_area_cluster) %>%
  sample_n(min(500, n())) %>%
  ungroup() %>%
  as.data.frame()

meta_sub <- as.data.frame(meta_sub)
rownames(meta_sub) <- meta_sub$sample_id

#load matrices
file_to_load <- list.files(here("data"), full.names = TRUE) %>%
  enframe() %>%
  filter(value %>% str_detect(region_oi)) 

mat <- readRDS(file_to_load$value[1])
mat <- mat[ , which(colnames(mat) %in% meta_sub$sample_id)]

for(i in 2:nrow(file_to_load)){
  print(str_c("loading ", i))
  mat_tmp <- readRDS(file_to_load$value[i])
  mat_tmp <- mat_tmp[ , which(colnames(mat_tmp) %in% meta_sub$sample_id)]
  mat <- cbind(mat, mat_tmp)
  gc()
}

dim(mat)  
gc()

#create seurat objects
seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta_sub)
seurat_obj <- NormalizeData(seurat_obj)

Idents(seurat_obj) <- seurat_obj$within_area_cluster

markers <- FindAllMarkers(seurat_obj, assay = "RNA", slot = "data", 
                          test.use = "wilcox", only.pos = TRUE) 


write.csv(markers, here("markers", str_c(region_oi, "_cluster_vs_all_markers.csv")))
gc()
