library(tidyverse)
library(Seurat)
library(here)
library(scattermore)
library(ggrepel)


#load reductions and metadata and combine to plot L5 ET UMAPs
meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))   
meta <- meta %>% 
  filter(region == "V1",
         within_area_subclass == "Sst") %>%
  as.data.frame()
rownames(meta) <- meta$sample_id

#load gene expression for features plots
file_to_load <- list.files(here("data"), full.names = T) %>%
  enframe() %>%
  filter(value %>% str_detect("mge_inh"),
         value %>% str_detect("V1"))


mat <- readRDS(file_to_load$value[1])
mat <- mat[ , which(colnames(mat) %in% meta$sample_id)]

dim(mat)  
gc()

#create seurat obj
seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta)
seurat_obj <- NormalizeData(seurat_obj)

Idents(seurat_obj) <- seurat_obj$within_area_cluster
markers <- FindAllMarkers(seurat_obj, slot = "data", test.use = "wilcox", only.pos = T, max.cells.per.ident = 500)

write.csv(markers, here("markers", "V1_SST_subtype_markers.csv"))
markers <- read.csv(here("markers", "V1_SST_subtype_markers.csv"))

cluster_order <- str_c("Sst_", c(27:22, 16, 15, 18, 14, 21, 20, 17, 19, 13, 11, 12, 8, 7, 9, 10, 6:3, 1, 2))



genes_use <- markers %>%
  as_tibble() %>%
  mutate(prop_diff = pct.1 - pct.2,
         cluster = cluster %>% as.factor() %>% fct_relevel(cluster_order)) %>%
  filter(p_val_adj < 0.05,
         avg_log2FC > 0.5,
         prop_diff > 0.4) %>% 
  arrange(cluster) %>%
  
  group_by(cluster) %>%
  slice(1:10) %>%
  ungroup() %>%
  distinct(gene) %>%
  pull()

#plot heatmap
table(seurat_obj$within_area_cluster)
seurat_obj_sub <- subset(seurat_obj, downsample = 50)
seurat_obj_sub <- ScaleData(seurat_obj_sub, features = genes_use)
levels(seurat_obj_sub) <- cluster_order

p1 <- DoHeatmap(seurat_obj_sub, slot = "scale.data", features = genes_use, raster = T, size = 2) +
  scale_fill_viridis_c(option = "B") +
  labs(caption = "p_val_adj < 0.05\nlog2FC > 0.5\nprop.diff > 0.4",
       fill = "Scaled Expression") +
  theme(aspect.ratio = 3,
        axis.text.y = element_text(size = 5)) 

pdf(file = here("plots", "SST_subtype_markers.pdf"),
    height = 8,
    width = 6,
    useDingbats = F)
print(p1)
dev.off()

#rearrange by tx similarity instead of cluster id number
sst_means <- seurat_obj_sub@assays$RNA@scale.data %>%
  as_tibble(rownames = "gene") %>%
  gather(key = "sample_id", value = "expr", -gene) %>%
  left_join(seurat_obj_sub@meta.data, by = "sample_id") %>%
  
  group_by(gene, within_area_cluster) %>%
  summarise(mean_expr = mean(expr)) %>%
  ungroup() %>%
  
  spread(key = within_area_cluster, value = mean_expr)

sst_means_mat <- as.matrix(sst_means[ , -1])
rownames(sst_means_mat) <- sst_means$gene

pheatmap::pheatmap(sst_means_mat)

new_order <- c("Sst_1", "Sst_3", "Sst_4", "Sst_5", "Sst_2", "Sst_6",
               "Sst_14", "Sst_18", "Sst_19", "Sst_15", "Sst_16", "Sst_17", "Sst_20", "Sst_12",
               "Sst_21", "Sst_26", "Sst_27", "Sst_25", "Sst_24", "Sst_23", "Sst_22", 
               "Sst_11", "Sst_13", "Sst_10", "Sst_9", "Sst_7", "Sst_8")


#rerun above with new order
genes_use <- markers %>%
  as_tibble() %>%
  mutate(prop_diff = pct.1 - pct.2,
         cluster = cluster %>% as.factor() %>% fct_relevel(new_order)) %>%
  filter(p_val_adj < 0.05,
         avg_log2FC > 0.5,
         prop_diff > 0.4) %>% 
  arrange(cluster) %>%
  
  group_by(cluster) %>%
  slice(1:10) %>%
  ungroup() %>%
  distinct(gene) %>%
  pull()

#plot heatmap
table(seurat_obj$within_area_cluster)
seurat_obj_sub <- subset(seurat_obj, downsample = 50)
seurat_obj_sub <- ScaleData(seurat_obj_sub, features = genes_use)
levels(seurat_obj_sub) <- new_order

DoHeatmap(seurat_obj_sub, slot = "scale.data", features = genes_use, raster = T, size = 2) +
  scale_fill_viridis_c(option = "B") +
  labs(caption = "p_val_adj < 0.05\nlog2FC > 0.5\nprop.diff > 0.4",
       fill = "Scaled Expression") +
  theme(aspect.ratio = 3,
        axis.text.y = element_text(size = 5)) 

