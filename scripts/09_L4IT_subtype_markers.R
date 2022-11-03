library(tidyverse)
library(Seurat)
library(here)
library(scattermore)
library(ggrepel)

#load reductions and metadata and combine to plot L5 ET UMAPs
meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))   
meta <- meta %>% 
  filter(region == "V1",
         within_area_subclass == "L4 IT") %>%
  as.data.frame()
rownames(meta) <- meta$sample_id

#load gene expression for features plots
file_to_load <- list.files(here("data"), full.names = T) %>%
  enframe() %>%
  filter(value %>% str_detect("it_type"),
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
write.csv(markers, here("markers", "V1_L4IT_subtype_markers.csv"))
markers <- read.csv(here("markers", "V1_L4IT_subtype_markers.csv"))

cluster_order <- str_c("L4 IT_", c(8, 7, 1:6))

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
p1
pdf(file = here("plots", "L4IT_subtype_markers.pdf"),
    height = 8,
    width = 6,
    useDingbats = F)
print(p1)
dev.off()
