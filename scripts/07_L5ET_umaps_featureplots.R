library(tidyverse)
library(Seurat)
library(here)
library(scattermore)

#load reductions and metadata and combine to plot L5 ET UMAPs
meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))        

reductions <- readRDS(here("umap_coords", "Neighborhood_UMAPS", "Deep_Exc_reductions.RDS"))

to_plot <- reductions$umap@cell.embeddings %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  filter(within_area_subclass == "L5 ET")

#make each within_region_cluster color a shade of the region color
(color_lookup <- to_plot %>% distinct(region, region_color, within_area_cluster))


color_oi <- color_lookup %>% filter(region == "A1") %>% distinct(region_color) %>% pull(region_color)
colfunc_A1 <- colorRampPalette(c(color_oi, "white"))
color_oi <- color_lookup %>% filter(region == "ACC") %>% distinct(region_color) %>% pull(region_color)
colfunc_ACC <- colorRampPalette(c(color_oi, "white"))
color_oi <- color_lookup %>% filter(region == "DLPFC") %>% distinct(region_color) %>% pull(region_color)
colfunc_DLPFC <- colorRampPalette(c(color_oi, "white"))
color_oi <- color_lookup %>% filter(region == "M1") %>% distinct(region_color) %>% pull(region_color)
colfunc_M1 <- colorRampPalette(c(color_oi, "white"))
color_oi <- color_lookup %>% filter(region == "MTG") %>% distinct(region_color) %>% pull(region_color)
colfunc_MTG <- colorRampPalette(c(color_oi, "white"))
color_oi <- color_lookup %>% filter(region == "S1") %>% distinct(region_color) %>% pull(region_color)
colfunc_S1 <- colorRampPalette(c(color_oi, "white"))
color_oi <- color_lookup %>% filter(region == "V1") %>% distinct(region_color) %>% pull(region_color)
colfunc_V1 <- colorRampPalette(c(color_oi, "white"))
color_oi <- color_lookup %>% filter(region == "ANG") %>% distinct(region_color) %>% pull(region_color)
colfunc_ANG <- colorRampPalette(c(color_oi, "white"))

color_lookup
color_lookup$cluster_color <- "none"
color_lookup$cluster_color[1] <- colfunc_A1(10)[1]
color_lookup$cluster_color[2] <- colfunc_A1(10)[5]
color_lookup$cluster_color[3] <- colfunc_ACC(10)[1]
color_lookup$cluster_color[4] <- colfunc_DLPFC(10)[1]
color_lookup$cluster_color[5] <- colfunc_DLPFC(10)[5]
color_lookup$cluster_color[6] <- colfunc_M1(10)[1]
color_lookup$cluster_color[7] <- colfunc_M1(10)[5]
color_lookup$cluster_color[8] <- colfunc_M1(10)[7]
color_lookup$cluster_color[9] <- colfunc_MTG(10)[1]
color_lookup$cluster_color[10] <- colfunc_MTG(10)[5]
color_lookup$cluster_color[11] <- colfunc_S1(10)[1]
color_lookup$cluster_color[12] <- colfunc_S1(10)[5]
color_lookup$cluster_color[13] <- colfunc_V1(10)[1]
color_lookup$cluster_color[14] <- colfunc_ANG(10)[1]
color_lookup$cluster_color[15] <- colfunc_ANG(10)[5]

color_lookup <- color_lookup %>% 
  mutate(region_cluster = str_c(region, "_", within_area_cluster))

region_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG","V1")

#make plot
colors_use <- color_lookup %>% select(region, region_color) %>% deframe()
p1 <- to_plot %>%
  mutate(region_cluster = str_c(region, "_", within_area_cluster),
         region = region %>% as_factor() %>% fct_relevel(region_order)) %>%
  filter(UMAP_1 < 0) %>%
  
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = region)) +
  scale_color_manual(values = colors_use) +
  theme_void() +
  theme(aspect.ratio = 1)

colors_use <- color_lookup %>% select(region_cluster, cluster_color) %>% deframe()
p2 <- to_plot %>%
  mutate(region_cluster = str_c(region, "_", within_area_cluster),
         region = region %>% as_factor() %>% fct_relevel(region_order)) %>%
  filter(UMAP_1 < 0) %>%
  
  ggplot() +
  geom_point(aes(x = UMAP_1, y = UMAP_2, color = region_cluster)) +
  facet_wrap(~region, ncol = 4) +
  scale_color_manual(values = colors_use) +
  labs(title = "within-region clusters") +
  theme_void() +
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = "white", color = "black"))

p1 + p2 


pdf(file = here("plots", "L5_ET_umap.pdf"),
    height = 5,
    width = 15,
    useDingbats = F)
print(p1 + p2)
dev.off()

#load gene expression for features plots
file_to_load <- list.files(here("data"), full.names = T) %>%
  enframe() %>%
  filter(value %>% str_detect("deep_exc"))


mat <- readRDS(file_to_load$value[1])
mat <- mat[ , which(colnames(mat) %in% to_plot$sample_id)]

for(i in 2:nrow(file_to_load)){
  print(str_c("loading ", i))
  mat_tmp <- readRDS(file_to_load$value[i])
  mat_tmp <- mat_tmp[ , which(colnames(mat_tmp) %in% to_plot$sample_id)]
  mat <- cbind(mat, mat_tmp)
  gc()
}

dim(mat)  
gc()

#create seurat object
meta_sub <- to_plot %>% as.data.frame()
rownames(meta_sub) <- meta_sub$sample_id

seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta_sub)
seurat_obj <- NormalizeData(seurat_obj)

#find markers of region clusters
seurat_obj$region_cluster <- str_c(seurat_obj$region, "_", seurat_obj$within_area_cluster)
Idents(seurat_obj) <- seurat_obj$region_cluster
markers <- FindAllMarkers(seurat_obj, slot = "data", test.use = "wilcox", only.pos = T)
write.csv(markers, here("markers", "within_L5ET_markers.csv"))

genes_oi <- markers %>%
  mutate(prop_diff = pct.1 - pct.2) %>%
  filter(p_val_adj < 0.05,
         avg_log2FC > 0.5,
         prop_diff > 0.5) %>% 
  
  # group_by(cluster) %>%
  # mutate(genes_oi = first(gene)) %>%
  # ungroup() %>%
  # distinct(genes_oi) %>%
  distinct(gene) %>%
  pull()

cell_include <- to_plot %>% filter(UMAP_1 < 0) %>% pull(sample_id)

seurat_obj@reductions <- reductions
FeaturePlot(seurat_obj, slot = "data", features = genes_oi[1:20], ncol = 4, raster = T)
FeaturePlot(seurat_obj, slot = "data", features = genes_oi[21:40], ncol = 4, raster = T)
FeaturePlot(seurat_obj, slot = "data", features = genes_oi[41:60], ncol = 4, raster = T)
FeaturePlot(seurat_obj, slot = "data", features = genes_oi[61:79], ncol = 4, raster = T)

curated_genes <- c("KCNG2", "MEIS3", "PROM1", "LINC00922", "DCN",
                   "LINC00836", "EYA4", "FREM2", "HTR1F", "HTR2A", "HTR2C",
                   "CALCRL", "ONECUT1", "PARD3", "ALDH1A2", "DOK5")
p1 <- FeaturePlot(seurat_obj, slot = "data", features = curated_genes, ncol = 4, raster = T, cells = cell_include, coord.fixed = 1) 

pdf(file = here("plots", "L5_ET_featureplot_full.pdf"),
    height = 15,
    width = 15,
    useDingbats = F)
print(p1)
dev.off()


curated_genes <- c("KCNG2", "PROM1", "DCN",
                   "LINC00836", "EYA4", "FREM2", "PARD3","DOK5")
p1 <- FeaturePlot(seurat_obj, slot = "data", features = curated_genes, ncol = 4, raster = T, cells = cell_include, coord.fixed = 1) 

pdf(file = here("plots", "L5_ET_featureplot_targeted.pdf"),
    height = 15,
    width = 15,
    useDingbats = F)
print(p1)
dev.off()
