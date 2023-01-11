library(here)
library(tidyverse)
library(Seurat)

region_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG","V1")
subclass_order <- c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
                    "L5 ET", "L5/6 NP", "L6b", "L6 CT",
                    "Lamp5 Lhx6", "Lamp5", "Sncg", "Vip", "Pax6",
                    "Chandelier", "Pvalb", "Sst", "Sst Chodl",
                    "OPC", "Oligo", "Astro", "Micro/PVM", "Endo", "VLMC")

#load metadta
meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))

meta_sub <- meta %>%
  group_by(region, within_area_subclass) %>%
  sample_n(min(500, n())) %>%
  ungroup() %>%
  as.data.frame()

meta_sub <- as.data.frame(meta_sub)
rownames(meta_sub) <- meta_sub$sample_id

mat_list <- list()
regions <- meta %>% distinct(region) %>% pull()

#loop across all regions and store in list
for(j in 1:length(regions)){
  
  region_oi <- regions[j]
  
  #load matrices
  file_to_load <- list.files(here("data"), full.names = TRUE) %>%
    enframe() %>%
    filter(value %>% str_detect(region_oi)) 
  
  mat <- readRDS(file_to_load$value[1])
  mat <- mat[ , which(colnames(mat) %in% meta_sub$sample_id)]
  
  for(i in 2:nrow(file_to_load)){
    print(str_c("loading ", j, "_", i))
    mat_tmp <- readRDS(file_to_load$value[i])
    mat_tmp <- mat_tmp[ , which(colnames(mat_tmp) %in% meta_sub$sample_id)]
    mat <- cbind(mat, mat_tmp)
    gc()
  }
  mat_list[[j]] <- mat
  gc()
}

gc()
mat <- cbind(mat_list[[1]], mat_list[[2]], mat_list[[3]], mat_list[[4]],
             mat_list[[5]], mat_list[[6]], mat_list[[7]], mat_list[[8]])

rm(mat_list, mat_tmp)
gc()

#create seurat objects
seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta_sub)
rm(mat, meta, meta_sub)
gc()

seurat_obj <- NormalizeData(seurat_obj)


#load and filter DEGs
marker_loc <- str_c(here("markers"), "/")
markers <- list.files(marker_loc) %>%
  enframe() %>%
  filter(value %>% str_detect("vs_all")) %>%
  mutate(data = map(str_c(marker_loc, value), read.csv)) %>%
  unnest() %>%
  
  filter(p_val_adj < 0.05,
         avg_log2FC > 0.5) %>%
  
  mutate(region = value %>% str_remove_all("_subclass_vs_all_markers.csv"))
  
#wrangle DEGs 
conserved_counts <- markers %>%
  group_by(cluster, gene) %>%
  count() %>%
  ungroup() %>%
  
  filter(n == 8) %>%
  
  group_by(cluster) %>%
  count() %>%
  ungroup() %>%
  
  mutate(cluster = cluster %>% as_factor() %>% fct_relevel(subclass_order))

marker_counts <- markers %>%
  group_by(region, cluster) %>%
  count() %>%
  ungroup() %>%
  
  mutate(cluster = cluster %>% as_factor() %>% fct_relevel(subclass_order))

#plot markers
p1 <- ggplot() +
  geom_boxplot(data = marker_counts, aes(x = cluster, y = n)) +
  geom_point(data = conserved_counts, aes(x = cluster, y = n), color = "cyan") +
  scale_y_continuous(breaks = c(0, 250, 500, 750, 1000),
                     labels = c("0", "250", "500", "750", "1000")) +
  expand_limits(y = c(0, 1000)) +
  labs(x = "",
       y = "Number of subclass markers",
       title = "Subclass vs all other subclasses\nwithin region",
       subtitle = "Pval_adj < 0.05\nLog2FC > 0.5",
       caption = "**Cyan point denote number of markers conserved across 8 regions") +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))

pdf(file = here("plots", "boxplot_conserved_subclass_DEGs.pdf"),
    height = 5,
    width = 5,
    useDingbats = F)
print(p1)
dev.off()  

#plot heatmap of conserved genes
genes_to_plot <- markers %>%
  group_by(gene, cluster) %>%
  mutate(gene_count = n(),
         mean_avg_log2FC = mean(avg_log2FC)) %>%
  ungroup() %>%
  
  filter(gene_count == 8) %>%
  distinct(cluster, gene, mean_avg_log2FC) %>%
  
  mutate(cluster = cluster %>% as_factor() %>% fct_relevel(subclass_order)) %>%
  arrange(cluster, desc(mean_avg_log2FC)) %>%
  distinct(gene) %>%
  pull(gene)


#downsample nuclei to 50 per region per subclass
seurat_obj$ds_label <- str_c(seurat_obj$region, "_", seurat_obj$within_area_subclass)
Idents(seurat_obj) <- seurat_obj$ds_label
seurat_obj <- ScaleData(seurat_obj, features = genes_to_plot)
seurat_obj_ds <- subset(seurat_obj, downsample = 50)

#subset to a single subclass
Idents(seurat_obj_ds) <- seurat_obj_ds$within_area_subclass

for(i in 1:length(subclass_order)){
  
subclass_oi <- subclass_order[i]
seurat_obj_ds_oi <- subset(seurat_obj_ds, idents = subclass_oi)

Idents(seurat_obj_ds_oi) <- seurat_obj_ds_oi$region
levels(seurat_obj_ds_oi) <- region_order
color_use <- unique(seurat_obj_ds_oi@meta.data$within_area_subclass_color)

p1 <- DoHeatmap(seurat_obj_ds_oi, slot = "scale.data", features = genes_to_plot, raster = T, label = F) +
  scale_fill_gradient2(low = "white", mid = "white", high = color_use) +
  ggtitle(subclass_oi) +
  theme(aspect.ratio = 6,
    legend.position = "none",
    axis.text.y = element_blank())

pdf(file = here("plots", "subclass_heatmaps", str_c("heatmap_", subclass_oi %>% 
                                                      str_remove_all(" ") %>%
                                                      str_remove_all("/") %>%
                                                      str_remove_all("-") %>%
                                                      str_remove_all("_"), ".pdf")),
    height = 6,
    width = 5,
    useDingbats = F)
print(p1)
dev.off()  
}
