library(tidyverse)
library(Seurat)
library(here)
library(scattermore)
library(ggrepel)

#load reductions and metadata and combine to plot L5 ET UMAPs
meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))        

reductions <- readRDS(here("umap_coords", "Neighborhood_UMAPS", "Glia_reductions.RDS"))

to_plot <- reductions$umap@cell.embeddings %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(meta, by = "sample_id") 


region_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG","V1")

#load gene expression for features plots
file_to_load <- list.files(here("data"), full.names = T) %>%
  enframe() %>%
  filter(value %>% str_detect("glia"))


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


calc_tau <- function(vec) {
  xi_vec = vec/max(vec)
  tau = sum(1 - xi_vec) / (length(vec) -1 )
  return(tau)
}

subclass_order <- c("OPC", "Oligo", "Astro", "Micro/PVM", "Endo", "VLMC")

genes_oi <- list.files(here("markers"), full.names = T) %>%
  enframe() %>%
  filter(value %>% str_detect("vs_all")) %>%
  mutate(data = map(value, read.csv)) %>%
  unnest() %>%
  
  mutate(prop_diff = pct.1 - pct.2) %>%
  filter(avg_log2FC > 1,
         p_val_adj < 0.05,
         prop_diff > 0.5,
         cluster %in% to_plot$within_area_subclass) %>%
  group_by(cluster, gene) %>%
  count() %>%
  ungroup() %>%
  filter(n == 8) %>%
  distinct(cluster, gene) %>%
  
  group_by(cluster) %>%
  slice(1:20) %>%
  ungroup() %>%
  
  mutate(cluster = cluster %>% as_factor() %>% fct_relevel(subclass_order)) %>%
  
  arrange(cluster) %>%
  distinct(gene) %>%
  pull(gene)



dot_plot_dat <- seurat_obj@assays$RNA@data[genes_oi, ] %>%
  as_tibble(rownames = "gene") %>%
  gather(key = "sample_id", value = "expr", -gene)

to_plot_dot <- dot_plot_dat %>%
  left_join(seurat_obj@meta.data, by = "sample_id") %>%
  mutate(expr_thresh = case_when(expr > 0 ~ 1,
                                 TRUE ~ 0)) 
  
to_plot_dot <- to_plot_dot %>%  
  group_by(gene, within_area_subclass, region) %>%
  mutate(total = n(),
         thresh_count = sum(expr_thresh),
         mean_expr = mean(expr)) %>%
  ungroup() %>%
  
  group_by(gene) %>%
  mutate(tau = calc_tau(mean_expr)) %>%
  ungroup() %>%
  
  distinct(gene, mean_expr, region, within_area_subclass, total, thresh_count, tau) %>%
  mutate(prop_count = thresh_count / total) %>%
 
  group_by(gene) %>%
  mutate(tau_prop = calc_tau(prop_count)) %>%
  ungroup()

genes_oi2 <- to_plot_dot %>% 
  filter(tau > 0.5, tau_prop > 0.7) %>% 
  distinct(gene) %>% pull()
genes_oi2 <- genes_oi %>% enframe() %>% filter(value %in% genes_oi2) %>% pull(value)

p1 <- to_plot_dot %>% 
  filter(tau > 0.5,
         tau_prop > 0.7) %>%
  mutate(region = region %>% as_factor() %>% fct_relevel(region_order),
         within_area_subclass = within_area_subclass %>% as_factor() %>% fct_relevel(subclass_order)) %>%
  ggplot() +
  geom_point(aes(x = region, y = gene %>% fct_relevel(genes_oi2) %>% fct_rev(), color = mean_expr, size = prop_count)) +
  scale_color_gradient(low = "white", high = "black") +
  labs(y = "Conserved subclass marker genes",
       x = "",
       color = "Mean Expression",
       size = "Proportion Expressing") +
  facet_wrap(~within_area_subclass, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        panel.background = element_rect(fill = "white", color = "grey"),
        panel.grid = element_blank())

p1
pdf(file = here("plots", "dotplot_conserved_subclass_glia_markers.pdf"),
    height = 15,
    width = 12,
    useDingbats = F)
print(p1)
dev.off()


#load astrocyte umap coords and plot genes of interest with A-P gradient
astro_genes <- c("PRRX1", "NR4A3", "LGR6", "PAX6", "SOX2", "NES", "S100B", "AQP4", "SLC1A3", "CRYAB", "HOPX")

seurat_obj@reductions <- reductions
cells_use <- to_plot %>% 
  filter(within_area_subclass == "Astro",
         UMAP_1 < 0,
         UMAP_2 < 5) %>%
  pull(sample_id)

seurat_obj <- ScaleData(seurat_obj, features = astro_genes)
p1 <- FeaturePlot(seurat_obj, reduction = "umap", features = astro_genes, cells = cells_use, raster = T, coord.fixed = 1)
p1
pdf(file = here("plots", "Astro_gradient_genes_of_interest.pdf"),
    height = 10,
    width = 10,
    useDingbats = F)
print(p1)
dev.off()

#astrocyte umap colored by putative subtypes
p1 <- seurat_obj@meta.data %>%
  filter(sample_id %in% cells_use) %>%
  mutate(astro_assignment = case_when(cross_area_cluster %in% c("Astro_2", "Astro_3") ~ "ILM",
                                      cross_area_cluster %in% c("Astro_4", "Astro_5") ~ "Fibrous",
                                      cross_area_cluster %in% c("Astro_1") ~ "Protoplasmic",
                                      TRUE ~ "none")) %>%
  filter(astro_assignment != "none") %>%
  ggplot() +
  geom_scattermore(aes(x = UMAP_1, y = UMAP_2, color = astro_assignment), pointsize = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  theme(aspect.ratio = 1)

pdf(file = here("plots", "Astro_assignment.pdf"),
    height = 5,
    width = 5,
    useDingbats = F)
print(p1)
dev.off()

#plot UMAP of consensus types on glial map
label_coords <- seurat_obj@meta.data %>% 
  group_by(cross_area_cluster) %>%
  summarise(mean_UMAP_1 = mean(UMAP_1),
            mean_UMAP_2 = mean(UMAP_2)) %>%
  ungroup()



p1 <-  ggplot() +
  geom_scattermore(data = seurat_obj@meta.data, aes(x = UMAP_1, y = UMAP_2, color = cross_area_cluster)) +
  geom_text_repel(data = label_coords, aes(x = mean_UMAP_1, y = mean_UMAP_2, label = cross_area_cluster), size = 2) +
  theme_void() +
  theme(aspect.ratio = 1) 

pdf(file = here("plots", "glia_umap_consensus_clusters.pdf"),
    height = 5,
    width = 5,
    useDingbats = F)
print(p1)
dev.off()  


#astro proportions plot
to_plot_prop <- seurat_obj@meta.data %>% 
  filter(cross_area_subclass == "Astro",
         dataset != "lein_10x_layer5_only") %>% 
  
  group_by(cross_area_cluster, region, donor) %>%
  count() %>%
  ungroup() %>%
  
  group_by(region, donor) %>%
  mutate(total_region = sum(n)) %>%
  ungroup() %>% 
  
  mutate(prop = n/total_region) %>%
  
  group_by(cross_area_cluster, region) %>%
  summarise(mean_prop = mean(prop),
            n_samp = n(),
            sd_prop = sd(prop)) %>%
  ungroup() %>%
  
  mutate(region = region %>% as_factor() %>% fct_relevel(region_order),
         se_prop = sd_prop / sqrt(n_samp)) 

p1 <- to_plot_prop %>% 
  
  ggplot() +
  geom_point(aes(x = region, y = mean_prop)) +
  geom_errorbar(aes(x = region, ymin = mean_prop - se_prop, ymax = mean_prop + se_prop)) +
  scale_y_log10() +
  labs(y = "Proportion of Astrocyte\n(mean ± se)",
       title = "Consensus astrocyte proportions across donors per region",
       x = "") +
  facet_wrap(~cross_area_cluster, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
p1  

pdf(file = here("plots", "astro_consensus_type_proportion_of_all_astro.pdf"),
    height = 5,
    width = 8,
    useDingbats = F)
print(p1)
dev.off()  
