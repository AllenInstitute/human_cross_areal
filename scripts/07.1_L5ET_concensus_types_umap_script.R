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

p1 <- to_plot %>%
 
  filter(UMAP_1 < 0) %>%
  
  ggplot() +
  geom_scattermore(aes(x = UMAP_1, y = UMAP_2, color = cross_area_cluster), pointsize = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  theme(aspect.ratio = 1)

pdf(file = here("plots", "L5_ET_consensus_types_umap.pdf"),
    height = 5,
    width = 5,
    useDingbats = F)
print(p1)
dev.off()
