library(tidyverse)
library(here)
library(scattermore)
library(Seurat)

meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))   
reductions <- readRDS(here("umap_coords", "Neighborhood_UMAPS", "MGE_reductions.RDS"))

to_plot <- reductions$umap@cell.embeddings %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(meta, by = "sample_id")

cluster_labs <- to_plot %>%
  group_by(cross_area_cluster) %>%
  summarise(mean_UMAP_1 = mean(UMAP_1),
            mean_UMAP_2 = mean(UMAP_2)) %>%
  ungroup()


 p1 <-  ggplot() +
  geom_scattermore(data = to_plot, aes(x = UMAP_1, y = UMAP_2, color = cross_area_cluster)) +
  geom_label(data = cluster_labs, aes(x = mean_UMAP_1, y = mean_UMAP_2, label = cross_area_cluster, fill = cross_area_cluster), size = 2) +
  theme_void() +
  theme(aspect.ratio = 1)

 pdf(file = here("plots", "UMAP_MGE_consensus_types.pdf"),
     height = 8,
     width = 10,
     useDingbats = F)
 print(p1)
 dev.off()
 
 
 #load gene expression and plot V1 specific SST genes
 #load gene expression for features plots
 file_to_load <- list.files(here("data"), full.names = T) %>%
   enframe() %>%
   filter(value %>% str_detect("mge"))
 
 genes_oi <- c(read.csv(here("markers", "V1_SST_subtype_markers.csv")) %>%
   as_tibble() %>%
   filter(p_val_adj < 0.05,
          avg_log2FC > 0.5) %>%
          #cluster %in% c("Sst_7", "Sst_8", "Sst)9")) %>%
   distinct(gene) %>% pull(), "SST")
 
 
 mat <- readRDS(file_to_load$value[1])
 mat <- mat[ , which(colnames(mat) %in% to_plot$sample_id)]
 mat <- CreateSeuratObject(mat)
 mat <- NormalizeData(mat)
 mat <- mat@assays$RNA@data[genes_oi, ]
 
 for(i in 2:nrow(file_to_load)){
   print(str_c("loading ", i))
   mat_tmp <- readRDS(file_to_load$value[i])
   mat_tmp <- mat_tmp[ , which(colnames(mat_tmp) %in% to_plot$sample_id)]
   mat_tmp <- CreateSeuratObject(mat_tmp)
   mat_tmp <- NormalizeData(mat_tmp)
   mat_tmp <- mat_tmp@assays$RNA@data[genes_oi, ]
   
   mat <- cbind(mat, mat_tmp)
   gc()
 }
 
 dim(mat)  
 gc()

 #make plot
 genes_oi 
 read.csv(here("markers", "V1_SST_subtype_markers.csv")) %>%
   as_tibble() %>%
   filter(p_val_adj < 0.05,
          avg_log2FC > 0.3) %>%
   filter(cluster %in% c("Sst_7", "Sst_8", "Sst_9") ) %>%
   count(gene) %>% filter(n == 3)
 
 gene_oi <- "LTBP1"
 
 to_plot2 <- mat[gene_oi, ] %>%
   enframe(name = "sample_id", value = "expr") %>%
   left_join(to_plot, by = "sample_id")
 
 
(p1 <-  ggplot() +
   geom_scattermore(data = to_plot2, aes(x = UMAP_1, y = UMAP_2, color = expr)) +
    scale_color_viridis_c(option = "B") +
    labs(title = gene_oi,
         color = "Log-Normalized Expression") +
   theme_void() +
   theme(aspect.ratio = 1))
  
 pdf(file = here("plots", str_c("MGE_featureplot_" , gene_oi, ".pdf")),
     height = 5,
     width = 5,
     useDingbats = F)
 print(p1)
 dev.off()
 