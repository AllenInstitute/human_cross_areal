library(tidyverse)
library(Seurat)
library(here)
library(scattermore)

#load reductions and metadata and combine to plot UMAPs
meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))        

reductions <- readRDS(here("umap_coords", "Neighborhood_UMAPS", "Glia_reductions.RDS"))

region_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG", "V1")

meta_sub <- meta %>%
  filter(within_area_subclass == "Astro") %>%
  as.data.frame()

rownames(meta_sub) <- meta_sub$sample_id

#load matrices
file_to_load <- list.files(here("data"), full.names = T) %>%
  enframe() %>%
  filter(value %>% str_detect("glia"))

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

#make seurat object
seurat_obj <- CreateSeuratObject(counts = mat, meta.data = meta_sub)
seurat_obj <- NormalizeData(seurat_obj)

#find markers between V1 prop and ILM subtypes
Idents(seurat_obj) <- seurat_obj$region
seurat_obj_sub <- subset(seurat_obj, idents = "V1")
Idents(seurat_obj_sub) <- seurat_obj_sub$within_area_cluster

proto_markers <- FindMarkers(seurat_obj_sub, ident.1 = c("Astro_1", "Astro_3"), ident.2 = c("Astro_2", "Astro_5", "Astro_4"), 
                             slot = "data", test.use = "wilcox", max.cells.per.ident = 200, only.pos = T)

ilm_markers <- FindMarkers(seurat_obj_sub, ident.1 = c("Astro_2", "Astro_5"), ident.2 = c("Astro_1", "Astro_3", "Astro_4"), 
                             slot = "data", test.use = "wilcox", max.cells.per.ident = 200, only.pos = T)

fibrous_markers <- FindMarkers(seurat_obj_sub, ident.1 = c("Astro_4"), ident.2 = c("Astro_1", "Astro_3", "Astro_2", "Astro_5"), 
                             slot = "data", test.use = "wilcox", max.cells.per.ident = 200, only.pos = T)


proto_subtype_markers <- FindMarkers(seurat_obj_sub, ident.1 = "Astro_1", ident.2 = "Astro_3", slot = "data", test.use = "wilcox",
            max.cells.per.ident = 200)

ilm_subtype_markers <- FindMarkers(seurat_obj_sub, ident.1 = "Astro_2", ident.2 = "Astro_5", slot = "data", test.use = "wilcox",
            max.cells.per.ident = 200)




major_markers <- bind_rows(
  proto_markers %>%
    mutate(type = "proto") %>%
    as_tibble(rownames = "gene"), 

  ilm_markers %>%
    mutate(type = "ilm") %>%
    as_tibble(rownames = "gene"),
  
  fibrous_markers %>%
    mutate(type = "fibrous") %>%
    as_tibble(rownames = "gene")) %>%
  
  filter(p_val_adj < 0.05,
         avg_log2FC > 0.5) 

seurat_obj_sub_ds <- subset(seurat_obj_sub, downsample = 200)
seurat_obj_sub_ds <- ScaleData(seurat_obj_sub_ds, features = major_markers$gene)

levels(seurat_obj_sub_ds) <- c("Astro_1", "Astro_3", "Astro_2", "Astro_5", "Astro_4")
DoHeatmap(seurat_obj_sub_ds, slot = "scale.data", features = major_markers$gene, raster = T)

#curate major markers for most binary
curated_markers <- major_markers %>%
  mutate(prop_diff = pct.1 - pct.2) %>%
  filter(prop_diff > 0.3)

curated_markers <- c("CABLES1", "NRP1", "KCTD8", "RIMS1", "DGKZ", #proto
                     
                     
                     "GRIA1", "LOC105370504", "LINC01411", "LOC105373158", #ilm
                     "FABP7", "SEZ6L", "LOC105373575", "FAM179A", "LMO2",
                    "AEBP1", "SH3GL2",
                    
                    "ADAMTSL3", "CD44", "CPAMD8", "AQP1", "ABI3BP", "SIDT1", #fibrous
                    "LOC105370803", "PLCB4", "ANGPTL4", "CCDC85A", "TTN",
                    "SMOC1", "GALNT13", "SORCS3", "LOC105378850", "PRR16",
                    "KCND2", "BNC2", "MOB3B", "HSPB1", "GRIP1", "TNR", "KSR2") 

#after running everything restrict gene list further
#add proto and ILM subtype markers from below
curated_markers<- c("PRRX1", "NR4A3", "LGR6", "PAX6", "SOX2", "NES", "S100B", "AQP4", "ALDH1L1", "SLC1A3", "CRYAB", "HOPX")

curated_markers <- c("GFAP", "GLUL", "SLC1A2", "SLC1A3", "AQP4",
                    
                      "CABLES1", "KCTD8", "RIMS1",  #proto
                     
                     "COX1", "COX2", "COX3", #Astro_1 (proto subtype- cell state)
                     
                     "LOC105370504", "LOC105373158", #ilm
                     "FABP7", "LOC105373575",
                     
                     "FAM19A1", "OAF", "CFAP47", "CFAP54", "HYDIN",
                     "SPARC", "GMPR", "MAN1A1", "NLGN1", "GUCY1A3", #ILM subtype markers
                    
                    "ADAMTSL3", "CD44",  "AQP1",  #fibrous
                     "TTN", "KCND2") 

curated_markers <- unique(c("GFAP", "PRRX1", "S100B", "CRYAB", "GLUL", "SLC1A2", "SLC1A3", "ALDH1L1", "AQP4",
                    
                     "NR4A3", "NES", 
                      "CABLES1", "KCTD8", "RIMS1",  #proto
                     
                     "COX1", "COX2", "COX3", #Astro_1 (proto subtype- cell state)
                     "CRYAB",
                     
                     "LOC105370504", "LOC105373158", #ilm
                     "FABP7", "LOC105373575",
                     
                     "FAM19A1", "OAF", "CFAP47", "CFAP54", "HYDIN",
                     "SPARC", "GMPR", "MAN1A1", "NLGN1", "GUCY1A3", #ILM subtype markers
                    "LGR6", "HOPX",
                     
                    "ADAMTSL3", "CD44",  "AQP1",  #fibrous
                     "TTN", "KCND2") )

seurat_obj_sub_ds <- ScaleData(seurat_obj_sub_ds, features = curated_markers)
DoHeatmap(seurat_obj_sub_ds, slot = "scale.data", features = curated_markers, raster = T)

#

to_plot <- seurat_obj_sub@meta.data %>%
  as_tibble() %>%
  left_join(
    seurat_obj_sub@assays$RNA@data[curated_markers, ] %>%
      as_tibble(rownames = "gene") %>%
      gather(key = "sample_id", value = "expr", -gene),
    by = "sample_id"
  )




p1 <- to_plot %>%
  
  mutate(within_area_cluster = within_area_cluster %>% as_factor() %>% fct_relevel(c("Astro_1", "Astro_3", "Astro_2", "Astro_5", "Astro_4")),
         gene = gene %>% as_factor() %>% fct_relevel(curated_markers)) %>%
  ggplot() +
  geom_violin(aes(x = within_area_cluster, y = expr, fill = within_area_cluster), scale = "width") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "",
       y= "Normalized expression",
       fill = "V1 astrocyte types") +
  facet_wrap(~gene, ncol = 1, scales = "free_y", strip.position = "right") +
  theme(aspect.ratio = 1/2,
        strip.text.y =  element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
p1
pdf(here("plots", "V1_Astro_markers_violin.pdf"),
    useDingbats = FALSE,
    height = 20,
    width = 6)
print(p1)
dev.off()

#plot horizontally
p1 <- to_plot %>%
  
  mutate(within_area_cluster = within_area_cluster %>% as_factor() %>% fct_relevel(c("Astro_1", "Astro_3", "Astro_2", "Astro_5", "Astro_4")),
         gene = gene %>% as_factor() %>% fct_relevel(curated_markers)) %>%
  ggplot() +
  geom_violin(aes(y = within_area_cluster, x = expr, fill = within_area_cluster), scale = "width") +
  scale_fill_brewer(palette = "Set3") +
  labs(y = "",
       x= "Normalized expression",
       fill = "V1 astrocyte types") +
  facet_wrap(~gene, nrow = 1, scales = "free_x") +
  theme(aspect.ratio = 2,
        strip.text.x =  element_text(angle = 270),
        axis.text.x = element_text(size = 4))
p1
pdf(here("plots", "V1_Astro_markers_violin_horizontal_plot2.pdf"),
    useDingbats = FALSE,
    height = 6,
    width = 20)
print(p1)
dev.off()

#now find proto and ilm subtype markers and add them above
proto_markers_sub <- proto_subtype_markers %>%
  as_tibble(rownames = "gene") %>%
  mutate(prop_thresh = abs(pct.1 - pct.2)) %>%
  filter(abs(avg_log2FC) > 0.2,
         p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))

ilm_markers_sub <- ilm_subtype_markers %>%
  as_tibble(rownames = "gene") %>%
  mutate(prop_thresh = abs(pct.1 - pct.2)) %>%
  filter(abs(avg_log2FC) > 0.5,
         p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))



to_plot <- seurat_obj_sub@meta.data %>%
  as_tibble() %>%
  left_join(
    seurat_obj_sub@assays$RNA@data[ilm_markers_sub$gene, ] %>%
      as_tibble(rownames = "gene") %>%
      gather(key = "sample_id", value = "expr", -gene),
    by = "sample_id"
  )




to_plot %>%
  
  mutate(within_area_cluster = within_area_cluster %>% as_factor() %>% fct_relevel(c("Astro_1", "Astro_3", "Astro_2", "Astro_5", "Astro_4")),
         gene = gene %>% as_factor() %>% fct_relevel(ilm_markers_sub$gene)) %>%
  
  filter(gene %in% ilm_markers_sub$gene[121:139]) %>%
  
  ggplot() +
  geom_violin(aes(x = within_area_cluster, y = expr, fill = within_area_cluster), scale = "width") +
  
  labs(x = "",
       y= "Normalized expression",
       fill = "V1 astrocyte types") +
  facet_wrap(~gene, ncol = 1, scales = "free_y", strip.position = "right") +
  theme(aspect.ratio = 1/2,
        strip.text.y =  element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
