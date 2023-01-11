library(tidyverse)
library(Seurat)
library(here)
library(scattermore)

#load reductions and metadata and combine to plot UMAPs
meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))        

reductions <- readRDS(here("umap_coords", "Neighborhood_UMAPS", "Deep_Exc_reductions.RDS"))

region_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG", "V1")

meta_sub <- meta %>%
  filter(within_area_subclass == "L6 CT") %>%
  as.data.frame()

rownames(meta_sub) <- meta_sub$sample_id

#load matrices
file_to_load <- list.files(here("data"), full.names = T) %>%
  enframe() %>%
  filter(value %>% str_detect("deep_exc"))

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

Idents(seurat_obj) <- seurat_obj$region
levels(seurat_obj) <- region_order
region_markers <- FindAllMarkers(seurat_obj, slot = "data", test.use = "wilcox", only.pos = T, max.cells.per.ident = 500)
write.csv(region_markers, here("markers", "L6CT_regional_markers.csv"))

genes_oi <- region_markers %>%
  as_tibble() %>%
  mutate(prop_diff = pct.1 - pct.2) %>%
  filter(avg_log2FC > 0.5,
         p_val_adj < 0.05) %>%
  mutate(cluster = cluster %>% as_factor())
  

genes_oi %>%
  count(cluster)

#plot genes_oi
seurat_obj_sub <- subset(seurat_obj, downsample = 200)
seurat_obj_sub <- ScaleData(seurat_obj_sub, features = genes_oi$gene)
DoHeatmap(seurat_obj_sub, features = genes_oi$gene, raster = T)

#plot heatmap plot 
umap_genes <- genes_oi %>%
  filter(cluster == "V1",
         gene != "XIST") %>% 
  distinct(gene) %>%
  arrange(gene) 


p1 <- DoHeatmap(seurat_obj_sub, features = umap_genes$gene, raster = T) +
  scale_fill_viridis_c(option = "B") +
  theme(aspect.ratio = 2,
        axis.text.y = element_blank())

pdf(here("plots", "L6CT_V1_marker_heatmap.pdf"),
    useDingbats = FALSE,
    height = 6,
    width = 6)
print(p1)
dev.off()


#plot violin
violin_genes <- genes_oi %>%
  filter(cluster == "V1",
         gene != "XIST",
         prop_diff > 0.3) %>% 
  distinct(gene) %>%
  arrange(gene) 

to_plot <- reductions$umap@cell.embeddings %>%
  as_tibble(rownames = "sample_id") %>%
  filter(sample_id %in% meta_sub$sample_id) %>%
  left_join(meta_sub, by = "sample_id") %>%
  left_join(
    seurat_obj@assays$RNA@data[violin_genes$gene, ] %>%
      as_tibble(rownames = "gene") %>%
      gather(key = "sample_id", value = "expr", -gene),
    by = "sample_id"
  )

colors_use <- to_plot %>% distinct(region, region_color) %>% deframe()

p1 <- to_plot %>%

  mutate(region = region %>% as_factor() %>% fct_relevel(region_order),
         gene = gene %>% as_factor() %>% fct_relevel(violin_genes$gene)) %>%
  ggplot() +
  geom_violin(aes(x = region, y = expr, fill = region), scale = "width") +
  scale_fill_manual(values = colors_use) +
  labs(x = "",
       y= "Normalized expression",
       fill = "Region") +
  facet_wrap(~gene, ncol = 2, scales = "free_y", strip.position = "right") +
  theme(aspect.ratio = 1/2,
        strip.text.y =  element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))
p1
pdf(here("plots", "L6CT_V1_markers_violin.pdf"),
    useDingbats = FALSE,
    height = 12,
    width = 12)
print(p1)
dev.off()

#run GO analysis
library(KEGGREST)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)


gene_universe <- names(which(rowSums(seurat_obj@assays$RNA@counts) > 1))

gene_universe <- mapIds(org.Hs.eg.db, gene_universe, "ENTREZID", "SYMBOL") %>%
  enframe() %>%
  set_names("SYMBOL", "ENTREZID") %>%
  filter(!(is.na(ENTREZID))) 

to_test <- umap_genes %>%
  distinct(gene) %>%
  filter(gene %in% gene_universe$SYMBOL) %>%
  left_join(gene_universe, by = c("gene" = "SYMBOL"))

go_results_bp <- enrichGO(
  gene = to_test$ENTREZID,
  universe = gene_universe$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  maxGSSize = 700,
  readable = T
)

go_results_cc <- enrichGO(
  gene = to_test$ENTREZID,
  universe = gene_universe$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  maxGSSize = 700,
  readable = T
)

go_results_mf <- enrichGO(
  gene = to_test$ENTREZID,
  universe = gene_universe$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  maxGSSize = 700,
  readable = T
)

#plot results
to_plot <- bind_rows(go_results_bp@result %>%
            as_tibble() %>%
            filter(p.adjust < 0.01) %>%
            
            arrange(p.adjust, desc(Count)) %>%
            dplyr::slice(1:10) %>%
            mutate(Description = Description %>% as_factor() %>% fct_rev(),
                   cat = "GO:BP"), 
            
          go_results_cc@result %>%
            as_tibble() %>%
            filter(p.adjust < 0.01) %>%
            
            arrange(p.adjust, desc(Count)) %>%
            dplyr::slice(1:10) %>%
            mutate(Description = Description %>% as_factor() %>% fct_rev(),
                   cat = "GO:CC"), 
          
          go_results_mf@result %>%
            as_tibble() %>%
            filter(p.adjust < 0.01) %>%
            
            arrange(p.adjust, desc(Count)) %>%
            dplyr::slice(1:10) %>%
            mutate(Description = Description %>% as_factor() %>% fct_rev(),
                   cat = "GO:MF"))
  
  p1 <- to_plot %>%
    group_by(cat) %>%
    arrange(cat, desc(Count)) %>%
    mutate(
      Description = as.character(Description),
      Description = Description %>% as_factor() %>% fct_rev()) %>%
    ungroup() %>%
    
    ggplot() +
    geom_col(aes(x = Count, y = Description, fill = cat)) +
    scale_fill_brewer(palette = "Set2") +
    facet_wrap(~cat, ncol = 1, scales = "free", strip.position = "right") +
    labs(fill = "",
         y = "",
         x = "Gene count (padj < 0.01)") +
    theme(aspect.ratio = 2,
          strip.text = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major.x = element_line(color = "grey"))
p1
  pdf(here("plots", "L6CT_GO_regional_genes.pdf"),
      useDingbats = FALSE,
      height = 10,
      width = 8)
print(p1)
dev.off()


