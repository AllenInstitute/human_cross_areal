library(tidyverse)
library(here)

meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))

#proportion of class across all subclasses
subclass_props <- meta %>% 
  filter(dataset != "lein_10x_layer5_only") %>%
  
  dplyr::select(region, within_area_subclass, class) %>%
  
  group_by(region, class) %>%
  mutate(class_count = n()) %>%
  ungroup() %>%
  
  group_by(region, within_area_subclass) %>%
  mutate(subclass_count = n()) %>%
  ungroup() %>% 
  
  distinct() %>%
  
  mutate(subclass_prop = subclass_count / class_count) %>%

  arrange(region, within_area_subclass)  

#proportion of neurons across all neuronal subclasses

subclass_props <- meta %>% 
  filter(dataset != "lein_10x_layer5_only",
         class != "non-neuronal") %>%
  
  dplyr::select(region, within_area_subclass) %>%
  
  group_by(region) %>%
  mutate(region_count = n()) %>%
  ungroup() %>%
  
  group_by(region, within_area_subclass) %>%
  mutate(subclass_count = n()) %>%
  ungroup() %>% 
  
  distinct() %>%
  
  mutate(subclass_prop = subclass_count / region_count) %>%

  arrange(region, within_area_subclass)  


#find spearman correlations across subclasses
subclasses <- subclass_props %>% distinct(within_area_subclass) %>% pull()
cor_mat <- matrix(0, nrow = length(subclasses), ncol = length(subclasses))
rownames(cor_mat) <- colnames(cor_mat) <- subclasses


for(i in 1:length(subclasses)){
  for(j in 1:length(subclasses)){
    cor_mat[j, i] <- cor(subclass_props %>% filter(within_area_subclass == colnames(cor_mat)[i]) %>% pull(subclass_prop),
    subclass_props %>% filter(within_area_subclass == rownames(cor_mat)[j]) %>% pull(subclass_prop), method = "spearman")
    
  }
}

pdf(file = here("plots", "subclass_spearman_cor_heatmap.pdf"),
    useDingbats = FALSE,
    height = 5,
    width = 5)
pheatmap::pheatmap(cor_mat)
dev.off()

pdf(file = here("plots", "subclass_spearman_cor_heatmap_neurons.pdf"),
    useDingbats = FALSE,
    height = 5,
    width = 5)
pheatmap::pheatmap(cor_mat)
dev.off()

