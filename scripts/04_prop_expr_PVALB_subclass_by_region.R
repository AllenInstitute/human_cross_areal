library(tidyverse)
library(here)

meta <- readRDS(here("data", "meta_10x_2_11_22.RDS"))

files_to_load <- list.files(here("data"), full.names = T) %>%
  enframe() %>%
  filter(!(value %>% str_detect("meta")))

#load all nuclei and their expression of PVALB
data_list <- list()
for(i in 1:nrow(files_to_load)){
  print(i)
  mat <- readRDS(files_to_load$value[i])
  data_list[[i]] <- mat["PVALB", ] %>%
    enframe()
  gc()
}

data <- data.table::rbindlist(data_list)
to_plot <- data %>%
  as_tibble() %>%
  set_names("sample_id", "expr") %>%
  filter(sample_id %in% meta$sample_id) %>%
  left_join(meta, by = "sample_id")


#plot expression of PVALB across regions and subclasses
area_order <- c("ACC", "DLPFC", "M1", "S1", "ANG", "A1", "MTG", "V1")


to_plot_use <- to_plot %>% 
  filter(dataset != "lein_10x_layer5_only") %>%
  select(within_area_subclass, region, expr) %>%
  mutate(gene_count = case_when(expr > 0 ~ 1,
                                TRUE ~ 0)) %>%
  group_by(within_area_subclass, region) %>%
  mutate(total = n(),
         thresh_count = sum(gene_count)) %>%
  ungroup() %>%
  
  mutate(prop_expr = thresh_count / total) %>%
  
  distinct(within_area_subclass, region, prop_expr) %>%
  
  mutate(region = region %>% as_factor() %>% fct_relevel(area_order)) 
  
p1 <- to_plot_use %>%
  ggplot() +
  geom_col(aes(x = within_area_subclass, y = prop_expr)) +
  facet_wrap(~region, ncol = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))

p2 <- to_plot_use %>%
  ggplot() +
  geom_tile(aes(x = within_area_subclass, y = region, fill = prop_expr)) +
  scale_fill_viridis_c(option = "B") +
  theme(aspect.ratio = 1/3,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))

pdf(here("plots", "PVALB_prop_expr_barplot.pdf"),
    useDingbats = F,
    height = 6,
    width = 4)
print(p1)
dev.off()

pdf(here("plots", "PVALB_prop_expr_heatmap.pdf"),
    useDingbats = F,
    height = 5,
    width = 5)
print(p2)
dev.off()
