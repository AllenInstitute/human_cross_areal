library(tidyverse)
library(here)

meta_ss <- readRDS(here("data", "meta_ss_3_29_22.RDS"))
area_order <-  c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "V1")

subclasses <- meta_ss %>% distinct(within_area_subclass) %>% pull()

layer_to_include <- meta_ss %>% 
  distinct(region, layer) %>% 
  mutate(layer_prop = 0,
         region_layer = str_c(region, "_", layer))

for(i in 1:length(subclasses)){
  
subclass_oi <- subclasses[i]  

to_plot <- meta_ss %>%
  filter(within_area_subclass == subclass_oi) %>%
  
  group_by(layer, region) %>%
  dplyr::count() %>%
  ungroup() %>%
  
  group_by(region) %>%
  mutate(total_n = sum(n)) %>%
  ungroup() %>%
  
  mutate(layer_prop = n / total_n) %>%
  dplyr::select(layer, region, layer_prop) %>%
  mutate(region_layer = str_c(region, "_", layer))

to_plot <- bind_rows(to_plot,
          layer_to_include %>%
            filter(!(region_layer %in% to_plot$region_layer)))

p1 <- to_plot %>%
  mutate(region = region %>% as_factor() %>% fct_relevel(area_order)) %>%
  ggplot() +
  geom_tile(aes(x = region,  y = layer %>% fct_rev(), fill = layer_prop)) +
  geom_text(aes(x = region, y = layer %>% fct_rev(), label = layer), color = "darkorange", size = 6) +
  labs(title = subclass_oi,
       y = "Layer",
       fill = "Layer prop.") +
  facet_wrap(~region, nrow = 1, scales = "free") +
  scale_fill_gradient(low = "white", high = "black") +
  theme_void() +
  theme(aspect.ratio = 3,
        text = element_text(face = "bold"))

pdf(file = here("plots", str_c(subclass_oi %>% str_remove_all("/") %>%
                                 str_remove_all("_") %>%
                                 str_remove_all("-") %>%
                                 str_remove_all(" "), "_layer_prop_heatmap.pdf")),
    useDingbats = F,
    height = 4,
    width = 6)
print(p1)
dev.off()
}
