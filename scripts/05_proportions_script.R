library(tidyverse)

meta_10x <- readRDS(here("data", "meta_10x_2_11_22.RDS"))
meta_ss <- readRDS(here("data", "meta_ss_3_29_22.RDS"))

region_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG","V1")
subclass_order <- c("L2/3 IT", "L4 IT", "L5 IT", "L6 IT", "L6 IT Car3",
                    "L5 ET", "L5/6 NP", "L6b", "L6 CT",
                    "Lamp5 Lhx6", "Lamp5", "Sncg", "Vip", "Pax6",
                    "Chandelier", "Pvalb", "Sst", "Sst Chodl",
                    "OPC", "Oligo", "Astro", "Micro/PVM", "Endo", "VLMC")


#make dataset proportion counts broken out by donor
meta_all <- bind_rows(meta_10x, meta_ss)

p1 <- meta_all %>%
  mutate(cleaned_dataset = case_when(dataset %>% str_detect("10x_all_layers") ~ "10x all layers",
                                     dataset %>% str_detect("lein_ss") ~ "SSv4",
                                     dataset %>% str_detect("layer5") ~ "10x layer 5 only")) %>%
  group_by(region, cleaned_dataset) %>%
  count() %>%
  ungroup() %>%
  
  mutate(cleaned_dataset = cleaned_dataset %>% as_factor() %>% fct_relevel(c("10x all layers", "10x layer 5 only", "SSv4")),
         region = region %>% as_factor() %>% fct_relevel(region_order)) %>%
  
  ggplot() +
  geom_col(aes(x = cleaned_dataset, y = n, fill = cleaned_dataset)) +
  scale_fill_manual(values = c(`10x all layers` = "#df7258",
                               `10x layer 5 only` = "#dd9a80",
                               SSv4 = "#78a4c9")) +
  scale_y_continuous(breaks = c(0, 50000, 100000),
                     labels = c("0", "50000", "100000")) +
  labs(y = "Number of nuclei",
       x = "",
       fill = "Dataset") +
  facet_wrap(~region, nrow = 1) +
  theme(aspect.ratio = 2,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        strip.background = element_rect(fill = "#2c3e50"),
        strip.text = element_text(color = "white"))
p1
pdf(file = here("plots", "dataset_nuclei_counts.pdf"),
    height = 4,
    width = 6,
    useDingbats = F)
print(p1)
dev.off()


p1 <- meta_all %>%
  mutate(cleaned_dataset = case_when(dataset %>% str_detect("10x_all_layers") ~ "10x all layers",
                                     dataset %>% str_detect("lein_ss") ~ "SSv4",
                                     dataset %>% str_detect("layer5") ~ "10x layer 5 only")) %>%
  distinct(region, cleaned_dataset, sex, donor) %>%
  group_by(region, cleaned_dataset, sex) %>%
  count() %>%
  ungroup() %>%
  
  mutate(cleaned_dataset = cleaned_dataset %>% as_factor() %>% fct_relevel(c("10x all layers", "10x layer 5 only", "SSv4")),
         region = region %>% as_factor() %>% fct_relevel(region_order))  %>%
  
  ggplot() +
  geom_col(aes(x = cleaned_dataset, y = n, fill = sex)) +
  scale_fill_manual(values = c(Female = "#bd3a8d",
                               Male = "#2b4086")) +
  labs(y = "Number of donors",
       x = "",
       fill = "") +
  facet_wrap(~region, nrow = 1) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
        strip.background = element_rect(fill = "#2c3e50"),
        strip.text = element_text(color = "white"))
p1
pdf(file = here("plots", "dataset_donor_counts.pdf"),
    height = 4,
    width = 6,
    useDingbats = F)
print(p1)
dev.off()


#stacked barplot 
colors_use <- meta_10x %>% distinct(within_area_subclass, within_area_subclass_color) %>% deframe()

  #10x
  p1 <- meta_10x  %>%
  
    filter(dataset != "lein_10x_layer5_only") %>%
    mutate(region = region %>% as_factor() %>% fct_relevel(region_order),
           within_area_subclass = within_area_subclass %>% as_factor() %>% fct_relevel(subclass_order)) %>%
    
    ggplot() +
    geom_bar(aes(x = region, fill = within_area_subclass), position = "fill") +
    scale_fill_manual(values = colors_use, breaks = subclass_order) +
    
    labs(title = "10x nuclei excluding L5 dissections") +
    
    facet_wrap(~class, ncol = 1) +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2, size = 6)
          )
  p1
  pdf(file = here("plots", "stacked_barplot_proportions_by_region_10x.pdf"),
      height = 8,
      width = 6,
      useDingbats = F)
  print(p1)
  dev.off()
  
  #ss
  meta_ss  %>%
  
    filter(dataset != "lein_10x_layer5_only") %>%
    mutate(region = region %>% as_factor() %>% fct_relevel(region_order),
           within_area_subclass = within_area_subclass %>% as_factor() %>% fct_relevel(subclass_order)) %>%
    
    ggplot() +
    geom_bar(aes(x = region, fill = within_area_subclass), position = "fill") +
    scale_fill_manual(values = colors_use, breaks = subclass_order) +
    
    labs(title = "SMARTseq nuclei") +
    
    facet_wrap(~class, ncol = 1) +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2, size = 6)
          )
  
  
#plot proportions of class
neighborhoods <- meta_10x %>% distinct(neighborhood) %>% pull()
neighborhood_use <- neighborhoods[1]    
  
  to_plot <- meta_10x %>% 
  
    filter(dataset != "lein_10x_layer5_only") %>%
    
    group_by(region, donor, class) %>%
    mutate(region_donor_class_n = n()) %>%
    ungroup() %>%
    
    group_by(region, donor, class, within_area_subclass) %>%
    mutate(region_donor_class_subclass_n = n()) %>%
    ungroup() %>%
    
    mutate(prop_class = region_donor_class_subclass_n / region_donor_class_n * 100) %>%
    
    filter(neighborhood == neighborhood_use) %>%
    
    group_by(region, class, within_area_subclass) %>%
    summarise(mean_prop = mean(prop_class),
              sd_prop = sd(prop_class)) %>%
    ungroup() 
  
  subclass_order_tmp <- subclass_order[subclass_order %in% to_plot$within_area_subclass]
  region_order_tmp <- region_order[region_order %in% to_plot$region]
  
  to_plot <- to_plot %>%
    mutate(region = region %>% as_factor() %>% fct_relevel(region_order_tmp),
           within_area_subclass = within_area_subclass %>% as_factor() %>% fct_relevel(subclass_order_tmp)) %>%
    
    arrange(within_area_subclass, region) %>%
    
    mutate(region_subclass = str_c(region, "_", within_area_subclass) %>% as_factor())
  
  colors_use <- meta_10x %>% distinct(within_area_subclass, within_area_subclass_color) %>% deframe()
  
p1 <- to_plot %>%
    ggplot() +
    geom_point(aes(x = region, y = mean_prop, color = within_area_subclass), size = 3) +
    geom_errorbar(aes(x = region, ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop)) +
    
    facet_wrap(~within_area_subclass, nrow = 1) +
    
    scale_color_manual(values = colors_use) +
    scale_y_log10() +
    
    labs(y = "Percent of class",
         x = "") +  
    
    theme(aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2),
          legend.position = "none",
          panel.background = element_rect(color = "black", fill = "white"),
          panel.grid.major.y = element_line(color = "lightgrey", size = 1, linetype = 2),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
p1

  pdf(file = here("plots", str_c("subclass_proportions_by_region_10x_", neighborhood_use %>% str_remove_all(" "), ".pdf")),
      height = 4,
      width = 6,
      useDingbats = F)
  print(p1)
  dev.off() 
  

  
#make E/I ratio plots by layer
meta_ss %>%
  group_by(region, layer, class) %>%
  mutate(class_n = n()) %>%
  ungroup() %>%
  
  distinct(region, layer, class, class_n) %>%
  spread(key = class,
         value = class_n) %>%
  
  mutate(ei_ratio = excitatory / inhibitory) %>%
  
  
  ggplot() +
  geom_tile(aes(x = 1, y = layer %>% as.factor() %>% fct_rev(), fill = ei_ratio)) +
  geom_text(aes(x = 1, y = layer %>% as.factor() %>% fct_rev(), 
                label = round(ei_ratio, 1)), size = 3) +
  
  scale_fill_gradientn(colours = rainbow(8), breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16)) +
  
  labs(y = "") +
  
  facet_wrap(~region, nrow = 1, scales = "free") +
  
  
  theme(aspect.ratio = 3,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())





meta_ss %>%
  group_by(region, layer, class, donor) %>%
  mutate(class_n = n()) %>%
  ungroup() %>%
  
  distinct(region, layer, class, class_n, donor) %>%
  spread(key = class,
         value = class_n) %>%
  
  mutate(ei_ratio = excitatory / inhibitory) %>%
  filter(excitatory > 3,
         inhibitory > 3) %>%
  
  
  ggplot() +
  geom_tile(aes(x = donor, y = layer %>% as.factor() %>% fct_rev(), fill = ei_ratio)) +
  geom_text(aes(x = donor, y = layer %>% as.factor() %>% fct_rev(), 
                label = round(ei_ratio, 1)), size = 3) +
  
  scale_fill_gradientn(colours = rainbow(8), breaks = c(0, 2, 4, 6, 8, 10, 12, 14, 16)) +
  
  labs(y = "") +
  
  facet_wrap(~region, nrow = 1, scales = "free") +
  
  
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))


#EI ratio barplots
colors_use <- meta_10x %>% distinct(region, region_color) %>% deframe()

p1 <- meta_10x %>%
  filter(!(dataset %>% str_detect("layer5"))) %>%
  filter(!(class %>% str_detect("non-neuronal"))) %>%
  
  group_by(region, donor, class) %>%
  count() %>%
  ungroup() %>%
  
  spread(key = class, value = n) %>%
  mutate(ei_ratio = excitatory / inhibitory) %>%
  
  group_by(region) %>%
  summarise(mean_ei = mean(ei_ratio),
            sd_ei = sd(ei_ratio)) %>%
  ungroup() %>%
  
  mutate(region = region %>% as_factor() %>% fct_relevel(region_order)) %>%
  
  ggplot() +
  geom_col(aes(x = region, y = mean_ei, fill = region)) +
  geom_errorbar(aes(x = region, ymin = mean_ei, ymax = mean_ei + sd_ei)) +
  scale_fill_manual(values = colors_use) +
  labs(y = "E/I ratio",
       x = "",
       fill = "Region") +
  theme(aspect.ratio = 2,
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))

pdf(file = here("plots", "ei_ratio_10x.pdf"),
    height = 5,
    width = 4,
    useDingbats = F)
print(p1)
dev.off() 

#run ANOVA on proportion stats
to_plot <- meta_10x %>% 
  
  filter(dataset != "lein_10x_layer5_only") %>%
  
  group_by(region, donor, class) %>%
  mutate(region_donor_class_n = n()) %>%
  ungroup() %>%
  
  group_by(region, donor, class, within_area_subclass) %>%
  mutate(region_donor_class_subclass_n = n()) %>%
  ungroup() %>%
  
  mutate(prop_class = region_donor_class_subclass_n / region_donor_class_n * 100) %>%
  
  distinct(region, donor, within_area_subclass, prop_class)
  
subclasses <- to_plot %>% distinct(within_area_subclass) 
subclasses$Pval <- 999
subclasses$Pval_adj <- 999

for(i in 1:nrow(subclasses)){
  (test_oi <- to_plot %>% filter(within_area_subclass == subclasses$within_area_subclass[i]))
  
  model <- aov(data = test_oi, prop_class ~ factor(region))
  model_res <- summary(model)
  subclasses$Pval[i] <- model_res[[1]]$`Pr(>F)`[1]  
  subclasses$Pval_adj[i] <- subclasses$Pval[i] * nrow(subclasses)
}

results <- subclasses %>% mutate(sig = case_when(Pval_adj >= 0.05 ~ "ns",
                                      Pval_adj < 0.001 ~ "***",
                                      Pval_adj < 0.01 ~ "**",
                                      Pval_adj < 0.05 ~ "*",
                                      TRUE ~ "other"
                                      )) %>%
  mutate(Pval = Pval %>% scales::comma(),
         Pval_adj = Pval_adj %>% scales::comma())
View(results)
write.csv(results, here("tables", "proportion_ANOVA_results.csv"))


