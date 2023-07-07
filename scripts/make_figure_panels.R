library(tidyverse)
library(broom)
library(pheatmap)
library(GGally)
library(ggtern)
library(janitor)
library(RColorBrewer)




theme_set(theme_bw())

dat_dir <- "data/"
output_dir <- "output/"
if (! dir.exists(output_dir)) dir.create(output_dir)


#### Figure 1 ####

# Figure 1E. Entropy
dat <- read_csv(file = "data/ca_entropy_analysis_100perm.csv") %>% 
  mutate(region = factor(region, levels = c("ACC", "DLPFC", "M1", "S1",
                                            "MTG", "A1", "ANG", "V1")))
cl <- read_csv(file = "data/region_subclass_cl_cnt.csv")

cl2 <- cl %>% 
  mutate(subclass = str_replace_all(within_area_subclass, "[ /-]", "")) %>% 
  select(c(region, subclass, cl_cnt:cl_norm))


dat2 <- dat %>%
  left_join(cl2, by = c("region", "subclass")) %>% 
  group_by(region, class, subclass, cl_cnt, nuc_cnt) %>% 
  summarize(ent_med = median(entropy)) %>% 
  mutate(region = factor(region, levels = c("ACC", "DLPFC", "M1", "S1",
                                            "MTG", "A1", "ANG", "V1")))

aov1 <- aov(ent_med ~ class, dat2) 
aov2 <- aov(ent_med ~ class + region, dat2) 
anova(aov1, aov2)
TukeyHSD(aov2, which = "class")

p1 <- dat2 %>% 
  ggplot(aes(x = region, y = ent_med, fill = class)) + 
  geom_boxplot() +
  facet_wrap(~ class, ncol = 1, scale = "free_y") +
  labs(x = "", y = "Entropy") +
  guides(x = guide_axis(angle = 90)) +
  theme(strip.text.x = element_blank(), 
        panel.grid = element_blank())
plot(p1)

ggsave(p1, filename = str_c(output_dir, "/region_entropy.pdf"), width = 2.75, height = 3)



# Figure 1F-H. Cluster counts

meta_10x <- readRDS("data/meta_10x_2_11_22.RDS")

meta_10x_rev <- meta_10x %>%
  mutate(
    class = factor(class),
    neighborhood = factor(neighborhood, levels = c("IT types", "Deep Exc", 
                                                   "MGE Inh", "CGE Inh", "non-neuronal")), 
    within_area_subclass = factor(within_area_subclass),
    within_area_subclass = fct_reorder(within_area_subclass,
                                       as.numeric(as.factor(neighborhood))),
    region = factor(region, levels = c("ACC", "DLPFC", "M1", "S1",
                                       "MTG", "A1", "ANG", "V1"))
  )


# Cluster counts
cl_cnt <- meta_10x_rev %>%
  group_by(region, class, neighborhood, within_area_subclass) %>%
  summarize(cl_cnt = n_distinct(within_area_cluster),
            nuc_cnt = n()) %>%
  group_by(within_area_subclass) %>%
  mutate(cl_norm = cl_cnt / max(cl_cnt))

write_csv(cl_cnt, file = str_c(output_dir, "/region_subclass_cl_cnt.csv"))


for (class1 in c("excitatory", "inhibitory", "non-neuronal")) {
  p1 <- meta_10x_rev %>%
    filter(class == class1) %>%
    group_by(region, class, neighborhood, within_area_subclass) %>%
    summarize(cl_cnt = n_distinct(within_area_cluster),
              nuc_cnt = n()) %>%
    group_by(within_area_subclass) %>%
    mutate(cl_norm = cl_cnt / max(cl_cnt)) %>% 
    ggplot(aes(x = region,
               y = fct_rev(within_area_subclass),
               size = cl_cnt,
               color = nuc_cnt)) +
    geom_point() +
    scale_size(range = c(1, 8)) +
    scale_color_continuous(trans = "log10") +
    labs(x = "", y = "") +
    guides(x = guide_axis(angle = -90))
  
  pdf(
    file = str_c(output_dir, "/region_subclass_cl_cnt_", class1, ".pdf"),
    width = 4,
    height = 3.5
  )
  plot(p1)
  dev.off()
}



#### Figure 2, S12 ####
meta <- readRDS("data/meta_10x_2_11_22.RDS")

prop <- meta %>%
  filter(layer == "All") %>%
  group_by(region, class) %>%
  mutate(class_count = n()) %>%
  group_by(class, within_area_subclass, region) %>%
  summarize(subclass_count = n(),
            subclass_prop = subclass_count / class_count) %>%
  unique() %>% 
  ungroup()

prop_cor <- prop %>%
  select(-subclass_count) %>% 
  pivot_wider(names_from = "within_area_subclass", 
              values_from = "subclass_prop") %>% 
  select(-region) %>% 
  group_nest(class) %>% 
  mutate(cors = map(data, corrr::correlate, method = "spearman"),
         stretch = map(cors, corrr::stretch)) %>%
  unnest(stretch) %>% 
  remove_missing() %>% 
  group_nest(class) %>% 
  mutate(data2 = map(data, ~ pivot_wider(.x, names_from = "y", values_from = "r")))


# Figure 2E, S12D. Subclass proportion correlations
for (class1 in prop_cor$class) {
  i <- match(class1, prop_cor$class)
  prop1 <- as.data.frame(prop_cor$data2[[i]][,-c(1:3)])
  row.names(prop1) <- prop_cor$data2[[i]]$x
  # prop1[is.na(prop1)] <- 1
  hc1 <- hclust(dist(prop1))
  hc_order <- hc1$labels[hc1$order]
  prop1 <- prop1[hc_order, hc_order]
  
  pdf(file = paste0(output_dir, "subclass_prop_cor_", class1, ".pdf"), 
      width = 4, height = 4)
  pheatmap(prop1,
           cluster_rows = FALSE, cluster_cols = FALSE,  
           main = class1)
  dev.off()
}



#### Figure 3 ####
meta <- readRDS(file = "data/meta_ss_3_29_22.RDS")

ap_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG", "V1")
ap_order <- intersect(ap_order, unique(meta$region))

ei_calc <- meta %>% 
  filter(! layer %in% c("L1", "WM")) %>%
  mutate(region = fct_relevel(region, ap_order)) %>% 
  mutate(layer2 = str_replace(layer, "[abc]+", ""), 
         layer2 = str_replace(layer2, "L[23]$", "L2_3")) %>% 
  group_by(region, layer2, layer, donor, class) %>% 
  count() %>% 
  filter(n >= 10) %>% 
  group_by(region, layer2, donor) %>% 
  mutate(ei_ratio2 = sum(n[class == "excitatory"]) / 
           sum(n[class == "inhibitory"])) %>% 
  group_by(region, layer2, ei_ratio2, layer, donor) %>% 
  summarize(ei_ratio = n[class == "excitatory"] / 
              n[class == "inhibitory"])

write_csv(ei_calc, file = paste0(output_dir, "ei_ratio.csv"))


# Figure 3D. By region
region_pal <- c("DLPFC" = "mediumblue", "ACC" = "gold2", 
                "M1" = "firebrick", "S1" = "salmon2", "A1" = "mediumpurple1",
                "MTG" = "paleturquoise3", "V1" = "slateblue4")  # "ANG" = "springgreen4"
region_pal <- region_pal[ap_order] 


p1 <- ei_calc %>% 
  ggplot(aes(x = ei_ratio, y = fct_rev(layer), fill = fct_rev(region))) +
  geom_boxplot() +
  scale_fill_manual(values = region_pal) +
  facet_wrap(~ region, ncol = 1, scales = "free_y") +
  labs(x = "E:I ratio", y = "Layer", fill = "Region")
plot(p1)

ggsave(p1, filename = str_c(output_dir, "ei_ratio_all_layers.pdf"), width = 4, height = 8)


# Figure 3B. By layer
p2 <- ei_calc %>% 
  ggplot(aes(x = ei_ratio2, y = fct_rev(region), fill = fct_rev(region))) +
  geom_boxplot() +
  scale_fill_manual(values = region_pal) +
  facet_wrap(~ layer2, ncol = 1) +
  labs(x = "E:I ratio", y = "", fill = "Region") +
  guides(fill = "none")
plot(p2)

ggsave(p2, filename = str_c(output_dir, "ei_ratio_simple_layers.pdf"), width = 3, height = 6)



#### Figure 4, S10 ####

# Load distances on cortical sheet
areal_dist <- read_csv("data/areal_dist.csv")

to_plot <- areal_dist %>%
  gather(key = "region_2",
         value = "dist",
         -region_1) %>%
  filter(!(is.na(dist))) %>%
  mutate(lookup = str_c(region_1, "_", region_2)) %>%
  mutate(dist = as.numeric(dist)) %>%
  left_join(
    all_results %>%
      mutate(lookup = str_c(region_1, "_", region_2)) %>%
      select(-c(region_1, region_2)),
    by = "lookup"
  )


# Add slopes
to_plot2 <- to_plot %>% 
  filter(! subclass %in% c("VLMC", "Endo")) %>% 
  mutate(
    pw_order = as.numeric(fct_reorder(lookup, -spearman_corr, mean)), 
    lookup = fct_reorder(lookup, pw_order), 
    v1_comp = grepl("V1", lookup)
  ) %>% 
  group_by(cellclass) %>% 
  nest() %>% 
  mutate(
    model = map(data, ~lm(spearman_corr ~ dist + v1_comp*subclass, data = .x)),
    cor_pred = map2(model, data, predict), 
    slope = map_dbl(model, ~ coef(.)["dist"]), 
    # int = map(model, ~ coef(.))
    int_v1 = map2(model, data, ~ as.numeric(predict(.x, newdata = data.frame(dist = rep(0, nrow(.y)), v1_comp = rep(TRUE, nrow(.y)), subclass = .y$subclass)))),
    int_nonv1 = map2(model, data, ~ as.numeric(predict(.x, newdata = data.frame(dist = rep(0, nrow(.y)), v1_comp = rep(FALSE, nrow(.y)), subclass = .y$subclass))))
    # int_v1 = map2(model, data, ~ predict(.x, newdata = data.frame(dist = 0, v1_comp = TRUE, subclass = .y$subclass))),
    # int_nonv1 = map2(model, data, ~ predict(.x, newdata = data.frame(dist = 0, v1_comp = FALSE, subclass = .y$subclass))),
    # v1_offset = int_v1 - int_nonv1
  ) %>% 
  unnest(c(data, cor_pred, int_v1, int_nonv1)) %>% 
  mutate(v1_offset = int_v1 - int_nonv1) %>% 
  select(-model) %>%
  ungroup()



# Plot classes separately, dist by class, color by V1, order subclass by V1 offset
classes <- c("Excitatory", "Inhibitory", "Non-neuronal")
p1 <- list()
for (class1 in classes) {
  p1[[class1]] <- to_plot2 %>%
    filter(cellclass == class1) %>%
    ggplot(aes(x = dist, y = spearman_corr, color = v1_comp)) +
    geom_point() +
    geom_line(aes(x = dist, y = cor_pred), size = 2) +
    facet_wrap( ~ fct_reorder(subclass, v1_offset), nrow = 1) +
    labs(x = "Distance", y = "Spearmann correlation")
}

# Figure 4E
pdf(file = str_c(output_dir, "subclass_cor_vs_dist_exc.pdf"),
    useDingbats = FALSE,
    height = 4, width = 9)
print(p1[[1]])
dev.off()

# Figure S10C
pdf(file = str_c(output_dir, "subclass_cor_vs_dist_inh.pdf"),
    useDingbats = FALSE,
    height = 4, width = 9)
print(p1[[2]])
dev.off()

# Figure S10D
pdf(file = str_c(output_dir, "subclass_cor_vs_dist_nn.pdf"),
    useDingbats = FALSE,
    height = 4, width = 5)
print(p1[[3]])
dev.off()



# Load expression means
expr_mean <- list()
expr_mean[["IT"]] <- readRDS(file = str_c(dat_dir, "it_types_means.RDS"))
expr_mean[["NonIT"]] <- readRDS(file = str_c(dat_dir, "deep_exc_means.RDS"))
expr_mean[["MGE"]] <- readRDS(file = str_c(dat_dir, "mge_inh_means.RDS"))
expr_mean[["CGE"]] <- readRDS(file = str_c(dat_dir, "cge_inh_means.RDS"))
expr_mean[["NN"]] <- readRDS(file = str_c(dat_dir, "glia_means.RDS"))

ap_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG", "V1")

expr_all <- NULL
for (class1 in names(expr_mean)) {
  subclass_mean <- NULL
  expr1 <- expr_mean[[class1]]
  expr1 <- expr1[ap_order]
  for (region1 in names(expr1)) {
    for (subclass1 in colnames(expr1[[region1]])) {
      mean1 <- data.frame(
        region = region1,
        neighborhood = class1,
        subclass = subclass1, 
        gene = row.names(expr1[[region1]]),
        expr = as.numeric(expr1[[region1]][, subclass1])
      )
      subclass_mean <- rbind(subclass_mean, mean1)
    }
  }
  expr_all <- rbind(expr_all, subclass_mean)
}

expr_all <- as_tibble(expr_all) %>% 
  filter(! subclass %in% c("VLMC", "Endo")) %>%
  mutate(region = fct_relevel(region, ap_order))


# Calc A-P, spec
set.seed(123)

calc_grad <- function(x, y) {
  cor1 <- cor(as.numeric(x), as.numeric(y), method = "spearman")
  if(is.na(cor1)) cor1 <- 0
  return(cor1)
}

expr_cor <- expr_all %>%
  pivot_wider(names_from = region, values_from = expr) %>%
  group_by(neighborhood, subclass) %>%
  nest() %>%
  mutate(
    cor_ap = map(data, ~ apply(.x[,-1], 1, function(y) {
      region_id <- 1:8
      calc_grad(y[region_id], region_id)
    })),
    cor_ap2 = map(data, ~ apply(.x[,-1], 1, function(y) {
      region_id <- 2:7
      calc_grad(y[region_id], region_id)
    })),
    cor_perms = map(data, ~ apply(.x[,-1], 1, function(y) {
      region_id <- 1:8
      region_id_perm <- sample(region_id, length(region_id), replace = FALSE)
      sapply(1:2, function(perm) {
        calc_grad(y[region_id], region_id_perm)
      })
    })),
    cor_perm = map(cor_perms, ~ rowMeans(apply(.x, 1, sort)))
  )


expr_calc <- expr_cor %>%
  mutate(
    expr_max = map(data, ~ apply(.x[,-1], 1, max)),
    expr_sd = map(data, ~ apply(.x[,-1], 1, sd)),
    region_max = map(data, ~ colnames(.x[, -1])[apply(.x[,-1], 1, which.max)]),
    expr_tau = map(data, ~ calc_tau(.x[, -1]))
  ) %>% 
  select(-cor_perms) %>% 
  unnest(c(data, expr_max, expr_sd, region_max, expr_tau, cor_ap, cor_ap2, cor_perm)) %>%
  group_by(subclass) %>%
  mutate(ks_stat = ks.test(abs(cor_ap), abs(cor_perm), alternative = "less")$p, 
         ap_sig = abs(cor_ap) > 0.7 & abs(cor_ap2) > 0.5 & 
           sign(cor_ap) == sign(cor_ap2) & expr_sd > 1, 
         cor_dir = recode_factor(sign(cor_ap), `-1` = "A-P", `0` = "none", `1` = "P-A")) %>%
  ungroup() %>% 
  mutate(neighborhood = factor(neighborhood,
                               levels = c("IT", "NonIT", "MGE", "CGE", "NN")))


# Figure S10A. Tau by region
p1 <- expr_calc %>% 
  filter(expr_max > 1) %>% 
  ggplot(aes(x = fct_reorder(region_max, expr_tau), y = expr_tau)) +
  geom_boxplot() +
  labs(x = "", y = "Expression specificity (tau)") +
  guides(x = guide_axis(angle = -45))

ggsave(p1, filename = str_c(output_dir, "expr_tau_region.pdf"), width = 3, height = 4)


# Figure S10B. Tau by subclass
p2 <- expr_calc %>% 
  filter(expr_max > 1) %>% 
  ggplot(aes(x = fct_reorder(subclass, expr_tau), y = expr_tau)) +
  geom_boxplot() +
  labs(x = "", y = "Expression specificity (tau)") +
  guides(x = guide_axis(angle = -45))

ggsave(p2, filename = str_c(output_dir, "expr_tau_subclass.pdf"), width = 6, height = 4)


# Figure 4D. Top tau genes
p4 <- expr_calc %>% 
  filter(expr_max > 1) %>% 
  filter(expr_tau > 0.8) %>% 
  # select(gene) %>% distinct(gene) %>% arrange(gene) %>% View()  # GO analysis
  group_by(subclass, region_max) %>% 
  count() %>% 
  group_by(subclass) %>% 
  mutate(n2 = sum(n)) %>% 
  ggplot(aes(x = n, y = fct_reorder(subclass, n2), 
             fill = fct_reorder(region_max, -n, sum))) +
  geom_bar(stat = "identity") +
  labs(x = "Number of genes", 
       y = "") +
  guides(fill = guide_legend("Region"))

ggsave(p4, filename = str_c(output_dir, "region_specific_count_bysubclass.pdf"), width = 4, height = 3)


# Figure S10I. Areal gradient histograms by subclass - SI
p5 <- expr_calc %>% 
  pivot_longer(cols = c(cor_ap, cor_perm), names_to = "cor_type", values_to = "cor") %>% 
  ggplot(aes(x = cor, color = cor_type)) + 
  # stat_ecdf(aes(x = abs(cor), color = cor_type)) +
  geom_density(adjust = 1.5) +
  facet_wrap(~ fct_reorder(subclass, ks_stat))

ggsave(p5, filename = str_c(output_dir, "ap_grad_densities_bysubclass.pdf"), width = 10, height = 10)


# Figure S10J. Gradient correlations
p6 <- to_plot %>% 
  ggplot(aes(x = gene, xend = gene, y = ap_min, yend = ap_max)) +
  geom_hline(yintercept = 0) +
  geom_segment(color = "grey") +
  geom_point(aes(x = gene, y = ap_mean, size = ap_cnt)) +
  geom_segment(data = to_plot %>% filter(ap_switch == TRUE),
               color = "red") +
  geom_point(data = to_plot %>% filter(ap_switch == TRUE),
             aes(x = gene, y = ap_mean, size = ap_cnt), 
             color = "red") +
  scale_size_continuous(range = c(0.25, 3)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x = "Gene", y = "A-P correlation") +
  guides(x = guide_axis(angle = 90))
plot(p6)

ggsave(p6, filename = str_c(output_dir, "ap_cor_range.pdf"), width = 6, height = 3)


# Figure 4G. Top A-P genes
p6 <- expr_calc %>%
  filter(ap_sig == TRUE) %>% 
  # filter(expr_tau < 0.8) %>%  # Exclude region-specific genes
  # select(gene) %>% distinct(gene) %>% arrange(gene) %>% View()  # GO analysis
  group_by(subclass, cor_dir) %>% 
  count() %>% 
  group_by(subclass) %>% 
  mutate(n2 = sum(n)) %>% 
  # filter(n > 0) %>% arrange(-n) %>% View()
  ggplot(aes(x = n, y = fct_reorder(subclass, n2), fill = cor_dir)) +
  geom_bar(stat = "identity") +
  labs(x = "Number of genes", 
       y = "") +
  guides(fill = guide_legend("Direction"))

ggsave(p6, filename = str_c(output_dir, "ap_grad_count_bysubclass.pdf"), width = 4, height = 3)


# Figure 4H. Plot A-P for example genes
plot_list <- list(subclass = list(), gene = list())
i <- 1; plot_list$subclass[[i]] <- c("L2/3 IT", "L4 IT"); plot_list$gene[[i]] <- "CBLN2"  # match fetal grad
i <- 2; plot_list$subclass[[i]] <- c("L5 IT"); plot_list$gene[[i]] <- c("CNTN5", "CNTN6")  # opp. grad related genes
i <- 3; plot_list$subclass[[i]] <- c("L5/6 NP", "Vip"); plot_list$gene[[i]] <- c("DCC")  # opposing grad

for (i in 1:length(plot_list$subclass)) {
  p1 <- expr_all %>% 
    filter(subclass %in% plot_list$subclass[[i]] & gene %in% plot_list$gene[[i]]) %>%  
    ggplot(aes(x = as.numeric(region), y = expr, color = subclass)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~ gene, nrow = 1, scales = "free_y") +
    labs(x = "", y = "Normalized expression") +
    scale_x_continuous(breaks = 1:8, labels = ap_order, guide = guide_axis(angle = 90)) +
    theme(panel.grid.minor = element_blank())
  
  plot(p1)
  ggsave(p1, filename = paste0(output_dir, "subclass_gene_trend_", i, ".pdf"), 
         width = 2*length(plot_list$gene[[i]]) + 1, height = 2.5)
}




# Subclass order
subclass <- tibble(celltype = c("L23IT","L4IT","L5IT","L6IT","L6ITCar3",
                                "L5ET","L56NP","L6CT","L6b",
                                "Lamp5Lhx6","Lamp5","Sncg","Vip","Pax6",
                                "Chandelier","Pvalb","Sst","SstChodl",
                                "OPC","Astro","Oligo","Endo","MicroPVM","VLMC"), 
                   neighborhood = c(rep("IT", 5), 
                                    rep("Non-IT", 4), 
                                    rep("CGE", 5), 
                                    rep("MGE", 4), 
                                    rep("NN", 6))
)
subclass_order <- subclass$celltype


dat <- read_csv(file = "data/all_variable_model.csv") %>%
  select(-row) %>%
  left_join(subclass, by = "celltype") %>% 
  filter(! celltype %in% c("Endo", "VLMC")) %>% 
  rename(
    medial_lateral = x,
    rostral_caudal = y,
    dorsal_ventral = z
  ) %>%
  pivot_longer(cols = donor:Residuals,
               names_to = "covar",
               values_to = "expr_var") %>%
  mutate(covar = fct_relevel(covar,
                             c("region", "rostral_caudal", "medial_lateral", "dorsal_ventral",
                               "donor", "region_within_area_cluster", "Residuals")
  ))


# Figure S10F. Var gene count barplot
p3 <- dat %>% 
  # filter(neighborhood %in% c("CGE", "MGE") & covar %in% c("rostral_caudal", "medial_lateral", "dorsal_ventral", "region")) %>%  # Inh only
  mutate(covar_type = ifelse(covar %in% c("rostral_caudal", "medial_lateral", 
                                          "dorsal_ventral", "region"), "Regional", "Other")) %>% 
  mutate(celltype = factor(celltype, levels = subclass_order)) %>% 
  group_by(celltype, covar, covar_type) %>% 
  summarize(var_genes = sum(expr_var > 0.05)) %>% 
  filter(! covar %in% c("Residuals")) %>%
  filter(var_genes > 0) %>% 
  #fct_reorder(celltype, var_genes, sum)
  ggplot(aes(x = celltype, y = var_genes, fill = covar)) + 
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~ covar_type, ncol = 1, scale = "free_y") +
  labs(x = "", y = "Number of genes") +
  guides(x = guide_axis(angle = 90)) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())
plot(p3)

ggsave(p3, filename = "output/var_prop_barplots_compact.pdf", width = 6, height = 3.5)


# Figure S10F. Median variance explained
p4 <- dat %>% 
  mutate(covar_type = ifelse(covar %in% c("rostral_caudal", "medial_lateral", 
                                          "dorsal_ventral", "region"), "Regional", "Other")) %>% 
  # filter(covar %in% c("rostral_caudal", "medial_lateral", "dorsal_ventral", "region")) %>%
  # filter(covar %in% c("donor", "region_within_area_cluster")) %>%
  filter(! covar %in% c("Residuals")) %>%
  mutate(celltype = factor(celltype, levels = subclass_order)) %>%
  filter(expr_var > 0.05) %>%
  group_by(celltype, covar, covar_type) %>% 
  summarize(var_median = median(expr_var), 
            var_genes = sum(expr_var > 0.05)) %>% 
  ggplot(aes(x = celltype, y = var_median, 
             color = covar)) + 
  geom_point() +
  # geom_line() +
  scale_color_brewer(palette = "Dark2") +
  # scale_color_manual(values = brewer.pal(4, "Dark2")) +
  # scale_color_manual(values = brewer.pal(6, "Dark2")[5:6]) +
  facet_wrap(~ covar_type, ncol = 1) +
  labs(x = "", y = "Median var.") +
  guides(x = guide_axis(angle = 90), color = "none") +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())
plot(p4)

ggsave(p4, filename = "output/var_summary.pdf", width = 6, height = 3.5)


# Figure S10H. Variance explained by gradients - subclass ternary plots
to_plot <- dat %>% 
  filter(covar %in% c("rostral_caudal", "medial_lateral", "dorsal_ventral")) %>% 
  # filter(celltype == "L4IT") %>% 
  group_by(celltype, gene) %>% 
  mutate(var_sum = sum(expr_var), 
         var_max = max(expr_var), 
         var_prop = expr_var/var_sum) %>% 
  filter(var_max > 0.05) %>% 
  pivot_wider(id_cols = c(celltype, neighborhood, gene, var_sum), names_from = covar, values_from = var_prop) %>%
  group_by(celltype) %>% 
  mutate(var_cnt = n_distinct(gene)) %>% 
  ungroup() %>% 
  mutate(celltype = factor(celltype, levels = intersect(subclass_order, unique(celltype))))

# By subclass 
p4 <- to_plot %>% 
  # filter(var_cnt > 5) %>%
  ggtern(aes(x = medial_lateral, y = rostral_caudal, z = dorsal_ventral, 
             size = var_sum)) +
  geom_point(shape = 1, color = "dark blue") +
  scale_size_continuous(limits=c(0, 1), breaks=seq(0, 0.6, by=0.15)) +
  labs(x = "M-L", y = "R-C", z = "D-V") +
  facet_wrap(~ celltype, ncol = 7) +
  theme_bw(base_size = 12)
plot(p4)

ggsave(p4, filename = "output/var_prop_ternary.pdf", width = 20, height = 16)



# Figure 4F. Variance explained by gradients - summary
p5 <- to_plot %>% 
  group_by(celltype, neighborhood, var_cnt) %>% 
  summarize(ml_mean = weighted.mean(medial_lateral, var_sum), 
            rc_mean = weighted.mean(rostral_caudal, var_sum), 
            dv_mean = weighted.mean(dorsal_ventral, var_sum)) %>% 
  ggtern(aes(x = ml_mean, y = rc_mean, z = dv_mean, 
             size = var_cnt, 
             color = factor(neighborhood, 
                            levels = c("IT", "Non-IT", "MGE", "CGE", "NN")))) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "M-L", y = "R-C", z = "D-V") +
  guides(size = guide_legend("No. of genes"), 
         color = guide_legend("Neighborhood")) +
  theme_bw()
plot(p5)

ggsave(p5, filename = "output/var_prop_ternary_summary.pdf", width = 5, height = 5)



#### Figure 5 ####

dend <- readRDS(file = "data/dend.RDS")
cl_order <- labels(dend)
ap_order <- c("ACC", "DLPFC", "M1", "S1", "MTG", "A1", "ANG", "V1")


dat_n <- read_csv(file = "data/neurons_region_cross_area_cluster_effect_size_p_gt_0.8.csv")
dat_g <- read_csv(file = "data/glia_region_cross_area_cluster_effect_size_p_gt_0.8.csv")

dat <- bind_rows(dat_n, dat_g) %>% 
  mutate(`Cell Type` = factor(`Cell Type`, levels = cl_order), 
         Covariate = factor(Covariate, levels = ap_order))



# Figure 5F. Areal differences in proportions
effect_th <- 1
p2 <- dat %>%
  mutate(prop_change = cut(`Effect Size`, breaks = c(-10, -1, -0.1, 0.1, 1, 10), 
                           labels = c("Large decr.", "Small decr.", "No change", 
                                      "Small incr.", "Large incr.")), 
         prop_change = fct_rev(prop_change)) %>%
  ggplot(aes(x = `Cell Type`, y = fct_rev(Covariate), fill = prop_change)) +
  geom_tile() +
  scale_fill_manual(values = brewer.pal(7, "RdBu")[c(1,3,4,5,7)]) +
  labs(x = "", y = "", fill = "Prop. change") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  guides(x = guide_axis(angle = 90))
p2

ggsave(p2, filename = str_c(output_dir, "sccoda_type_prop_change_thresh.pdf"), width = 16, height = 1.5)




#### Figure 6 ####
meta <- readRDS(file = "data/meta_10x_2_11_22.RDS")
dend <- readRDS(file = "data/v1_dend.RDS")

within_vs_cross <- meta %>% 
  group_by(region, within_area_cluster, cross_area_cluster) %>% 
  summarize(count = n()) %>% 
  group_by(region, within_area_cluster) %>% 
  mutate(within_cl_size = sum(count)) %>% 
  group_by(cross_area_cluster, region) %>% 
  mutate(region_cross_cl_size = sum(count)) %>%
  group_by(cross_area_cluster) %>% 
  mutate(cross_cl_size = sum(count)) %>%
  ungroup() %>% 
  mutate(cross_cl_prop = count / cross_cl_size, 
         within_cl_prop = count / within_cl_size, 
         cross_cl_specificity = region_cross_cl_size / cross_cl_size) %>%
  group_by(region, within_area_cluster) %>% 
  mutate(within_cl_specificity = sum(within_cl_prop * cross_cl_specificity)) %>% 
  ungroup()
select(-contains("size"))

write_csv(within_vs_cross, file = paste0(output_dir, "within_cross_area_cluster_mapping.csv"))


# Figure 6A. Plot V1 specificity
v1_dend_order <- labels(dend)

v1_spec <- within_vs_cross %>% 
  filter(region == "V1") %>% 
  select(within_area_cluster, within_cl_specificity) %>% 
  unique() %>% 
  arrange(match(within_area_cluster, v1_dend_order))


pdf(file = paste0(output_dir, "v1_spec_bar.pdf"), width = 6, height = 3)
barplot(v1_spec$within_cl_specificity, ylim = c(0, 1), col = "slateblue4", las = 1)
dev.off()




#### Figure 1 ####