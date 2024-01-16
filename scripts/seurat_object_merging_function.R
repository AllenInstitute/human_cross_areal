library(Seurat)
library(tidyverse)
library(scrattch.hicat)

seurat_obj <- seurat_obj %>% 
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeature = 3000) %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors()

seurat_obj <- FindClusters(seurat_obj, resolution = 6) #over-cluster the data in 'metacells'

results <- merge_clusters(seurat_obj, starting_clusters = "seurat_clusters", reduction = "pca")

#Run the below code to activate the merging function...
# consider which reduction you are using
# consider changing lines 175-179 to modify DEG criteria

#merge cluster function
merge_clusters <- function(seurat_obj, n_cells_thresh = 20, n_degs_thresh = 8, ds_to_n_cells = 200, 
                           starting_clusters, reduction){
  # n_cells_thresh <- 20 #min cluster size (number of cells) that must exist in cluster of interest to prevent merge
  # n_degs_thresh <- 8 #number of DEGs between cluster of interest and NN that must exist to prevent merge
  # ds_to_n_cells <- 200 #number of cells to downsample to
  # starting_clusters <- "seurat_clusters" #some field in the seurat object metadata that identifies clusters to merge
  # reduction <- "harmony" #PC space to find NN over.. options are "harmony" or "pca"
  
  #downsample dataset to max of 200 cells per cluster
  Idents(seurat_obj) <- seurat_obj@meta.data %>% 
    as_tibble(rownames = "cell_ids") %>%
    dplyr::select(cell_ids, starting_clusters) %>%
    deframe()
  
  seurat_obj_sub <- subset(seurat_obj, downsample = ds_to_n_cells)
  
  #lookup table to store original cluster assignments and iteratively merged cl assignments
  lookup_table <- seurat_obj_sub@meta.data %>% 
    as_tibble(rownames = "cell_ids") %>%
    dplyr::select(cell_ids, starting_clusters) %>%
    mutate(current_cl_assignment = get(starting_clusters)) %>%
    mutate_all(as.character)
  
  #full PC space with cell assignments to find NN for merging
  if(reduction == "harmony"){
    harmony_pcs <- seurat_obj_sub@reductions$harmony@cell.embeddings %>% 
      as_tibble(rownames = "cell_ids") %>%
      dplyr::select(cell_ids, str_c("harmony_", 1:30)) 
  }
  
  if(reduction == "pca"){
    PCA_pcs <- seurat_obj_sub@reductions$pca@cell.embeddings %>% 
      as_tibble(rownames = "cell_ids") %>%
    dplyr::select(cell_ids, str_c("PC_", 1:30)) 
  }
  
  #loop over the next steps a max of 500 times, setting a new seed on each iteration to randomly sample
  #the cluster each time. If nothing changes over the course of an iteration, keep the same cls order and progress through
  #until a change occurs. Stop the loop once all clusters are iterated over and no change has occured.
  n_cl <- lookup_table %>% #initial n_cl... will update if merge occurs
    distinct(current_cl_assignment) %>%
    nrow()
  stop_counter <- 0
  cl_counter <- 1
  seeds <- 1
  
  cls <- lookup_table %>% #initial cluster ordering
    distinct(current_cl_assignment) %>% 
    sample_n(n()) %>%
    pull() %>% 
    as.character()
  
  while(cl_counter < n_cl){ #this condition ensures current cluster ordering is held constant until a merge occurs.. 
    
    if(stop_counter > 0){ #if merge occured resample and a new random cl ordering should be used
      set.seed(seeds)
      cls <- lookup_table %>% 
        distinct(current_cl_assignment) %>% 
        sample_n(n()) %>%
        pull() %>% 
        as.character()
      
      cl_counter <- 1 #reset cluster counter to 1
      stop_counter <- 0 #reset stop counter
    }
    
    #find nearest neighbors of current clusters' means for a given cluster of interest
    
    cl_oi <- cls[cl_counter] # set to random
    
    if(reduction == "harmony"){
      current_means_tbl <- harmony_pcs %>%
        left_join(lookup_table, by = "cell_ids") %>%
        dplyr::select(current_cl_assignment, contains("harmony")) %>%
        group_by(current_cl_assignment) %>%
        summarise_all(mean) %>%
        ungroup() 
    }
    
    if(reduction == "pca"){
      current_means_tbl <- PCA_pcs %>%
        left_join(lookup_table, by = "cell_ids") %>%
        dplyr::select(current_cl_assignment, contains("PC")) %>%
        group_by(current_cl_assignment) %>%
        summarise_all(mean) %>%
        ungroup() 
    }
    
    current_means_mat <- as.matrix(current_means_tbl[ , -1])
    rownames(current_means_mat) <- current_means_tbl$current_cl_assignment
    
    nn_oi <- dist(current_means_mat, method = "euclidean") %>%
      as.matrix() %>%
      as_tibble(rownames = "group_1") %>%
      gather(key = "group_2", value = "dist", -group_1) %>%
      filter(group_1 != group_2,
             group_1 == cl_oi) %>%
      arrange(dist) %>%
      dplyr::select(group_2) %>%
      dplyr::slice(1) %>%
      pull()
    
    #print(str_c("testing ", cl_oi, " against ", nn_oi))
    #abundance testing
    #if cl_oi has less than n cells, merge with NN
    n_cells <- lookup_table %>%
      filter(current_cl_assignment == cl_oi) %>%
      nrow()
    
    if(n_cells < n_cells_thresh){
      print(str_c("n cells in ", cl_oi, " are under threshold... merging with ", nn_oi, " and testing next cluster"))
      
      lookup_table <- lookup_table %>%
        mutate(current_cl_assignment = case_when(current_cl_assignment == cl_oi ~ nn_oi,
                                                 TRUE ~ current_cl_assignment))
      
      stop_counter <- stop_counter + 1 #this condition means a merge occurred 
      
      n_cl <- lookup_table %>% #update if merge occurs
        distinct(current_cl_assignment) %>%
        nrow()
      
      seeds <- seeds + 1 #change seed for next random ordering
      
    }else{
      #print("n cells in cl_oi are over threshold... move onto DEG test")
      
      ##DEG testing
      # Differentially expressed genes for this use case are defined as being expressed 
      # in >50% of a metacell, have 2 times the normalized expression of the 
      # comparator metacell, and have a proportion expressed differential of 0.3 or greater.
      
      #pull out cell ids in cl_oi and nn_oi
      comparison_ids <- lookup_table %>%
        filter(current_cl_assignment %in% c(cl_oi, nn_oi)) %>%
        dplyr::select(cell_ids, current_cl_assignment) %>%
        deframe()
      
      #find prop stats
      cl_props <- get_cl_prop(seurat_obj_sub@assays$RNA@data[ , names(comparison_ids)], cl = comparison_ids) %>%
        as_tibble(rownames = "features") %>%
        mutate(prop_diff = get(cl_oi) - get(nn_oi)) %>%
        dplyr::rename(c("cl_oi_prop" = cl_oi, 
                        "nn_oi_prop" = nn_oi))
      
      #find mean stats
      cl_means <- get_cl_means(seurat_obj_sub@assays$RNA@data[ , names(comparison_ids)], cl = comparison_ids) %>%
        as_tibble(rownames = "features") %>%
        dplyr::rename(c("cl_oi_mean" = cl_oi, 
                        "nn_oi_mean" = nn_oi))
      
      #bind results and find n_degs
      deg_stats <- cl_props %>% 
        left_join(cl_means, by = "features") %>%
        mutate(cl_oi_deg = case_when(cl_oi_prop >= 0.5 & prop_diff >= 0.3 & cl_oi_mean > 2 * nn_oi_mean ~ TRUE,
                                     TRUE ~ FALSE),
               nn_oi_deg = case_when(nn_oi_prop >= 0.5 & prop_diff <= -0.3 & nn_oi_mean > 2 * cl_oi_mean ~ TRUE,
                                     TRUE ~ FALSE)
        ) 
      
      n_degs <- deg_stats %>% filter(cl_oi_deg) %>% nrow() +         
        deg_stats %>% filter(nn_oi_deg) %>% nrow()          
      
      
      #if n_degs < threshold, merge metacells
      if(n_degs < n_degs_thresh){
        print(str_c("n degs between ", cl_oi, " and ", nn_oi, " are under threshold... merging with NN and testing next cluster"))
        
        lookup_table <- lookup_table %>%
          mutate(current_cl_assignment = case_when(current_cl_assignment == cl_oi ~ nn_oi,
                                                   TRUE ~ current_cl_assignment))
        
        stop_counter <- stop_counter + 1 #this condition means a merge occurred 
        
        n_cl <- lookup_table %>% #update if merge occurs
          distinct(current_cl_assignment) %>%
          nrow()
        
        seeds <- seeds + 1 #change seed for next random ordering
        
        
      }else{
        #print("n degs between cl_oi and nn_oi are over threshold... testing next cluster")
        cl_counter <- cl_counter + 1 #move on to next cluster in cluster vector
        #print(str_c("cl_counter: ", cl_counter, " stop counter: ", stop_counter))
      }
    }
  }
  return(lookup_table)
}




