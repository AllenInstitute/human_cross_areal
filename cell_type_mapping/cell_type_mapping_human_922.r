####################### load all data and process #####################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(gridExtra)
  library(Matrix)
  library(matrixStats)
  library(scrattch.hicat)
  library(tibble)
  library(patchwork)
  library(dplyr)
  library(MASS)
  library(viridis)
  library(dplyr)
  library(cowplot)
  library(ggthemes)
  library(data.table)
  library(doMC)
  library(anndata)
  library(uwot)
  library(spdep)
  library(feather)
})

options(stringsAsFactors = FALSE)
options(scipen=999)
options(mc.cores = 20)

########################## code to process cell by gene table #######################

### Placeholders ###
inputFolder <- SOURCE_DATA_FOLDER 
run_tracker <- HUMAN_IMAGING_DATA_TRACKER
# load files for mapping by loading map_knn script from https://github.com/AllenInstitute/bioinformatics/blob/main/Mapping/map_knn.R
metadata_rnaseq <- ANNOTATION_DATA 
cl.means <- READ_COUNT_MATRIX # read count matrix for human MTG rnaseq data

run_tracker_human <- run_tracker %>% 
  dplyr::filter(Species == "Human") %>% 
  drop_na(`Experiment Name`) 

min_genes <- 3
min_total_reads <- 30
min_vol <- 100

upper_bound_reads <- 4000
upper_genes_read <- 130

gene_panel <- "VZG167a"

#filter out smart-seq data from metadata
metadata_rnaseq <- metadata_rnaseq %>% filter(species_tech=="human_10x")
  
cl          = metadata_rnaseq$cluster
names(cl)   = metadata_rnaseq$sample_id

# generate anno file
train.cl.df <- metadata_rnaseq %>% 
  dplyr::select(cluster,
                cluster_color,
                subclass,
                subclass_color,
                neighborhood,
                cross_species_cluster,
                cross_species_cluster_color,
                class) %>% 
  distinct(cluster, .keep_all = TRUE)

folder <- list.dirs(inputFolder, full.names = TRUE, recursive = FALSE)

counter <- 0 

for (inputFolder in folder) {
  file_name <- basename(inputFolder)
  file_name_parts <- str_split(file_name, pattern = "_")
  experiment <- file_name
  idx <- which(run_tracker_human$'Run Name' == file_name | run_tracker_human$'Experiment Name' == file_name | run_tracker_human$'Experiment Name' == file_name_parts[[1]][2])
  if (length(idx) == 0) {
	idx <- counter
  }
  species <- 'Human'
  section <- run_tracker_human[[idx, 1]]
  merscope <- run_tracker_human[[idx, 32]]
  collection_year <- substring(run_tracker_human[[idx, 30]],1,4)
  source <- substring(run_tracker_human[[idx, 45]],4,5)
  
  #cell by gene table, from experiment
  cbg <- read.csv(paste0(inputFolder,"/region_0/cell_by_gene.csv"), 
				  header=TRUE, 
				  row.names = 1,
				  check.names = FALSE)
  
  #cell metadata, from experiment
  metadata <- read.csv(paste0(inputFolder,"/region_0/cell_metadata.csv"), 
					   header=TRUE,
					   row.names = 1,
					   check.names = FALSE)
  
  #clean and combine cell data
  metadata <- metadata[match(rownames(cbg),rownames(metadata)),]
  
  blanks <- dplyr::select(cbg,contains("Blank"))
  cbg <- dplyr::select(cbg,-contains("Blank"))
  
  metadata$genes_detected <- rowSums(cbg!=0)
  metadata$total_reads <- rowSums(cbg)
  
  dir.create(file.path(paste0(inputFolder,"/processed/")), showWarnings = FALSE)
  
  metadata$species <- species
  metadata$collection_year <- collection_year
  metadata$source <- source
  metadata$merscope <- merscope
  metadata$gene_panel <- gene_panel
  metadata$section <- section

  upper_bound_area <- 3*(median(metadata$volume[metadata$volume>100]))

  #calculate qc metrics for cell detection
  metadata <- metadata %>% mutate(cell_qc = if_else(genes_detected < min_genes |
													total_reads < min_total_reads|
													volume < min_vol |
													volume > upper_bound_area |
													genes_detected > upper_genes_read |
													total_reads > upper_bound_reads,
													"Low","High"))

  metadata$min_genes <- min_genes
  metadata$min_total_reads <- min_total_reads
  metadata$min_vol <- min_vol
  metadata$upper_bound_area <- upper_bound_area
  metadata$upper_bound_reads <- upper_bound_reads
  metadata$upper_genes_read <- upper_genes_read
  
  #prepare for mapping
  cbg <- as.matrix(cbg)
  cbg <- Matrix(cbg, sparse = TRUE)
  #normalize counts per volume
  cbg_cpum <- (cbg / metadata$volume*1000)
  cbg_norm <- log2(cbg_cpum+1)

  # start mapping
  vizgen.dat <- t(cbg_cpum)
  vizgen.dat <- log2(vizgen.dat+1)

  genes <- rownames(vizgen.dat)
  useGenes <- intersect(rownames(cl.means),genes)
  vizgen.dat = vizgen.dat[useGenes,]
  train.cl.dat <- cl.means[useGenes,]
  
  index.bs = build_train_index_bs(train.cl.dat, method="cor",fn = "fb.index") 
  map.result = map_cells_knn_bs(vizgen.dat, train.index.bs=index.bs, method="cor", mc.cores=20) #from map_knn script mentioned at top
  best.map.df = map.result$best.map.df
  cl.anno = best.map.df %>% left_join(train.cl.df,by = c("best.cl"="cluster"))
   
  cl.list = with(cl.anno, list(cl=setNames(best.cl, sample_id),subclass=setNames(subclass, sample_id),subclass=setNames(best.cl, sample_id)))
  z.score=z_score(cl.list, val=with(best.map.df, setNames(avg.cor,sample_id)), min.sample=100)
  cl.anno$z.score = z.score[cl.anno$sample_id]

  #add metadata to mapped cells
  anno_mfish <- merge(metadata,
					 cl.anno,
					 by.x = 0,
					 by.y="sample_id"
					 )

  anno_mfish <- column_to_rownames(anno_mfish, var = "Row.names")
  anno_mfish <- anno_mfish[match(rownames(cbg_cpum),rownames(anno_mfish)),]
  colnames(anno_mfish)[29] <- "cluster"

  #filter cells with low qc metrics
  to_keep <- intersect(rownames(cbg),rownames(subset(anno_mfish,cell_qc %in% "High")))
  cbg_filtered <- cbg_norm[to_keep,]
  cbg_filtered <- as.matrix(cbg_filtered)

  #generate coordinates for cirro mapping
  coordinates <- anno_mfish %>%
	dplyr::filter(cell_qc == "High") %>%
	dplyr::select(center_x,center_y)

  coordinates_cirro <- anno_mfish %>%
	dplyr::filter(cell_qc == "High") %>%
	dplyr::select(center_x,center_y)

  coordinates_cirro$center_y = -(coordinates_cirro$center_y - mean(coordinates_cirro$center_y))
  coordinates_cirro$center_x = coordinates_cirro$center_x - min(coordinates_cirro$center_x)
  coordinates_cirro$center_x = coordinates_cirro$center_x + (counter*10000)
  counter = counter+1

  # prepare data for anndata format
  uns <- c("species","merscope","gene_panel","min_genes","min_total_reads","min_vol","upper_bound_area","upper_bound_reads","upper_genes_read")
 
  anno_mfish_subset <- anno_mfish %>%
	dplyr::filter(cell_qc == "High") %>%
	dplyr::select(cluster,
				  subclass,
				  neighborhood,
				  class, 
				  merscope,
				  avg.cor,
				  genes_detected,
				  total_reads,
				  prob,
				  volume)
  umap_mfish <- umap(cbg_filtered,n_neighbors = 25,n_components = 2,metric = "euclidean",min_dist = 0.4,pca = 50)

  blanks_filtered <- log2(blanks[to_keep,]+1)
  
  #save anndata file
  ad <- AnnData(
	X = cbg_filtered,
	obs = anno_mfish_subset,
	obsm = list(
	  blanks = as.matrix(blanks_filtered),
	  spatial = as.matrix(coordinates),
	  spatial_cirro = as.matrix(coordinates_cirro),
	  X_umap = umap_mfish
	  ),
	uns = list(
	  species = species,
	  section = section,
	  merscope = merscope,
	  gene_panel = gene_panel,
	  min_genes = min_genes,
	  min_total_reads = min_total_reads,
	  min_vol = min_vol,
	  upper_bound_area = upper_bound_area,
	  upper_bound_reads = upper_bound_reads,
	  upper_genes_read = upper_genes_read
	  )
	)
	write_h5ad(ad,paste0(inputFolder,"/processed/",file_name,".h5ad"))
}
