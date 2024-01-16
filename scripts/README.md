## Scripts

seurat_object_merging_function:

This is a generic function for taking a Seurat object, clustering it into meta cells, then merging those meta cells into clusters with specified criteria. This approach has been implemented to identify defined neuronal and non-neuronal (glial) clusters in recent publications by the Allen Institute.

General steps:
	1.	Place dataset into a Seurat object
	2.	Create a reduced dimension embedding on which to perform clustering and nearest neighbor comparisons 
	⁃	Note: embedding can be an integrated or non-integrated PCA space 
	3.	Over-cluster the data (above the number of cell types you expect) to produce small clusters, or meta cells
	4.	(OPTIONAL- for speed) downsample the dataset to a maximum number of cells per meta cell
	5.	Set merging criteria thresholds - 
	⁃	n_cells_thresh = minimum number of cells required in a meta cell (else merge), default 20 cells
	⁃	n_degs_thresh = minimum number of DEGs that must exist between meta cell and nearest neighbor (else merge), default 8 DEGs
	6.	DEG criteria - 
	⁃	Differentially expressed genes for this use case are defined as being expressed in >50% of a meta cell, have 2 times the normalized expression of the comparator meta cell, and have a proportion expressed differential of 0.3 or greater

Considerations:
	•	Parameters will need to be adjusted based on the number of cells in your dataset and how deeply they are sequenced.
	•	We used an integrated embedding between SMARTseq and 10x, but used the 10x data for meta cell comparisons and merging logic.
	•	We used the merging function after subsetting the dataset into 5 neighborhoods, integrating each neighborhood across modalities/donors, then overclustering each neighborhood to ~100 meta cells.
	•	Glial subtypes require a lower number of DEGs to resolve subtypes (ie different types of microglia or astrocytes) compared to neurons. 
	•	After applying this approach, assessing cluster quality and composition are still necessary steps.