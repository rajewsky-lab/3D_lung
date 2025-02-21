# 3D_lung
Code avaliability for the 'Combining spatial transcriptomics and ECM imaging in 3D for mapping cellular interactions the tumor microenvironment' paper.

Recommended: Please refer to the R notebook '' to reproduce the paper analyses and figures starting from the 'preprocessed_objects' downloaded from Zenodo.

Optional: The '' python script and the 'cosmx_flat_files' (heavy), 'SHG' and 'stimwrap_files' are shared for full reproducibility.

3D data can also be interactively explored at https://lung-3d-browser.mdc-berlin.de. In the browser one can:
1. Visualise and interact with the 3D volumetric rendering of multicellular niches. To hide them, click in the dots of all the ‘Mesh’ items in the ‘Pipeline’ menu top left and leave only ‘PC_Model’ unchecked to visualise the single cells

2. Visualise different cell-level annotations that can be selected from ‘Annotation’ in the ‘Active Model’ menu. Besides x/y/z spatial coordinates, cells can also be plotted both in gene expression and neighbourhood composition UMAP space from the ‘Coordinates’ menu.

3. Visualize gene expression by selecting specific genes from the ‘Genes’ menu. Both SCT (log-normalised, recommended) and raw counts can be selected from the ‘Matrix’ menu. 

4. Subset the data according to specific metadata by selecting the relevant annotation and the ’Subset’ menu. This can be handy to e.g. plot only 1 section, one 1 or multiple cell types or only cells inside or outside the EMT niche. To go back to the full dataset just click ‘Reload'
