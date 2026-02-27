import anndata as ad

import scanpy as sc
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

# read in data
#adata = ad.read_h5ad('/lustre/compbio/projects/glioma_scRNAseq/data/Seurat_objects/hypoxia/microglia_0702.h5ad')
#adata = ad.read_h5ad('/lustre/compbio/projects/glioma_scRNAseq/data/Seurat_objects/hypoxia/macros.h5ad')
adata = ad.read_h5ad('/scratch/svc_td_cri/projects/glioma_scs_st/Narvi_glioma_scRNAseq/data/Seurat_objects/hypoxia/anndata_objs/TAM_check.h5ad')

# filter data to contain either macrophages or microglia
#adata = adata[adata.obs["preliminary_cell_types"].isin(['macrophages (M2)', 'macrophages (M1)'])]
adata = adata[adata.obs["preliminary_cell_types"].isin(['microglial cells (M2)', 'microglial cells (M1)'])]

# add 0/starting point cell to anndata, a cd163+ macrophage with lowest hypoxia score
# In Microglia this is CD163- microglia with lowest hypoxia score

# Filter cells with the celltype label 'M2 macrophage'
#m2_macrophages = adata[adata.obs['preliminary_cell_types'] == 'macrophages (M2)']
m2_macrophages = adata[adata.obs['preliminary_cell_types'] == 'microglial cells (M1)']

# Find the cell with the lowest hypoxia score
min_hypoxia_index = np.argmin(m2_macrophages.obs['hypoxia_score'])
cell_with_lowest_hypoxia = m2_macrophages[min_hypoxia_index]

# Create a new variable in the AnnData object
adata.obs['hypoxia_label'] = 'other'
adata.obs.loc[cell_with_lowest_hypoxia.obs_names, 'hypoxia_label'] = 'lowest_hypoxia_m2'

#Filter data

adata.obs["seurat_clusters"] = pd.Categorical(adata.obs["seurat_clusters"])
adata = adata[adata.obs["seurat_clusters"] != 45]

#adata = adata[adata.obs["filttered_data_set_annotation"] != "Dendritic cells"]

# neighbor graph and clustering
sc.pp.neighbors(adata, use_rep="X_harmony")
sc.tl.leiden(adata,resolution=1) # adjust res if needed
#adata = adata[adata.obs["filttered_data_set_annotation"] != "Plasmacytoid Dendritic cells"]
adata = adata[adata.obs["leiden"] != "13"]

# graph abstraction
sc.tl.draw_graph(adata, random_state=36)
#sc.pl.draw_graph(adata)
sc.pl.draw_graph(adata, color={"leiden", "hypoxia_classification", 'preliminary_cell_types'}, legend_loc="on data", save="_leiden_hypoxia1510_seed36_check_no45_nocluster13_MICROS.pdf")

#adata = adata[adata.obs["leiden"] != "15"]
#adata = adata[adata.obs["leiden"] != "12"]

sc.tl.paga(adata, groups="leiden")
sc.pl.paga(adata, color=['leiden', 'hypoxia_classification', 'preliminary_cell_types'], save="_leiden_hypoxia1510_no45_nocluster13_MICROS.pdf")

sc.tl.draw_graph(adata, init_pos="paga", random_state=36)

# dpt pseudotime
adata.uns["iroot"] = np.flatnonzero(adata.obs["hypoxia_label"] == "lowest_hypoxia_m2")[0] # define root cluster (where to start pseudotime)
#adata.uns["iroot"] = np.flatnonzero(adata.obs["leiden"] == "2")[0] # define root cluster (where to start pseudotime)
#adata.uns["iroot"] = np.flatnonzero(adata.obs["hypoxia_classification"] == "normoxia")[0] # define root cluster (where to start pseudotime)

sc.tl.dpt(adata)

sc.pl.draw_graph(adata, color=["dpt_pseudotime", "hypoxia_classification", "preliminary_cell_types", "leiden", "author", "orig.ident", "Phase", "filttered_data_set_annotation", "seurat_clusters"], vmax=1, legend_loc="on data", save="_pseudotime_vmax1_graph_set0_1510_seed36_check_no45_nocluster13_MICROS.pdf")
sc.pl.draw_graph(adata, color=["dpt_pseudotime", "hypoxia_classification", "preliminary_cell_types", "leiden", "author", "orig.ident", "Phase", "filttered_data_set_annotation", "seurat_clusters"], vmax=0.5, legend_loc="on data", save="_pseudotime_vmax0.5_graph_set0_1510_seed36_check_no45_nocluster13_MICROS.pdf")
#sc.pl.draw_graph(adata, color=["hypoxia_classification", "leiden", "preliminary_cell_types", "dpt_pseudotime", "CD68", "CD163", "MRC1", "TMEM119", "P2RY12", "SALL1", "ITGAX", "ITGAM", "CD4"], vmax=1,legend_loc="on data", save="_pseudotime_graph_hypoxia_markers_0702_Nocluster11_seed36_vmax1.pdf")

# export
pd.DataFrame(adata.obs[['dpt_pseudotime']]).to_csv("/scratch/svc_td_cri/projects/glioma_scs_st/Narvi_glioma_scRNAseq/data/Seurat_objects/hypoxia/pseudotime/dpt_pseudotimes/MICROS_dpt_1510_check_no45_nocluster13.csv")