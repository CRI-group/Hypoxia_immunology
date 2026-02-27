library(infercnv)
library(Matrix)

df <- read.csv('/path/to/Neftel_filttered_counts.csv', row.names = 1, check.names = FALSE)

if ("X" %in% colnames(df)) {
  df <- df[, colnames(df) != "X"]
}

celltypes <- read.table('/path/to/Neftel_filttered_celltypes.csv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
rownames(celltypes) <- celltypes$V1
celltypes <- celltypes[, -1, drop = FALSE]

colnames(df) <- rownames(celltypes)

genepos <- read.table('/path/to/gene_postion.txt', header = FALSE, stringsAsFactors = FALSE)
rownames(genepos) <- genepos$V1
genepos <- genepos[, 2:4]

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = as.matrix(df),
  annotations_file = celltypes,
  gene_order_file = genepos,
  ref_group_names = NULL
)

# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = '/path/to/infercnv_results/Neftel_filttered_k10',
  denoise = TRUE,
  HMM = TRUE,
  k_obs_groups = 10, #Neftel dataset k=10, Antunes and Ricahrds k=5
  window_length = 301
)