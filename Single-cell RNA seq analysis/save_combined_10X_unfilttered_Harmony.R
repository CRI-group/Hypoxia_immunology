## This script merges the datasets and performs Harmony integration 

## Combine datasets
library(Seurat)
library(Matrix)
library(harmony)

##Read GSE140819 in
new_idents_table <- read.table('/lustre/compbio/projects/glioma_scRNAseq/Infercnv_results/GSE140819_filttered_NULL_k5/infercnv.observation_groupings.txt')
GSE140819 <- read.csv('/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/filttered_count_matrixes/GSE140819_filttered0606.csv')
rownames(GSE140819) <- GSE140819[,1]
GSE140819 = subset(GSE140819, select = -c(X))
listnames <- gsub('\\.', '-', colnames(GSE140819))
colnames(GSE140819) <- listnames

new_idents_table <- new_idents_table[,1:2]
temp <- new_idents_table[,1]
temp[temp == 2] <- 'no CNV'
temp[temp == 4 | temp==1 | temp==3 | temp==5] <- 'CNV'
new_idents_table[,2] <- temp
colnames(new_idents_table) <- c('infercnv.cluster','infercnv.annotation')
new_idents_table$dataset.name <- c('Slyper_GSE140819')

seurat_object_140 <- CreateSeuratObject(GSE140819, project = 'GSE140819_MGH125')
seurat_object_140 <- AddMetaData(seurat_object_140, new_idents_table)
rm(GSE140819)

##Read GSE163120 in
df1 <- read.csv('/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/GSE163120/GSM4972211_Human.GBM.ND1_2_3_4_5_6_7.filtered.gene.bc.matrix.csv')
df2 <- read.csv('/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/GSE163120/GSM4972210_Human.GBM.R1_2_3_4_4nc.filtered.gene.bc.matrix.csv')

rownames(df1) <- df1$X
rownames(df2) <- df2$X

df1 <- df1[,-1]
df2 <- df2[,-1]

seurat_object1_163 <- CreateSeuratObject(df1, project = 'GSE163120_ND', min.cells = 3, min.features = 200) 
samples <- read.csv('/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/GSE163120/GSM4972211_annot.Human.GBM.ND1_2_3_4_5_6_7.csv')
samples <- samples[,5:6]
rownames(samples) <- samples$cell; samples = subset(samples, select = -c(cell))
samples$orig.ident <- paste('GSE163120_', samples$sample, sep = '')
samples <- subset(samples, select = -c(sample))
rnames <- gsub('\\-', '.', rownames(samples))
rownames(samples) <- rnames
samples$infercnv.annotation <- c('no CNV')
samples$dataset.name <- c('Antunes_GSE163120_ND')

seurat_object1_163 <- AddMetaData(seurat_object1_163, samples)

seurat_object2_163 <- CreateSeuratObject(df2, project = 'GSE163120_R', min.cells = 3, min.features = 200)
samples <- read.csv('/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/GSE163120/GSM4972210_annot.Human.GBM.R1_2_3_4_4nc.csv')
samples <- samples[,5:6]
rownames(samples) <- samples$cell; samples = subset(samples, select = -c(cell))
samples$orig.ident <- paste('GSE163120_', samples$sample, sep = '')
samples <- subset(samples, select = -c(sample))
rnames <- gsub('\\-', '.', rownames(samples))
rownames(samples) <- rnames
samples$infercnv.annotation <- c('no CNV')
samples$dataset.name <- c('Antunes_GSE163120_R')

seurat_object2_163 <- AddMetaData(seurat_object2_163, samples) 
rm(df1)
rm(df2)
rm(samples)
rm(rnames)

## Merge GSE140819 and GSE163120 and save table


#REad GSE131928
df <- ReadMtx('/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/GSE131928/IDHwtGBM.processed.10X.counts.mtx', '/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/GSE131928/cells1.proc.tsv', '/lustre/compbio/projects/glioma_scRNAseq/data/opened_datasets/GSE131928/genes1.tsv')
new_idents_table <- read.table('/lustre/compbio/projects/glioma_scRNAseq/Infercnv_results/GSE131928_1_filttered_NULL_k10_0607/infercnv.observation_groupings.txt')
new_idents_table <- new_idents_table[,1:2]
temp <- new_idents_table[,1]
temp[temp==4 | temp==3] <- 'no CNV'
temp[temp==7 | temp == 9 | temp==2 | temp==1 | temp==5 | temp==10 | temp==6 | temp==8] <- 'CNV'
new_idents_table[,2] <- temp
colnames(new_idents_table) <- c('infercnv.cluster','infercnv.annotation')
new_idents_table$dataset.name <- c('Neftel_GSE131928')

seurat_object_131 <- CreateSeuratObject(df, project = 'GSE131928_1', min.cells = 3, min.features = 200) 
seurat_object_131 <- AddMetaData(seurat_object_131, new_idents_table)
rm(df)

## Read Richards scRNA, not GSC scRNA
new_idents_table <- read.table('/lustre/compbio/projects/glioma_scRNAseq/Infercnv_results/Richards_scRNA_filttered_NULL_k5/infercnv.observation_groupings.txt')
new_idents_table <- new_idents_table[,1:2]
df <- read.csv('/lustre/compbio/projects/glioma_scRNAseq/data/packaged_datasets/Richards_2021/Richards_NatureCancer_GBM_scRNAseq_counts.csv.gz')
rownames(df) <- df[,1]
df = subset(df, select = -c(X) )

temp <- new_idents_table[,1]
temp[temp == 4 | temp== 2] <- 'no CNV'
temp[temp == 3 | temp==1 | temp==5] <- 'CNV'
new_idents_table[,2] <- temp
colnames(new_idents_table) <- c('infercnv.cluster','infercnv.annotation')
new_idents_table$dataset.name <- c('Richards')

seurat_object_Richards <- CreateSeuratObject(df, project = 'Richards', min.cells = 3, min.features = 200) 
seurat_object_Richards <- AddMetaData(seurat_object_Richards, new_idents_table)
rm(df)

print('all done, harmony next')
#Merge Richards, 131, 140 and 163
# normalize and identify variable features for each dataset independently
seurat_object_harmony <- merge(seurat_object_140, y = c(seurat_object1_163, seurat_object2_163, seurat_object_131, seurat_object_Richards), add.cell.ids = c("GSE140", "GSE163_1", "GSE163_2", "GSE131", 'Richards'), project = "Harmony_scRNA_set")
print('merge done')

seurat_object_harmony[["percent.mt"]] <- PercentageFeatureSet(seurat_object_harmony, pattern = "^MT-") 
seurat_object_harmony <- subset(seurat_object_harmony, subset = percent.mt < 5) 

seurat_object_harmony <- NormalizeData(seurat_object_harmony, verbose = F)
seurat_object_harmony <- FindVariableFeatures(seurat_object_harmony, selection.method = "vst", nfeatures = 2000, verbose = F)
seurat_object_harmony <- ScaleData(seurat_object_harmony, verbose = F)
seurat_object_harmony <- RunPCA(seurat_object_harmony, npcs = 50, verbose = F)
seurat_object_harmony <- RunHarmony(seurat_object_harmony, group.by.vars = "orig.ident")

saveRDS(seurat_object_harmony, "/lustre/compbio/projects/glioma_scRNAseq/data/combined_10X_unfilttered_CNV_final_mt_removed_Harmony.RDS")










