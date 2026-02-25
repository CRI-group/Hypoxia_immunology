library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(plyr)
library(ggplot2)
library(rstatix)
library(gridExtra)
library(ggpubr)
library(gplots)
library(IlluminaHumanMethylation450kprobe)
library(IlluminaHumanMethylation450k.db)
library(TCGAbiolinks)
library(RColorBrewer)

########################################
# download TCGA data
query_legacy_gbm_gene_expression_illumina = GDCquery(project="TCGA-GBM",
                                                    data.category="Gene expression",
                                                    data.type="Gene expression quantification",
                                                    platform="Illumina HiSeq",
                                                    file.type="normalized_results",
                                                    legacy=T)
GDCdownload(query_legacy_gbm_gene_expression_illumina)
gbm_gene_hg19_illumina = GDCprepare(query_legacy_gbm_gene_expression_illumina)

save(gbm_gene_hg19_illumina, file="TCGA/RNA-seq/gbm_expr_hg19_all_illumina_normalized.RData")

query_legacy_lgg_gene_expression_illumina = GDCquery(project="TCGA-LGG",
                                                     data.category = "Gene expression",
                                                     data.type = "Gene expression quantification",
                                                     platform="Illumina HiSeq",
                                                     file.type = "normalized_results",
                                                     legacy=T)
GDCdownload(query_legacy_lgg_gene_expression_illumina)
lgg_gene_hg19_illumina = GDCprepare(query_legacy_lgg_gene_expression_illumina)

save(lgg_gene_hg19_illumina, file="TCGA/RNA-seq/lgg_expr_hg19_all_illumina_normalized.RData")

gbm_matrix = SummarizedExperiment::assay(gbm_gene_hg19_illumina, "normalized_count")
lgg_matrix = SummarizedExperiment::assay(lgg_gene_hg19_illumina, "normalized_count")
expr_matrix = cbind(lgg_matrix, gbm_matrix)
expr_log2 = log2(expr_matrix + 1)

# methylation data
query_legacy_gbm_meth_illumina450 = GDCquery(project="TCGA-GBM",
                                             data.category="DNA methylation",
                                             platform="Illumina Human Methylation 450",
                                             legacy=T)
GDCdownload(query_legacy_gbm_meth_illumina450)
gbm_meth_hg19_illumina450 = GDCprepare(query_legacy_gbm_meth_illumina450)

save(gbm_meth_hg19_illumina450, file="../TCGA/Methylation/gbm_meth_hg19_illumina450.RData")

query_legacy_lgg_meth_illumina450 = GDCquery(project="TCGA-LGG",
                                             data.category="DNA methylation",
                                             platform="Illumina Human Methylation 450",
                                             legacy=T)
GDCdownload(query_legacy_lgg_meth_illumina450, files.per.chunk = 10)
lgg_meth_hg19_illumina450 = GDCprepare(query_legacy_lgg_meth_illumina450)

save(lgg_meth_hg19_illumina450, file="lgg_meth_hg19_illumina450.RData")

gbm_beta = as.data.frame(SummarizedExperiment::assay(gbm_meth_hg19_illumina450))
gbm_beta = cbind(SummarizedExperiment::rowData(gbm_meth_hg19_illumina450), gbm_beta)
gbm_beta_coldata = SummarizedExperiment::colData(gbm_meth_hg19_illumina450)

lgg_beta = as.data.frame(SummarizedExperiment::assay(lgg_meth_hg19_illumina450))
lgg_beta = cbind(SummarizedExperiment::rowData(lgg_meth_hg19_illumina450), lgg_beta)
lgg_hg19_coldata = SummarizedExperiment::colData(lgg_meth_hg19_illumina450)


####################################
# create annotation dataframe

tumor_classes = as.data.frame(colnames(expr_matrix))
colnames(tumor_classes) = "sample"

tumor_classes[,"idh_codel_subtype"] = NA
tumor_classes[,"grade"] = NA
tumor_classes[,"survival"] = NA
tumor_classes[,"meth_cluster"] = NA
tumor_classes[,"chr7_chr10"] = NA
tumor_classes[,"TERT_promoter"] = NA
tumor_classes[,"EGFR_cn"] = NA
tumor_classes[,"CDKN2A_cn"] = NA
tumor_classes[,"CDKN2A_expr"] = NA

# CNAs for selected genes downloaded manually from cbioportal.org
cnvs = read.table("TCGA/all_cna.txt", sep="\t", header=T, stringsAsFactors = F)

for (i in 1:nrow(tumor_classes)) {
  sample = tumor_classes[i,"sample"]
  tumor_classes[i, "CDKN2A_expr"] = expr_matrix["CDKN2A", sample]

  sample = substr(sample, 1, 12)  
  if (any(grep(sample, sampleinfo$Case))) {
    tumor_classes[i, "idh_codel_subtype"] = sampleinfo[sampleinfo$Case==sample, "IDH/codel.subtype"]
    tumor_classes[i, "survival"] = sampleinfo[sampleinfo$Case==sample, "Survival.(months)"]
    tumor_classes[i, "grade"] = sampleinfo[sampleinfo$Case==sample, "Grade"]
    tumor_classes[i, "meth_cluster"] = sampleinfo[sampleinfo$Case==sample, "Supervised.DNA.Methylation.Cluster"]
    tumor_classes[i, "chr7_chr10"] = sampleinfo[sampleinfo$Case==sample, "Chr.7.gain/Chr.10.loss"]
    tumor_classes[i, "TERT_promoter"] = sampleinfo[sampleinfo$Case==sample, "TERT.promoter.status"]
  } 
  if (any(grep(sample, cnvs$SAMPLE_ID))) {
    tumor_classes[i, "CDKN2A_cn"] = as.numeric(cnvs[grep(sample, cnvs$SAMPLE_ID), "CDKN2A"])  
    tumor_classes[i, "EGFR_cn"] = as.numeric(cnvs[grep(sample, cnvs$SAMPLE_ID), "EGFR"])
  }
}

tumor_classes$CDKN2A_expr_log2 = log2(tumor_classes$CDKN2A_expr)
tumor_classes$survival = as.numeric(tumor_classes$survival)
tumor_classes$IDH_codel_grade = paste0(tumor_classes$idh_codel_subtype, "_", tumor_classes$grade)
#708 cases
tumor_classes = tumor_classes[complete.cases(tumor_classes[,c("sample", "idh_codel_subtype", "grade")]),]

# update classification according to WHO 2021
keep_rows = tumor_classes$TERT_promoter=="Mutant" | tumor_classes$EGFR_cn > 0 | tumor_classes$chr7_chr10=="Gain chr 7 & loss chr 10"
keep_rows[is.na(keep_rows)] = FALSE
tumor_classes_IDHwt_pass = tumor_classes[(tumor_classes$idh_codel_subtype!="IDHwt" | tumor_classes$grade=="G4" | keep_rows), ]

tumor_classes_IDHwt_pass$grade_2021 = tumor_classes_IDHwt_pass$grade
tumor_classes_IDHwt_pass[tumor_classes_IDHwt_pass$IDH_codel_grade=="IDHwt_G2", ]$grade_2021 = "G4"
tumor_classes_IDHwt_pass[tumor_classes_IDHwt_pass$IDH_codel_grade=="IDHwt_G3", ]$grade_2021 = "G4"

astro_G4_rows = (tumor_classes_IDHwt_pass$IDH_codel_grade=="IDHmut-non-codel_G2" | tumor_classes_IDHwt_pass$IDH_codel_grade=="IDHmut-non-codel_G3") & tumor_classes_IDHwt_pass$CDKN2A_cn==-2
astro_G4_rows[is.na(astro_G4_rows)] = FALSE
tumor_classes_IDHwt_pass[astro_G4_rows, ]$grade_2021 = "G4"
tumor_classes_IDHwt_pass$IDH_codel_grade_2021 = paste0(tumor_classes_IDHwt_pass$idh_codel_subtype, "_", tumor_classes_IDHwt_pass$grade_2021)

colnames(tumor_classes_IDHwt_pass)[colnames(tumor_classes_IDHwt_pass)=="grade"] = "grade_old"
colnames(tumor_classes_IDHwt_pass)[colnames(tumor_classes_IDHwt_pass)=="IDH_codel_grade"] = "IDH_codel_grade_old"
colnames(tumor_classes_IDHwt_pass)[colnames(tumor_classes_IDHwt_pass)=="grade_2021"] = "grade"
colnames(tumor_classes_IDHwt_pass)[colnames(tumor_classes_IDHwt_pass)=="IDH_codel_grade_2021"] = "IDH_codel_grade"

# 602 cases
tumor_classes = tumor_classes_IDHwt_pass

# keep only cases with expression data
tumor_classes = tumor_classes[complete.cases(tumor_classes[,c("sample", "idh_codel_subtype", "grade", "CDKN2A_expr")]),]

# keep those recurrents where primary grade is 4, those which only have recurrent sample
# keep TCGA-DU-6404-02A and -02B, DU-6407-02A and -02B, FG-5963, FG-5965, 06-0125, 06-0190, 06-0210, 06-0211, 14-1034
# keep TCGA-06-0152, 06-0171, 06-0221, 19-0957, 19-1389
# remove DU-5870-02A, DU-5872-02A, FG-A4MT-02A, TM-A7CF-02A, DH-A669-02A, DU-6397-02A, DU-7304-02A
tumor_classes = tumor_classes[-grep("DU-5870-02A|DU-5872-02A|FG-A4MT-02A|TM-A7CF-02A|DH-A669-02A|DU-6397-02A|DU-7304-02A", tumor_classes$sample),]

# previous radiation
# yes: 06-0152, 06-0171, 06-0221, 19-0957, 19-1389, 06-0125, 14-1034, DU-6404, DU-6407, FG-5963
# no: 06-0210, 06-0211
# NA: FG-5965
tumor_classes$previous_radiation = NA
tumor_classes[grepl("06-0152|06-0171|06-0221|19-0957|19-1389|06-0125|14-1034|DU-6404|DU-6407|FG-5963", tumor_classes$sample) & 
                grepl("02A|02B", tumor_classes$sample), ]$previous_radiation = "yes"
tumor_classes[grepl("06-0210|06-0211", tumor_classes$sample) & grepl("02A|02B", tumor_classes$sample), ]$previous_radiation = "no"

# 595 cases
tumor_classes = tumor_classes

tumor_classes$primary_recurrent = "primary"
tumor_classes[grep("02A|02B", tumor_classes$sample), ]$primary_recurrent = "recurrent"
tumor_classes$idh_codel_primary_recurrent = paste0(tumor_classes$idh_codel_subtype, "_", tumor_classes$primary_recurrent)


##################################################
# Methylation probe selection

# load probe annotations
data(IlluminaHumanMethylation450kprobe)
IlluminaHumanMethylation450kprobe = IlluminaHumanMethylation450kprobe[,1:5]

beta_CA9 = cbind(as.data.frame(as.matrix(gbm_beta[which(gbm_beta$Gene_Symbol=="CA9"), ]), stringsAsFactors = F),
                    as.data.frame(as.matrix(lgg_beta[which(lgg_beta$Gene_Symbol=="CA9"), 3:536]), stringsAsFactors=F))
beta_CA9[, 3:ncol(beta_CA9)] = lapply(beta_CA9[,3:ncol(beta_CA9)], as.numeric)
beta_CA9 = beta_CA9[rowSums(is.na(beta_CA9[3:ncol(beta_CA9)])) != ncol(beta_CA9)-2, ]

beta_ADM = cbind(as.data.frame(as.matrix(gbm_beta[which(gbm_beta$Gene_Symbol=="ADM"), ]), stringsAsFactors = F),
                   as.data.frame(as.matrix(lgg_beta[which(lgg_beta$Gene_Symbol=="ADM"), 3:536]), stringsAsFactors=F))
beta_ADM[, 3:ncol(beta_ADM)] = lapply(beta_ADM[,3:ncol(beta_ADM)], as.numeric)
beta_ADM = beta_ADM[rowSums(is.na(beta_ADM[3:ncol(beta_ADM)])) != ncol(beta_ADM)-2, ]

# expression of the cases with methylation values
expr_CA9 = matrix(nrow=ncol(beta_CA9)-2, ncol=2)
rownames(expr_CA9) = colnames(beta_CA9)[3:ncol(beta_CA9)]
colnames(expr_CA9) = c("expression", "expression_log2")
for (i in 1:nrow(expr_CA9)) {
  sample = rownames(expr_CA9)[i]
  sample = substr(sample, 1, 16)
  expr_CA9[i,"expression"] = median(unlist(expr_matrix["CA9", grep(sample, colnames(expr_matrix))]))
  expr_CA9[i,"expression_log2"] = median(unlist(expr_log2["CA9", grep(sample, colnames(expr_matrix))]))
}

expr_ADM = matrix(nrow=ncol(beta_ADM)-2, ncol=2)
rownames(expr_ADM) = colnames(beta_ADM)[3:ncol(beta_ADM)]
colnames(expr_ADM) = c("expression", "expression_log2")
for (i in 1:nrow(expr_ADM)) {
  sample = rownames(expr_ADM)[i]
  sample = substr(sample, 1, 16)
  expr_ADM[i,"expression"] = median(unlist(expr_matrix["ADM", grep(sample, colnames(expr_matrix))]))
  expr_ADM[i,"expression_log2"] = median(unlist(expr_log2["ADM", grep(sample, colnames(expr_matrix))]))
}

median_CA9 = apply(beta_CA9[, 3:ncol(beta_CA9)], 2, median, na.rm=T)
median_CA9 = data.frame(median_CA9)
median_ADM = apply(beta_ADM[, 3:ncol(beta_ADM)], 2, median, na.rm=T)
median_ADM = data.frame(median_ADM)

CA9_tss = 35673925
ADM_tss = 10326620

expr_ADM = as.data.frame(expr_ADM)
expr_CA9 = as.data.frame(expr_CA9)
expr_ADM$subtype = NA
expr_CA9$subtype = NA

for (i in 1:nrow(expr_ADM)) {
  sample = rownames(expr_ADM)[i]
  sample = substr(sample,1 ,16)
  subtype = tumor_classes[grep(sample, tumor_classes$sample), "idh_codel_subtype"]
  if (length(subtype)==1) {
  expr_ADM$subtype[i] = subtype
  expr_CA9$subtype[i] = subtype
  expr_VEGFA$subtype[i] = subtype
  expr_PDK1$subtype[i] = subtype

  } else if (length(subtype)>1) {
    print(i)
  }
}

plot_heatmap = function(df, name, tss, median_data, expression_data) {
  
  max_dist = 10000
  min_var = 0.02
  
  df = cbind(expression_data, t(df[,3:ncol(df)]))
  df = df[complete.cases(df[,c(1:2, 4:ncol(df))]),]
  
  median_data = as.data.frame(median_data[rownames(df), 1:ncol(median_data), drop=F])
  colnames(median_data)[1] = "all_median"
  
  probe_annot = IlluminaHumanMethylation450kprobe[IlluminaHumanMethylation450kprobe$Probe_ID %in% colnames(df[,4:ncol(df)]),]
  probe_annot$dist_TSS = apply(abs(cbind(tss-probe_annot$start, tss-probe_annot$end)), 1, FUN=min)
  
  #print(probe_annot)
  
  keep_probes = probe_annot[probe_annot$dist_TSS<max_dist, "Probe_ID"]
  df_close = df[, c("expression", "expression_log2", "subtype", keep_probes)]
  
  if (length(keep_probes) > 0) { 
    
    var_filter = apply( df_close[,4:ncol(df_close), drop=F], MARGIN=2, FUN=var, na.rm=T)
    #print(var_filter)
    df_var = df_close[, names(which(var_filter > min_var)), drop=F]
    #print(df[1:10])
    
    median_closest = apply(df_var, 1, median)
    median_write = data.frame(as.list(median_closest), stringsAsFactors = F, check.names=F)
    # use 10000 dist and 0.02 variance
    write.table(t(median_write), paste0("TCGA/HM450/", name, "_methylation_median_filtered.txt"), sep="\t", quote=F, row.names = T)
    
    beta_col = colorRamp2(c(0,0.15,0.5, 0.65,1), c("royalblue4", "royalblue1", "white","salmon", "red"))
    expr_col = colorRamp2(c(0,7,14), c("forestgreen", "lightyellow", "orangered"))
    subtype_col = c("IDHmut-codel"="forestgreen", "IDHmut-non-codel"="purple", "IDHwt"="orange")
    
    probe_annot_var = IlluminaHumanMethylation450kprobe[IlluminaHumanMethylation450kprobe$Probe_ID %in% colnames(df_var),]
    probe_annot_var$dist_TSS = apply(abs(cbind(tss-probe_annot_var$start, tss-probe_annot_var$end)), 1, FUN=min)
    
    if (ncol(df_close)>3) {
      cluster_columns=T
    } else {
      cluster_columns=F
    }
    
    if(ncol(df_var)>1) {
      cluster_columns_var=T
    } else {
      cluster_columns_var=F
    }
    
    if (ncol(median_data)==1) {
      ha_col = list(project=c("TCGA-GBM" = "cadetblue3", "TCGA-LGG" = "cadetblue4"), 
                    subtype=c("IDHmut-codel"="darkolivegreen2", "IDHmut-non-codel"="violet", "IDHwt"="sienna"),
                    all_median=beta_col, median_closest=beta_col)
    } else if (ncol(median_data)==3) {
      ha_col = list(project=c("TCGA-GBM" = "cadetblue3", "TCGA-LGG" = "cadetblue4"), 
                    subtype=c("IDHmut-codel"="darkolivegreen2", "IDHmut-non-codel"="violet", "IDHwt"="sienna"),
                    all_median=beta_col, cluster1_median=beta_col, 
                    cluster2_median=beta_col, median_closest=beta_col)
    }
    
 
    gbms = grep("TCGA-[0-9]", rownames(df_close))
    lggs = grep("TCGA-[A-Z]", rownames(df_close))
    
    ha_top = HeatmapAnnotation(project=c(rep("TCGA-GBM", length(gbms)), rep("TCGA-LGG", length(lggs))),
                               subtype=df_close$subtype,
                               df = median_data,
                               median_closest = median_closest,
                               col = ha_col,
                               show_legend = c(T,T,rep(F, ncol(median_data)),F))
    
    dist_col = colorRamp2(c(0,max(probe_annot$dist_TSS)), c("white", "slateblue4"))
    
    ha_row = rowAnnotation(dist_TSS = probe_annot[keep_probes, ]$dist_TSS, show_legend=T, col=list(dist_TSS=dist_col))
    ha_row_var = rowAnnotation(dist_TSS=probe_annot_var$dist_TSS, show_legend=T, col=list(dist_TSS=dist_col))
    
    hm_expr = Heatmap(as.matrix(t(df_close[,1])),
                      cluster_rows=F, show_column_names=F,
                      heatmap_legend_param=list(title="Expression"))
    
    hm_expr_log2 = Heatmap(as.matrix(t(df_close[,2])),
                           cluster_rows=F,  show_column_names=F,
                           heatmap_legend_param=list(title="Expression log2"))
    
    hm = Heatmap(as.matrix(t(df_close[,4:ncol(df_close)])), top_annotation = ha_top, right_annotation=ha_row, 
                 cluster_rows=T, cluster_columns = cluster_columns, 
                 show_column_names = F, col=beta_col,
                 clustering_distance_rows = "spearman", clustering_distance_columns = "spearman",
                 name=paste(name, "beta"), column_title=paste(name, "all probes pearson clustering"))
    
    hm_list = hm %v% hm_expr %v% hm_expr_log2
    draw(hm_list)
    
    hm_var = Heatmap(as.matrix(t(df_var)), top_annotation = ha_top, right_annotation=ha_row_var,
                     cluster_rows=T, cluster_columns = cluster_columns_var, 
                     show_column_names = F, col=beta_col,
                     clustering_distance_rows = "spearman", clustering_distance_columns = "spearman",
                     name=paste(name, "beta"), column_title=paste(name, "variance filtered", min_var,"probes spearman clustering"))
    
    hm_var_list = hm_var %v% hm_expr  %v% hm_expr_log2
    draw(hm_var_list)
    
    #dev.off()
    
    pdf(paste0("TCGA/HM450/", name, "_methylation_expression_scatter_max_dist_", max_dist, "_min_var_", min_var,"_subtypes.pdf"), height=4, width=4)
    
    
    df_scatter = df[complete.cases(df),]
    
    cp_all_idhmutcodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-codel",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHmut-codel",]),1], method="pearson")
    cs_all_idhmutcodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-codel",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHmut-codel",]),1], method="spearman", exact=F)
    
    cp_closest_idhmutcodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-codel",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHmut-codel",])], method="pearson")
    cs_closest_idhmutcodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-codel",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHmut-codel",])], method="spearman", exact=F)
    
    cp_all_idhmutnoncodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-non-codel",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHmut-non-codel",]),1], method="pearson")
    cs_all_idhmutnoncodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-non-codel",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHmut-non-codel",]),1], method="spearman", exact=F)
    
    cp_closest_idhmutnoncodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-non-codel",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHmut-non-codel",])], method="pearson")
    cs_closest_idhmutnoncodel = cor.test(df_scatter[df_scatter$subtype=="IDHmut-non-codel",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHmut-non-codel",])], method="spearman", exact=F)
    
    cp_all_idhwt = cor.test(df_scatter[df_scatter$subtype=="IDHwt",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHwt",]),1], method="pearson")
    cs_all_idhwt = cor.test(df_scatter[df_scatter$subtype=="IDHwt",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHwt",]),1], method="spearman", exact=F)
    
    cp_closest_idhwt = cor.test(df_scatter[df_scatter$subtype=="IDHwt",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHwt",])], method="pearson")
    cs_closest_idhwt = cor.test(df_scatter[df_scatter$subtype=="IDHwt",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHwt",])], method="spearman", exact=F)
    
    plot(df_scatter[df_scatter$subtype=="IDHmut-codel",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHmut-codel",]),1], 
         ylab="beta", xlab="expression", 
         main=paste0("IDHmut-codel - all median\n pearson: ", 
                     format(round(cp_all_idhmutcodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_all_idhmutcodel$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_all_idhmutcodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_all_idhmutcodel$p.value, 5), nsmall=5, scientific=T)))
    
    plot(df_scatter[df_scatter$subtype=="IDHmut-codel",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHmut-codel",])], 
         ylab="beta", xlab="expression",  
         main=paste0("IDHmut-codel - median closest\n pearson: ", 
                     format(round(cp_closest_idhmutcodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_closest_idhmutcodel$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_closest_idhmutcodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_closest_idhmutcodel$p.value, 5), nsmall=5, scientific=T)))
    
    plot(df_scatter[df_scatter$subtype=="IDHmut-non-codel",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHmut-non-codel",]),1], 
         ylab="beta", xlab="expression", 
         main=paste0("IDHmut-non-codel - all median\n pearson: ", 
                     format(round(cp_all_idhmutnoncodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_all_idhmutnoncodel$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_all_idhmutnoncodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_all_idhmutnoncodel$p.value, 5), nsmall=5, scientific=T)))
    
    plot(df_scatter[df_scatter$subtype=="IDHmut-non-codel",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHmut-non-codel",])], 
         ylab="beta", xlab="expression",  
         main=paste0("IDHmut-non-codel - median closest\n pearson: ", 
                     format(round(cp_closest_idhmutnoncodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_closest_idhmutnoncodel$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_closest_idhmutnoncodel$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_closest_idhmutnoncodel$p.value, 5), nsmall=5, scientific=T)))
    
    plot(df_scatter[df_scatter$subtype=="IDHwt",2], median_data[rownames(df_scatter[df_scatter$subtype=="IDHwt",]),1], 
         ylab="beta", xlab="expression", 
         main=paste0("IDHwt - all median\n pearson: ", 
                     format(round(cp_all_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_all_idhwt$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_all_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_all_idhwt$p.value, 5), nsmall=5, scientific=T)))
    
    plot(df_scatter[df_scatter$subtype=="IDHwt",2], median_closest[rownames(df_scatter[df_scatter$subtype=="IDHwt",])], 
         ylab="beta", xlab="expression",  
         main=paste0("IDHwt - median closest\n pearson: ", 
                     format(round(cp_closest_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_closest_idhwt$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_closest_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_closest_idhwt$p.value, 5), nsmall=5, scientific=T)))
    
    #dev.off()
        
    plot(df_scatter[df_scatter$subtype!="IDHmut-codel",2], median_data[rownames(df_scatter[df_scatter$subtype!="IDHmut-codel",]),1], 
         ylab="beta", xlab="expression", 
         main=paste0("GBM+astro - all median\n pearson: ", 
                     format(round(cp_all_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_all_idhwt$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_all_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_all_idhwt$p.value, 5), nsmall=5, scientific=T)))
    
    plot(df_scatter[df_scatter$subtype!="IDHmut-codel",2], median_closest[rownames(df_scatter[df_scatter$subtype!="IDHmut-codel",])], 
         ylab="beta", xlab="expression",
         main=paste0("GBM+astro - median closest\n pearson: ", 
                     format(round(cp_closest_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cp_closest_idhwt$p.value, 5), nsmall=5, scientific=T), "\nspearman: ",
                     format(round(cs_closest_idhwt$estimate, 5), nsmall=5), " p: ", 
                     format(round(cs_closest_idhwt$p.value, 5), nsmall=5, scientific=T)))
    
    dev.off()
    
  } else {
    print(probe_annot)
  }
}


plot_heatmap(beta_CA9, "CA9", CA9_tss, median_CA9, expr_CA9)
plot_heatmap(beta_ADM, "ADM", ADM_tss, median_ADM, expr_ADM)

#######################################
# add hypoxia genes to annotation dataframe

tumor_classes = tumor_classes[c("sample", "idh_codel_subtype", "grade_old", "survival", "meth_cluster", 
                                "chr7_chr10", "TERT_promoter", "EGFR_cn",
                                "grade", "IDH_codel_grade", "previous_radiation", "primary_recurrent")]

sampleinfo = read.xlsx("TCGA/GBM_LGG_sampleinfo.xlsx", sheet=1, startRow=2)

# hypoxia correlating signature defined below
hypoxia_signature = read.table("TCGA/RNA-seq/GSVA/hypoxia_correlating_genes_TCGA.txt", sep="\t", header=T)
hypoxia_expr = expr_log2[unlist(hypoxia_signature), tumor_classes$sample]

calculate_score = function(df) {
  z_score = t(scale(t(df)))
  activity.z = colSums(z_score/sqrt(nrow(df)))
  return(activity.z)
}

sc_hypoxia_all = calculate_score(hypoxia_expr)
sc_hypoxia_astro = calculate_score(hypoxia_expr[, tumor_classes[tumor_classes$idh_codel_subtype=="IDHmut-non-codel", "sample"]])
sc_hypoxia_astro_G4 = calculate_score(hypoxia_expr[, tumor_classes[tumor_classes$IDH_codel_grade=="IDHmut-non-codel_G4" , "sample"]])
sc_hypoxia_astro_primary = calculate_score(hypoxia_expr[, tumor_classes[tumor_classes$idh_codel_subtype=="IDHmut-non-codel" & tumor_classes$primary_recurrent=="primary", "sample"]])
sc_hypoxia_astro_recurrent = calculate_score(hypoxia_expr[, tumor_classes[tumor_classes$idh_codel_subtype=="IDHmut-non-codel" & tumor_classes$primary_recurrent=="recurrent", "sample"]])


tumor_classes$vital_status = NA
tumor_classes$age = NA

tumor_classes$TMEM119_expr = NA
tumor_classes$TMEM119_expr_log2 = NA
tumor_classes$CA9_expr = NA
tumor_classes$CA9_expr_log2 = NA
tumor_classes$VEGFA_expr = NA
tumor_classes$VEGFA_expr_log2 = NA
tumor_classes$ADM_expr = NA
tumor_classes$ADM_expr_log2 = NA
tumor_classes$PDK1_expr = NA
tumor_classes$PDK1_expr_log2 = NA

tumor_classes$CA9_beta_cg06908460 = NA
tumor_classes$CA9_beta_cg09566069 = NA
tumor_classes$CA9_beta_cg13849253 = NA
tumor_classes$CA9_beta_cg13938361 = NA
tumor_classes$CA9_beta_cg14563831 = NA
tumor_classes$CA9_beta_cg19257550 = NA
tumor_classes$CA9_beta_cg20610181 = NA

tumor_classes$ADM_beta_cg08259810 = NA
tumor_classes$ADM_beta_cg05149386 = NA
tumor_classes$ADM_beta_cg19315081 = NA
tumor_classes$ADM_beta_cg20741987 = NA
tumor_classes$ADM_beta_cg03548673 = NA
tumor_classes$ADM_beta_cg02843237 = NA
tumor_classes$ADM_beta_cg10044466 = NA
tumor_classes$ADM_beta_cg06875754 = NA
tumor_classes$ADM_beta_cg05727225 = NA

tumor_classes$hypoxia_score_all = NA

for (i in 1:nrow(tumor_classes)) {
  sample = tumor_classes$sample[i]
  case = substr(sample, 1, 12)
  tumor_classes$vital_status[i] = sampleinfo[sampleinfo$Case==case, "Vital.status.(1=dead)"]
  tumor_classes$age[i] = sampleinfo[sampleinfo$Case==case, "Age.(years.at.diagnosis)"]
  
  tumor_classes$TMEM119_expr[i] = expr_matrix["TMEM119", sample]
  tumor_classes$TMEM119_expr_log2[i] = expr_log2["TMEM119", sample]
  tumor_classes$CA9_expr[i] = expr_matrix["CA9", sample]
  tumor_classes$CA9_expr_log2[i] = expr_log2["CA9", sample]
  tumor_classes$VEGFA_expr[i] = expr_matrix["VEGFA", sample]
  tumor_classes$VEGFA_expr_log2[i] = expr_log2["VEGFA", sample]
  tumor_classes$ADM_expr[i] = expr_matrix["ADM", sample]
  tumor_classes$ADM_expr_log2[i] = expr_log2["ADM", sample]
  tumor_classes$PDK1_expr[i] = expr_matrix["PDK1", sample]
  tumor_classes$PDK1_expr_log2[i] = expr_log2["PDK1", sample]
  
  case2 = substr(sample, 1, 16)
  
  if (length(grep(case2, colnames(beta_CA9)))==1) {
    tumor_classes$CA9_beta_cg06908460[i] = beta_CA9["cg06908460",grep(case2, colnames(beta_CA9))]
    tumor_classes$CA9_beta_cg09566069[i] = beta_CA9["cg09566069",grep(case2, colnames(beta_CA9))]
    tumor_classes$CA9_beta_cg13849253[i] = beta_CA9["cg13849253",grep(case2, colnames(beta_CA9))]
    tumor_classes$CA9_beta_cg13938361[i] = beta_CA9["cg13938361",grep(case2, colnames(beta_CA9))]
    tumor_classes$CA9_beta_cg14563831[i] = beta_CA9["cg14563831",grep(case2, colnames(beta_CA9))]
    tumor_classes$CA9_beta_cg19257550[i] = beta_CA9["cg19257550",grep(case2, colnames(beta_CA9))]
    tumor_classes$CA9_beta_cg20610181[i] = beta_CA9["cg20610181",grep(case2, colnames(beta_CA9))]
    
    tumor_classes$ADM_beta_cg08259810[i] = beta_ADM["cg08259810",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg05149386[i] = beta_ADM["cg05149386",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg19315081[i] = beta_ADM["cg19315081",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg20741987[i] = beta_ADM["cg20741987",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg03548673[i] = beta_ADM["cg03548673",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg02843237[i] = beta_ADM["cg02843237",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg10044466[i] = beta_ADM["cg10044466",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg06875754[i] = beta_ADM["cg06875754",grep(case2, colnames(beta_ADM))]
    tumor_classes$ADM_beta_cg05727225[i] = beta_ADM["cg05727225",grep(case2, colnames(beta_ADM))]
    
  } else if (length(grep(case2, colnames(beta_CA9)))>1) {
    print(c(sample))
  }
  tumor_classes$hypoxia_score_all[i] = sc_hypoxia_all[sample]
}


###########################
# violin plots

col_rad = c("yes"="blue", "no"="red", "NA"="black")
grade_col = c("G2"=rgb(62,183,99, maxColorValue = 255), "G3"=rgb(150,131,189, maxColorValue = 255), "G4"="#FDC086")
col_rec = c("recurrent"="red", "primary"="black")

plot_violin = function(data, column, datatype = c("expression", "score", "methylation"), 
                        stat.test = NULL, stat.test.grade = NULL, stat.test.subtype = NULL, 
                        gene, label = "p.adj.signif", hide.ns = TRUE, outdir = "TCGA/hypoxia/") {
  
  datatype <- match.arg(datatype)
  
  ylab_text <- switch(datatype, expression = "RNA expression", score = "Activity score", methylation = "Methylation beta")
  ylim_vals <- switch(datatype, expression = c(0, 17.5), score = NULL, methylation = c(0, 1))
  file_suffix <- switch(datatype, expression = "expr_log2", score = "score", methylation = "beta")
  
  if (!is.null(stat.test)) {
    stat.test <- stat.test %>% add_xy_position(x = "idh_codel_subtype", step.increase = 0.2)
    stat.test.grade <- stat.test.grade %>% add_xy_position(x = "idh_codel_subtype", step.increase = 0.2)
    stat.test.subtype <- add_xy_position(test = stat.test.subtype, group = "grade", x = "idh_codel_subtype")
  }

  p <- ggplot(data, aes(x = idh_codel_subtype, y = .data[[column]])) +
    geom_violin(position = position_dodge(width = 0.8), aes(fill = grade)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.25, seed = 1, dodge.width = 0.8),
                alpha = 0.5, aes(group = grade, fill = grade, col = as.factor(primary_recurrent))) +
    scale_color_manual(values = col_rec, na.value = "black", name = "Primary/recurrent") +
    scale_fill_manual(values = grade_col, name = "Grade") +
    labs(x = "Tumor subtype", y = ylab_text) +
    theme_classic() +
    guides(x = guide_axis(angle = 45)) +
    ggtitle(gene)
  
  if (!is.null(ylim_vals)) {
    p <- p + coord_cartesian(ylim = ylim_vals)
  }
  
  if (!is.null(stat.test)) {
    p_sig <- p + stat_pvalue_manual(stat.test, label = label, hide.ns = hide.ns)
    p_sig_grade <- p + stat_pvalue_manual(stat.test.grade, label = label, hide.ns = hide.ns)
    p_sig_subtype <- p + stat_pvalue_manual(stat.test.subtype, label = label, hide.ns = hide.ns)
    
    ggsave(paste0(outdir, gene, "_", file_suffix, "_violin.pdf"), p_sig, width = 5.7, height = 5)
    ggsave(paste0(outdir, gene, "_", file_suffix, "_violin_grade.pdf"), p_sig_grade, width = 5.7, height = 5)
    ggsave(paste0(outdir, gene, "_", file_suffix, "_violin_subtype.pdf"), p_sig_subtype, width = 5.7, height = 5)
  } else {
    ggsave(paste0(outdir, gene, "_", file_suffix, "_violin.pdf"), p, width = 5, height = 5)
  }
}

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(TMEM119_expr_log2 ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(TMEM119_expr_log2 ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(TMEM119_expr_log2~idh_codel_subtype)
plot_violin(tumor_classes,"TMEM119_expr_log2", "expression", stat.test, stat.test.grade, stat.test.subtype, "TMEM119", label="p.adj", hide.ns=F)

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(CA9_expr_log2 ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(CA9_expr_log2 ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(CA9_expr_log2~idh_codel_subtype)
plot_violin(tumor_classes,"CA9_expr_log2", "expression", stat.test, stat.test.grade, stat.test.subtype, "CA9", label="p.adj", hide.ns=F)

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(VEGFA_expr_log2 ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(VEGFA_expr_log2 ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(VEGFA_expr_log2~idh_codel_subtype)
plot_violin(tumor_classes,"VEGFA_expr_log2", "expression", stat.test, stat.test.grade, stat.test.subtype, "VEGFA", label="p.adj", hide.ns=F)

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(ADM_expr_log2 ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(ADM_expr_log2 ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(ADM_expr_log2~idh_codel_subtype)
plot_violin(tumor_classes,"ADM_expr_log2", "expression", stat.test, stat.test.grade, stat.test.subtype, "ADM", label="p.adj", hide.ns=F)

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(PDK1_expr_log2 ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(PDK1_expr_log2 ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(PDK1_expr_log2~idh_codel_subtype)
plot_violin(tumor_classes,"PDK1_expr_log2", "expression", stat.test, stat.test.grade, stat.test.subtype, "PDK1", label="p.adj", hide.ns=F)

# hypoxia score
stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(hypoxia_score_all ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(hypoxia_score_all ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(hypoxia_score_all~idh_codel_subtype)
plot_violin(tumor_classes,"hypoxia_score_all", "score", stat.test, stat.test.grade, stat.test.subtype, "hypoxia_score", label="p.adj", hide.ns=F)

# beta values
# calculate median of variable probes
tumor_classes$CA9_median_var_probes = apply(tumor_classes[c("CA9_beta_cg20610181", "CA9_beta_cg09566069", "CA9_beta_cg06908460", "CA9_beta_cg13849253")], 1, median)
tumor_classes$VEGFA_median_var_probes = apply(tumor_classes[c("VEGFA_beta_cg25343661", "VEGFA_beta_cg12279019")], 1, median)
tumor_classes$ADM_median_var_probes = apply(tumor_classes[c("ADM_beta_cg08259810", "ADM_beta_cg05149386", "ADM_beta_cg19315081", "ADM_beta_cg20741987",
                                                            "ADM_beta_cg03548673", "ADM_beta_cg02843237", "ADM_beta_cg10044466", "ADM_beta_cg06875754",
                                                            "ADM_beta_cg05727225")], 1, median)
tumor_classes$PDK1_median_var_probes = apply(tumor_classes[c("PDK1_beta_cg17679246", "PDK1_beta_cg21044834")], 1, median)

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(CA9_median_var_probes ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(CA9_median_var_probes ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(CA9_median_var_probes~idh_codel_subtype)
plot_violin(tumor_classes,"CA9_median_var_probes", "methylation", stat.test, stat.test.grade, stat.test.subtype, "CA9_beta_filtered_median", label="p.adj", hide.ns=F)

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(ADM_median_var_probes ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(ADM_median_var_probes ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(ADM_median_var_probes~idh_codel_subtype)
plot_violin(tumor_classes,"ADM_median_var_probes", "methylation", stat.test, stat.test.grade, stat.test.subtype, "ADM_beta_filtered_median", label="p.adj", hide.ns=F)

stat.test.grade = tumor_classes[tumor_classes$idh_codel_subtype %in% c("IDHmut-codel", "IDHmut-non-codel"),] %>% group_by(idh_codel_subtype) %>% wilcox_test(CA9_beta_cg20610181 ~ grade)
stat.test.subtype = tumor_classes %>% group_by(grade) %>% wilcox_test(CA9_beta_cg20610181 ~ idh_codel_subtype)
stat.test = tumor_classes %>% wilcox_test(CA9_beta_cg20610181~idh_codel_subtype)
plot_violin(tumor_classes,"CA9_beta_cg20610181", "methylation", stat.test, stat.test.grade, stat.test.subtype, "CA9_beta_cg20610181", label="p.adj", hide.ns=F)

################################
# gene set filtering
# example for hypoxia gene set, other processed similarly

hypoxia = c("CA9", "VEGFA", "ADM", "PDK1", "ENO1", "HK2", "SLC2A1", "AK3", "PFKFB3", "CCNG2")

hypoxia = expr_log2[unlist(hypoxia), tumor_classes[tumor_classes$idh_codel_subtype!="IDHmut-codel", "sample"]]


get_correlating_genes = function(correlations, cutoff) {
  n_over_cutoff = rowSums(correlations>cutoff)
  if (var(n_over_cutoff) == 0 ){
    return(correlations)
  } else {
    min_n = which(n_over_cutoff==min(n_over_cutoff))
    if (length(min_n)==1) {
      min_n_rowsum = names(min_n)
    } else if (length(min_n)>1) {
      min_n_rowsum = names(which.min(rowSums(correlations[min_n, ])))
    }
    
    gfd = correlations[!rownames(correlations)==min_n_rowsum, !colnames(correlations)==min_n_rowsum]
    return(get_correlating_genes(gfd, cutoff))
  }
}

plot_correlation_heatmap = function(df, title="", method="pearson", usecollabel=F, print_cor=F, printgenes=T, stop_in_return=F) {
  
  df <- df[apply(df, 1, sd) > 0, ]
  Correlations = cor(t(as.matrix(df)), method=method)
  
  # set all values below 0.4 to 0
  Correlations[Correlations<=0.4] = 0
  
  colrows = rep("black", nrow(Correlations))
  colcols = colrows
  if (print_cor){
    cor_text = round(Correlations, 2)
    
  } else {
    cor_text=matrix("", ncol=ncol(Correlations), nrow=nrow(Correlations))
  }
  
  if (printgenes) {
    genenames = rownames(Correlations)
  } else {
    genenames = ""
  }

  hm = heatmap.2(Correlations, scale="none", trace="none", density.info="none",
                 cellnote=cor_text, notecol="black", notecex = 0.4, # 0.7
                 colRow = colrows, colCol = colcols,
                 breaks=seq(-1,1, 0.02),
                 margins=c(5,5),
                 labRow=genenames, labCol=genenames,
                 cexCol=0.7, cexRow=0.7,
                 col=bluered(100), key.par=list(cex=0.5), keysize = 0.75,
                 main=paste0(title, " - ", method)) #, breaks=seq(-1, 1, length.out=101))
  
  if (stop_in_return) {
    return(list(Correlations, hm))
  }
  
  hr_all = hclust(as.dist(1-Correlations))
  
  if (usecollabel) {
    usecollabel=colnames(df)
    usecollabel = substr(usecollabel, 1, 16)
  }
  
  hm2 = heatmap.2(df, density.info = "none", trace = "none", scale = "row",
                  Rowv=as.dendrogram(hr_all),
                  #Rowv=F,
                  Colv=T, 
                  colRow = colrows,
                  labCol = usecollabel, cexCol = 0.5,
                  cexRow=0.7,
                  margins=c(5, 5),
                  col=bluered(100), key.par=list(cex=0.5), keysize = 0.75)#,
                  #main=title)
  
  return(hm)
}

plot_correlation_heatmap(hypoxia, "TCGA IDHmut astro + GBM", print_cor=T, method="spearman")
plot_correlation_heatmap(hypoxia, "TCGA IDHmut astro + GBM", print_cor=T, method="pearson")



get_cluster_corr_genes = function(geneset, method1, method2) {
  l1 = plot_correlation_heatmap(geneset, "TCGA IDHmut astro + GBM", print_cor=T, method=method1, stop_in_return=T)
  l2 = plot_correlation_heatmap(geneset, "TCGA IDHmut astro + GBM", print_cor=T, method=method2, stop_in_return=T)

  cor1 = l1[[1]]
  hm1 = l1[[2]]
  cor2 = l2[[1]]
  hm2 = l2[[2]]

  hm1mat = hm1$carpet
  hm2mat = hm2$carpet

  cluster_threshold = round(0.05*nrow(hm1mat))
  if(cluster_threshold>20) {
    cluster_threshold=20
  } else if (cluster_threshold<4) {
    cluster_threshold=4
  }

  selected_genes = list()

  for(i in 1:nrow(hm1mat)) {
    diagonalgene = rownames(hm1mat)[i]
    # find cluster around
    correlated_ix = c(i)
  
    forward = i+1
    backward = i-1
    while(forward <= nrow(hm1mat) && hm1mat[i, forward] > 0.4) {
      correlated_ix = c(correlated_ix, forward)
      forward = forward+1
    }
    while(backward >= 1 && hm1mat[i, backward] > 0.4) {
      correlated_ix = c(correlated_ix, backward)
      backward = backward-1
    }
  
    correlated_mat = hm1mat[correlated_ix, correlated_ix]
    
  
    if(all(correlated_mat > 0.4) && length(correlated_ix) >= cluster_threshold) {
      print(i)
      print(c("seedcluster size: ", length(correlated_ix))) # IF this seed cluster is bigger than X, plot heatmaps separately!
      correlated_genes = rownames(correlated_mat)
      # take genes that also correlate with other correlation method
      corr2 = hm2mat[correlated_genes,correlated_genes]
      corr2genes =  rownames(corr2)[apply(corr2, 1, function(row) all(row > 0.4))]
    
      if (length(corr2genes) >= cluster_threshold) {
        print(i)
        # take all other genes that correlate with all of these > 0.4
        c = hm1mat[!(rownames(hm1mat) %in% corr2genes), corr2genes]
        other_genes = rownames(c)[apply(c, 1, function(row) all(row > 0.4))]
      
        # take genes that correlate with other correlation method
        c2 = hm2mat[other_genes, corr2genes]
        if (length(other_genes)==1) { # is a vector
          other_genes = if(all(c2>0.4)) other_genes else NULL
        } else if (length(other_genes)>0) { # is a matrix
          other_genes = rownames(c2)[apply(c2, 1, function(row) all(row > 0.4))]
          print(c("other genes: ", length(other_genes)))
        } else {
          other_genes = NULL
        }
        # add if there are any genes
        if(length(other_genes)>0) {
          selected_genes[[diagonalgene]] = c(corr2genes, other_genes)
        } else {
          selected_genes[[diagonalgene]] = corr2genes
        } 
      }
    
    }
  
  }
  
  return(selected_genes)
}


sel_genes_s_hypoxia = get_cluster_corr_genes(hypoxia, method1="spearman", method2="pearson")
sel_genes_p_hypoxia = get_cluster_corr_genes(hypoxia, method1="pearson", method2="spearman")
sel_genes_hypoxia = intersect(unlist(sel_genes_s_hypoxia), unlist(sel_genes_p_hypoxia))

hm_s = plot_correlation_heatmap(hypoxia[sel_genes_hypoxia,], "TCGA IDHmut astro + GBM", print_cor=T, method="spearman")
hm_p = plot_correlation_heatmap(hypoxia[sel_genes_hypoxia,], "TCGA IDHmut astro + GBM", print_cor=T, method="pearson")



#####################################
# Capper data 

load("TCGA/hypoxia/Capper_original_cohort_filtered_for_Ca9.RData")
load("TCGA/hypoxia/Capper_original_cohort_sample_anno.RData")

# e.g. for CA9 probes
df = anno[, c("title","geo_accession","characteristics_ch1", "methylation class:ch1")]
df$cg13849253 = NA
df$cg06908460 = NA
df$cg09566069 = NA
df$cg20610181 = NA

for (i in 1:nrow(df)) {
  df[i, "cg13849253"] = Capper_original_filtered[grep(df$geo_accession[i], rownames(Capper_original_filtered)), "cg13849253"]
  df[i, "cg06908460"] = Capper_original_filtered[grep(df$geo_accession[i], rownames(Capper_original_filtered)), "cg06908460"]
  df[i, "cg09566069"] = Capper_original_filtered[grep(df$geo_accession[i], rownames(Capper_original_filtered)), "cg09566069"]
  df[i, "cg20610181"] = Capper_original_filtered[grep(df$geo_accession[i], rownames(Capper_original_filtered)), "cg20610181"]
}


df$median_variable_probes = apply(df[,5:8], 1, median)
df$mean_variable_probes = apply(df[,5:8], 1, mean)
df$`methylation class:ch1`

classes = c("O IDH","A IDH","A IDH, HG","CONTR, ADENOPIT","CONTR, CEBM",
            "CONTR, HEMI","CONTR, HYPTHAL","CONTR, PINEAL","CONTR, PONS",
            "CONTR, WM","GBM, MES","GBM, RTK I","GBM, RTK II")

df = df[df$`methylation class:ch1` %in% classes, ]

df$`methylation class:ch1` = factor(df$`methylation class:ch1`, levels=classes) 


p = ggplot(df, aes(y=median_variable_probes, x=`methylation class:ch1`))+
  geom_violin()+
  geom_point(position=position_jitter(seed=1,width=0.15), alpha=0.4, size=1)+
  tehem_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# statistical tests
df$group = "control"
df[df$'methylation class:ch1' %in% c("A IDH", "A IDH, HG", "O IDH"), "group"] = "IDHmut"
df[df$'methylation class:ch1' %in% c("GBM, MES", "GBM, RTK I", "GBM, RTK II"), "group"] = "GBM"

stat.test = df %>% wilcox_test(median_variable_probes~group)
