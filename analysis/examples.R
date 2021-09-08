library(cola)
library(GetoptLong)

################### Golub_leukemia ########################

library(golubEsets)
data(Golub_Merge)
m = exprs(Golub_Merge)
colnames(m) = paste0("sample_", colnames(m))
anno = pData(Golub_Merge)

m[m <= 1] = NA
m = log10(m)

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
rn = rownames(m)
m = normalize.quantiles(m)
colnames(m) = cn
rownames(m) = rn

rh = hierarchical_partition(m, cores = 4, 
	anno = anno[, c("ALL.AML"), drop = FALSE],
    anno_col = c("ALL" = "red", "AML" = "blue"))
saveRDS(rh, file = "Golub_leukemia_cola_rh.rds")

#################### Ritz_ALL #########################

library(ALL)
data(ALL)

m = exprs(ALL)
anno = pData(ALL)
anno = anno[, c("sex", "age", "BT")]

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
rn = rownames(m)
m = normalize.quantiles(m)
colnames(m) = cn
rownames(m) = rn

rh = hierarchical_partition(m, cores = 4, anno = anno)
saveRDS(rh, file = "Ritz_ALL_cola_rh.rds")

##################### TCGA_GBM_microarray #######################

library(RColorBrewer)

m = read.table("https://jokergoo.github.io/cola_examples/TCGA_GBM/unifiedScaled.txt", header = TRUE, row.names = 1, check.names = FALSE)
m = as.matrix(m)

subtype = read.table("https://jokergoo.github.io/cola_examples/TCGA_GBM/TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, 
    check.names = FALSE, stringsAsFactors = FALSE)
subtype = structure(unlist(subtype[1, -(1:2)]), names = colnames(subtype)[-(1:2)])
subtype_col = structure(seq_len(4), names = unique(subtype))

m = m[, names(subtype)]
m = adjust_matrix(m)


library(preprocessCore)
cn = colnames(m)
rn = rownames(m)
m = normalize.quantiles(m)
colnames(m) = cn
rownames(m) = rn

rh = hierarchical_partition(m, cores = 4, anno = subtype, anno_col = subtype_col)
saveRDS(rh, file = "TCGA_GBM_microarray_cola_rh.rds")


###################### HSMM_single_cell ######################
library(RColorBrewer)
library(HSMMSingleCell)
data(HSMM_expr_matrix)
data(HSMM_sample_sheet)

anno = HSMM_sample_sheet[, c("Hours", "Media", "State")]
anno_col = list(
    Hours = structure(brewer.pal(9, "Blues")[c(2, 4, 6, 8)], names = c("0", "24", "48", "72")),
    Media = c("GM" = "orange", "DM" = "purple"),
    State = c("1" = "red", "2" = "blue", "3" = "green"))

m = adjust_matrix(log10(HSMM_expr_matrix + 1))
gt = readRDS(url("https://jokergoo.github.io/cola_examples/HSMM_single_cell/gene_type_gencode_v17.rds"))
m = m[gt[rownames(m)] == "protein_coding", , drop = FALSE]

rh = hierarchical_partition(m, cores = 4,
	anno = anno,
    anno_col = anno_col
)
saveRDS(rh, file = "HSMM_single_cell_cola_rh.rds")


####################### PBMC #######################

library(Seurat)
pbmc.data = Read10X(data.dir = "filtered_gene_bc_matrices")
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc = NormalizeData(pbmc)

pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

Seurat_class = as.character(Idents(pbmc))

mat = as.matrix(GetAssayData(pbmc))
mat = adjust_matrix(mat)

rh = hierarchical_partition(mat, subset = 500, cores = 4, 
	anno = data.frame(Seurat_class = Seurat_class))
saveRDS(rh, file = "PBMC_cola_rh.rds")


####################################################
library(scRNAseq)

process_counts = function(data, column = NULL) {
    mat = assays(data)$counts
    mat = as.matrix(mat)
    s = colSums(mat)
    fa = s/mean(s)
    for(i in 1:ncol(mat)) mat[, i]/fa[i]
    mat = adjust_matrix(log2(mat + 1))

    anno = NULL
    if(!is.null(column)) {
        anno = colData(data)
        anno = as.data.frame(anno)
        anno = anno[, column, drop = FALSE]
    }

    list(mat = mat, anno = anno)
}

######## BaronPancreas_mouse ############
data = BaronPancreasData('mouse')
lt = process_counts(data)
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "BaronPancreas_mouse_cola_rh.rds")

############ BuettnerESC #############
data = BuettnerESCData()
lt = process_counts(data)
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "BuettnerESC_cola_rh.rds")

############ DarmanisBrain ###############
data = DarmanisBrainData()
lt = process_counts(data, c("age", "cell.type"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "DarmanisBrain_cola_rh.rds")

############ GrunHSC ###########
data = GrunHSCData()
lt = process_counts(data, c("sample"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "GrunHSC_cola_rh.rds")


############ GrunPancreas #################
data = GrunPancreasData()
lt = process_counts(data, c("donor"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "GrunPancreas_cola_rh.rds")


############# KolodziejczykESC #################
data = KolodziejczykESCData()
lt = process_counts(data)
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "KolodziejczykESC_cola_rh.rds")

############ LaMannoBrain_human_es #################
data = LaMannoBrainData('human-es')
lt = process_counts(data, c("Cell_type", "Timepoint"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "LaMannoBrain_human_es_cola_rh.rds")

########### LaMannoBrain_human_embryo ##################
data = LaMannoBrainData('human-embryo')
lt = process_counts(data, c("Cell_type", "Timepoint"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "LaMannoBrain_human_embryo_cola_rh.rds")

########## LaMannoBrain_human_ips #################
data = LaMannoBrainData('human-ips')
lt = process_counts(data, c("Cell_type", "Timepoint"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "LaMannoBrain_human_ips_cola_rh.rds")

######### LaMannoBrain_mouse_adult ###############
data = LaMannoBrainData('mouse-adult')
lt = process_counts(data, c("Cell_type"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "LaMannoBrain_mouse_adult_cola_rh.rds")


######## LaMannoBrain_mouse_embryo #############
data = LaMannoBrainData('mouse-embryo')
lt = process_counts(data, c("Cell_type", "Timepoint"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "LaMannoBrain_mouse_embryo_cola_rh.rds")

########### LawlorPancreas ################
data = LawlorPancreasData()
lt = process_counts(data, c("cell.type", "disease"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "LawlorPancreas_cola_rh.rds")

########### LengESC ##############
data = LengESCData()

mat = assays(data)$normcounts
mat = log2(mat + 1)
mat = adjust_matrix(mat)

anno = colData(data)
anno = as.data.frame(anno)

rh = hierarchical_partition(mat, subset = 500, cores = 4, anno = anno)
saveRDS(rh, file = "LengESC_cola_rh.rds")

########### ReprocessedTh2 ###############
data = ReprocessedTh2Data()
mat = assays(data)$rsem_tpm
mat = log2(mat + 1)
mat = adjust_matrix(mat)

anno = colData(data)[, c("NREADS")]
anno = as.data.frame(anno)

rh = hierarchical_partition(mat, subset = 500, cores = 4, anno = anno)
saveRDS(rh, file = "ReprocessedTh2_cola_rh.rds")


############ MessmerESC ################
data = MessmerESCData()
lt = process_counts(data, c("phenotype"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "MessmerESC_cola_rh.rds")


############ MuraroPancreas #################
data = MuraroPancreasData()
lt = process_counts(data, c("label"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "MuraroPancreas_cola_rh.rds")

########### NestorowaHSC ##############
data = NestorowaHSCData()
lt = process_counts(data)
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4)
saveRDS(rh, file = "NestorowaHSC_cola_rh.rds")

############ PollenGlia ################
data = PollenGliaData()
lt = process_counts(data, c("Age", "Inferred.Cell.Type"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "PollenGlia_cola_rh.rds")

########### RichardTCell ###################
data = RichardTCellData()
lt = process_counts(data, c("age"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "RichardTCell_cola_rh.rds")

############ RomanovBrain #################
data = RomanovBrainData()
lt = process_counts(data, "level1.class")
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "RomanovBrain_cola_rh.rds")

############# SegerstolpePancreas ################
data = SegerstolpePancreasData()
lt = process_counts(data, c("cell.type", "disease", "age"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "SegerstolpePancreas_cola_rh.rds")

############## TasicBrain #################
data = TasicBrainData()
lt = process_counts(data, c("broad_type"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "TasicBrain_cola_rh.rds")

######## ReprocessedAllen #################
data = ReprocessedAllenData()
mat = assays(data)$rsem_tpm
mat = log2(mat + 1)
mat = adjust_matrix(mat)

anno = colData(data)[, c("driver_1_s", "dissection_s", "Core.Type", "Primary.Type", "Secondary.Type")]
anno = as.data.frame(anno)

rh = hierarchical_partition(mat, subset = 500, cores = 4, anno = anno)
saveRDS(rh, file = "ReprocessedAllen_cola_rh.rds")

########### XinPancreas ################
data = XinPancreasData()
mat = as.matrix(assays(data)$rpkm)
mat = log2(mat + 1)
mat = adjust_matrix(mat)

anno = colData(data)[, c("age", "cell.type")]
anno = as.data.frame(anno)

rh = hierarchical_partition(mat, subset = 500, cores = 4, anno = anno)
saveRDS(rh, file = "XinPancreas_cola_rh.rds")

############ ZeiselBrain #################
data = ZeiselBrainData()
lt = process_counts(data, c("level1class", "level2class"))
rh = hierarchical_partition(lt$mat, subset = 500, cores = 4, anno = lt$anno)
saveRDS(rh, file = "ZeiselBrain_cola_rh.rds")

######## Fluidigm ################
data = ReprocessedFluidigmData()

anno = colData(data)[, c("Biological_Condition", "Coverage_Type", "Cluster1", "Cluster2")]
anno = as.data.frame(anno)

mat = assays(data)$rsem_tpm[, anno$Coverage_Type == "High"]
count = assays(data)$rsem_counts[, anno$Coverage_Type == "High"]
l = apply(count, 1, function(x) sum(x > 0)/length(x) > 0.1)
mat = log2(mat[l, ] + 1)
mat = adjust_matrix(mat)

anno = anno[anno$Coverage_Type == "High", ]

rh = hierarchical_partition(mat, cores = 4, anno = anno)
saveRDS(rh, file = "Fluidigm_cola_rh.rds")

######### GSE90496 ##############

tumor_type = read.table("GSE90496_tumor_type.csv", sep = ";", header = TRUE, comment.char = "")
tumor_col = structure(unique(tumor_type$color), names = unique(tumor_type$tumor_type))
tumor_type = structure(tumor_type$tumor_type, names = tumor_type$subtype)

lt = readRDS("GSE90496_islands_processed.rds")
mat = as.matrix(lt$mat)
anno = lt$anno
anno$subtype2 = gsub("\\s+", " ", anno$subtype2)
anno$subtype2 = gsub("\\s+$", "", anno$subtype2)

anno = data.frame(meth_class = anno$subtype2, tumor_type = tumor_type[anno$subtype2])

anno_col = list(
    tumor_type = tumor_col
)

mat = adjust_matrix(mat)

library(matrixStats)
rh = hierarchical_partition(mat, cores = 4, 
    top_value_method = c("SD", "ATC"), max_k = 8, 
    partition_method = c("kmeans", "skmeans"),
    scale_rows = FALSE, anno = anno, top_n = 1000, 
    subset = 500, group_diff = 0.25, min_n_signatures = 1000,
    filter_fun = function(mat) {
        s = rowSds(mat)
        order(-s)[1:30000]
    })

saveRDS(rh, file = "GSE90496_cola_rh.rds")


# # generate code for TCGA 450K array datasets ###########
# files = list.files("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data", full.names = TRUE)
# cohort = gsub("\\.sampleMap.*$", "", basename(files))
# cohort = gsub("\\.", "_", cohort)

# for(i in seq_along(files)) {
#     mat = read.table(files[i], header = TRUE, row.names = 1, nrows = 10)
#     if(ncol(mat) > 50) {
#         qqcat('
# ######### @{cohort[i]}_methylation ############
# library(cola)
# library(matrixStats)

# mat = read.table("@{files[i]}", header = TRUE, row.names = 1)
# mat = adjust_matrix(as.matrix(mat))
# rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"), 
#                             max_k = 8, partition_method = c("kmeans", "skmeans"),
#                             scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
#                             filter_fun = function(mat) {
#                                 s = rowSds(mat)
#                                 order(-s)[1:30000]
#                             })
# saveRDS(rh, file = "@{cohort[i]}_methylation_cola_rh.rds")
# ')}
# }


######### TCGA_ACC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.ACC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_ACC_methylation_cola_rh.rds")

######### TCGA_BLCA_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.BLCA.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_BLCA_methylation_cola_rh.rds")

######### TCGA_BRCA_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.BRCA.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_BRCA_methylation_cola_rh.rds")

######### TCGA_CESC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.CESC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_CESC_methylation_cola_rh.rds")

######### TCGA_COAD_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.COAD.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_COAD_methylation_cola_rh.rds")

######### TCGA_COADREAD_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.COADREAD.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_COADREAD_methylation_cola_rh.rds")

######### TCGA_ESCA_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.ESCA.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_ESCA_methylation_cola_rh.rds")

######### TCGA_GBM_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.GBM.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_GBM_methylation_cola_rh.rds")

######### TCGA_GBMLGG_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.GBMLGG.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_GBMLGG_methylation_cola_rh.rds")

######### TCGA_HNSC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.HNSC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_HNSC_methylation_cola_rh.rds")

######### TCGA_KICH_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.KICH.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_KICH_methylation_cola_rh.rds")

######### TCGA_KIRC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.KIRC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_KIRC_methylation_cola_rh.rds")

######### TCGA_KIRP_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.KIRP.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_KIRP_methylation_cola_rh.rds")

######### TCGA_LAML_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.LAML.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_LAML_methylation_cola_rh.rds")

######### TCGA_LGG_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.LGG.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_LGG_methylation_cola_rh.rds")

######### TCGA_LIHC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.LIHC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_LIHC_methylation_cola_rh.rds")

######### TCGA_LUAD_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.LUAD.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_LUAD_methylation_cola_rh.rds")

######### TCGA_LUNG_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.LUNG.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_LUNG_methylation_cola_rh.rds")

######### TCGA_LUSC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.LUSC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_LUSC_methylation_cola_rh.rds")

######### TCGA_MESO_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.MESO.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_MESO_methylation_cola_rh.rds")

######### TCGA_PAAD_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.PAAD.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_PAAD_methylation_cola_rh.rds")

######### TCGA_PCPG_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.PCPG.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_PCPG_methylation_cola_rh.rds")

######### TCGA_PRAD_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.PRAD.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_PRAD_methylation_cola_rh.rds")

######### TCGA_READ_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.READ.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_READ_methylation_cola_rh.rds")

######### TCGA_SARC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.SARC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_SARC_methylation_cola_rh.rds")

######### TCGA_SKCM_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.SKCM.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_SKCM_methylation_cola_rh.rds")

######### TCGA_STAD_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.STAD.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_STAD_methylation_cola_rh.rds")

######### TCGA_TGCT_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.TGCT.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_TGCT_methylation_cola_rh.rds")

######### TCGA_THCA_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.THCA.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_THCA_methylation_cola_rh.rds")

######### TCGA_THYM_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.THYM.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_THYM_methylation_cola_rh.rds")

######### TCGA_UCEC_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.UCEC.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_UCEC_methylation_cola_rh.rds")

######### TCGA_UCS_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.UCS.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_UCS_methylation_cola_rh.rds")

######### TCGA_UVM_methylation ############
library(cola)
library(matrixStats)

mat = read.table("/omics/groups/OE0246/internal/guz/cola_hc/data/TCGA_methylation/data/TCGA.UVM.sampleMap__HumanMethylation450.gz", header = TRUE, row.names = 1)
mat = adjust_matrix(as.matrix(mat))
rh = hierarchical_partition(mat, cores = 4, top_value_method = c("SD", "ATC"),
                            max_k = 8, partition_method = c("kmeans", "skmeans"),
                            scale_rows = FALSE, top_n = 1000, subset = 500, group_diff = 0.25, min_n_signatures = 1000,
                            filter_fun = function(mat) {
                                s = rowSds(mat)
                                order(-s)[1:30000]
                            })
saveRDS(rh, file = "TCGA_UVM_methylation_cola_rh.rds")
