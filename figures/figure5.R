

library(cola)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(Seurat)
pbmc.data = Read10X(data.dir = "~/workspace/cola_hc/filtered_gene_bc_matrices")
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
nc = length(unique(Seurat_class))
Seurat_col = structure(Polychrome::alphabet.colors(nc), names = unique(Seurat_class))
Seurat_col = Seurat_col[sort(names(Seurat_col))]
mat = as.matrix(GetAssayData(pbmc))
mat = adjust_matrix(mat)


### perform hierarchical partitioning
set.seed(123)
if(file.exists("~/workspace/cola_hc/PBMC_cola_hierarchical_partition.rds")) {
	rh = readRDS("~/workspace/cola_hc/PBMC_cola_hierarchical_partition.rds")
	res = readRDS("~/workspace/cola_hc/PBMC_cola_consensus_partition.rds")
} else {
	rh = hierarchical_partition(mat, subset = 500, cores = 4,
		anno = data.frame(Seurat_class = Seurat_class), anno_col = Seurat_col,
		min_n_signatures = 200)
	saveRDS(rh, file = "~/workspace/cola_hc/PBMC_cola_hierarchical_partition.rds")

	res = consensus_partition_by_down_sampling(mat, subset = 500, cores = 4, max_k = 10,
		anno = data.frame(Seurat_class = Seurat_class), anno_col = Seurat_col)
	saveRDS(res, file = "~/workspace/cola_hc/PBMC_cola_consensus_partition.rds")
}

cola_cl = get_classes(res, k = 4)[, 1]
rh@list[[1]]@full_anno$cola_4_groups = cola_cl
rh@list[[1]]@anno$cola_4_groups = rep(NA, 500)
rh@list[[1]]@anno_col$cola_4_groups = cola_opt$color_set_2[as.character(1:4)]

### t-SNE plot
p1 = ~{set.seed(123); dimension_reduction(rh, method = "t-SNE", cex = 0.5, top_n = 1000)}
p2 = ~{set.seed(123); 
	# just to change the legend title
	rh2 = rh
	rh2@list[[1]]@full_anno$Class = rh2@list[[1]]@full_anno$Seurat_class
	rh2@list[[1]]@anno$Class = rep(NA, 500)
	rh2@list[[1]]@anno_col$Class = {x <- rh@list[[1]]@anno_col$Seurat_class; x[sort(names(x))]}
	dimension_reduction(rh2, method = "t-SNE", cex = 0.5, top_n = 1000, color_by = "Class")
}
p3 = ~{
	rh2 = rh
	rh2@list[[1]]@full_anno$Class = rh2@list[[1]]@full_anno$cola_4_groups
	rh2@list[[1]]@anno$Class = rep(NA, 500)
	rh2@list[[1]]@anno_col$Class = {x <- rh@list[[1]]@anno_col$cola_4_groups; x[sort(names(x))]}
	set.seed(123); dimension_reduction(rh2, method = "t-SNE", cex = 0.5, top_n = 1000, color_by = "Class")
}


### signature heatmap
sig_tb = get_signatures(rh, plot = FALSE)

cl = get_classes(rh)
ht = Heatmap(t(scale(t(mat[sig_tb$which_row, ]))), name = "Z-score", column_title = qq("@{length(unique(cl))} groups, @{nrow(sig_tb)} signatures under FDR < 0.05"),
	col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
	top_annotation = HeatmapAnnotation(Class = cl, col = list(Class = rh@subgroup_col),
		annotation_legend_param = list(Class = list(ncol = 2))),
	bottom_annotation = HeatmapAnnotation(Seurat_class = Seurat_class, cola_4_groups = cola_cl,
		col = list(Seurat_class = Seurat_col, cola_4_groups = cola_opt$color_set_2),
		annotation_legend_param = list(Seurat_class = list(ncol = 2), cola_4_groups = list(ncol = 2))),
	show_row_names = FALSE, show_column_names = FALSE, show_row_dend = FALSE,
	cluster_columns = cola:::calc_dend(rh, mat = mat), 
	column_split = length(unique(cl)), row_km = 7, row_title = paste0("km", 1:7), row_title_rot = 0
)
p4 = grid.grabExpr(ht <- draw(ht, merge_legend = TRUE))

### stability

library(clue)
rh_list = readRDS("PBMC_cola_hierarchical_partition_rh_list.rds")
cl = cl_ensemble(list = lapply(rh_list, as.cl_hard_partition))
clr = cl_class_ids(cl_consensus(cl))
clt = lapply(rh_list, function(x) relabel_class(x, clr, return_map = FALSE))
clm = do.call(cbind, clt)

map = relabel_class(get_classes(rh), clr)
map = map[grepl("^0", names(map))]

nl = unique(as.vector(clm))
col = structure(Polychrome::alphabet.colors(length(nl)), names = nl)
col[map] = rh@subgroup_col[names(map)]
p_stability = grid.grabExpr(
	draw(Heatmap(clm, row_order = do.call(order, clt), col = col,
		cluster_columns = hclust(cl_dissimilarity(cl_ensemble(list = lapply(clt, as.cl_hard_partition)))),
		show_column_dend = FALSE, show_heatmap_legend = FALSE,
		column_title = "Stability of classification\nfrom 20 runs")
	)
)

### correspondance

tb = data.frame(cola_class = as.character(get_classes(rh)), Seurat_class = as.character(Seurat_class))

foo = table(tb$Seurat_class, tb$cola_class)
foo = as.matrix(foo)

foo = foo[rev(c("2", "0", "4", "5", "7", "1", "3", "6")),
          rev(c("01", "041", "0211", "0212", "022", "03", "042"))]

w1 = rowSums(foo)
w2 = colSums(foo)

cum1 = cumsum(w1)
cum2 = cumsum(w2)

p_agreement = grid.grabExpr({
	grid.newpage()
	pushViewport(viewport(yscale = c(0, sum(w1)), height = unit(1, "npc") - unit(4, "mm") - unit(1, "cm"), y = unit(2, "mm"), just = "bottom"))
	for(i in seq_len(nrow(foo))) {
	    for(j in seq_len(ncol(foo))) {
	        if(foo[i, j] == 0) next
	        x11 = cum1[i] - w1[i] + sum(foo[i, 1:j]) - foo[i, j]
	        x12 = x11 + foo[i, j]

	        x21 = cum2[j] - w2[j] + sum(foo[1:i, j]) - foo[i, j]
	        x22 = x21 + foo[i, j]

	        pt1 = circlize:::get_bezier_points(x11, 0.2, x21, 0.6, xlim = c(0, sum(w1)), ylim = c(0, 1))
	        pt2 = circlize:::get_bezier_points(x12, 0.2, x22, 0.6, xlim = c(0, sum(w1)), ylim = c(0, 1))

	        grid.polygon(
	            c(0.2, 0.2, pt2[, 2], rev(pt1[, 2]), 0.6), default.unit = "native",
	            c(x11, x12, pt2[, 1], rev(pt1[, 1]), x11),
	            gp = gpar(fill = add_transparency(Seurat_col[rownames(foo)][i], 0.25), col = NA)
	        )

	    }
	}

	grid.rect(y = cum1, x = unit(0.2, "native"), height = w1, width = unit(5, "mm"), just = c("right", "top"), default.units = "native",
	    gp = gpar(fill = Seurat_col[rownames(foo)], col = NA))
	grid.text("Seurat classification", x = unit(0.2, "native") - unit(7, "mm"), y = 0.5, rot = 90, just = "bottom")

	grid.rect(y = cum2, x = 0.6, height = w2, width = unit(5, "mm"), just = c("left", "top"), default.units = "native",
	    gp = gpar(fill = rh@subgroup_col[colnames(foo)], col = NA))
	grid.text("cola HCP classification", x = unit(0.6, "native") + unit(7, "mm"), y = 0.5, rot = 90, just = "top")

	rh_col = rh@subgroup_col[colnames(foo)]
	lgd = packLegend(
		Legend(title = "Seurat", at = names(Seurat_col), legend_gp = gpar(fill = Seurat_col)),
		Legend(title = "cola HCP", at = names(rh_col), legend_gp = gpar(fill = rh_col))
	)

	draw(lgd, x = unit(0.6, "npc") + unit(15, "mm"), just = "left")

})

### puting all the figure together

library(cowplot)

pdf("figure5.pdf", width = 14, height = 10)
plot_grid( plot_grid(p1, p2, p3, nrow = 1), 
	plot_grid(p_agreement, p4, p_stability, rectGrob(gp = gpar(fill = NA, col = NA)), 
		nrow = 1, rel_widths = c(1.5, 2, 0.6, 1)), rel_heights = c(1, 1.2), nrow = 2)
dev.off()



concor = function(x, y) {
	sum(x == y)/length(x)
}

s = NULL
for(i in 1:19) {
	for(j in (i+1):20) {
		s = c(s, concor(clt[[i]], clt[[j]]))
	}
}


# function `overall_classification_agreement()` is from figure6.R
overall_classification_agreement(tb$cola_class, tb$Seurat_class)
