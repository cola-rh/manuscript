
fig_width = 20
fig_height = 12

library(cola)
library(ComplexHeatmap)
library(GetoptLong)
library(circlize)
library(HSMMSingleCell)
data(HSMM_expr_matrix)
data(HSMM_sample_sheet)


# `HSMM_expr_matrix` is a FPKM matrix
m = adjust_matrix(log10(HSMM_expr_matrix + 1), verbose = FALSE)
gt = readRDS("~/workspace/cola_hc/gene_type_gencode_v17.rds")
m = m[gt[rownames(m)] == "protein_coding", , drop = FALSE]


s = ATC(m, cores = 4)

m2 = m[order(-s)[1:1000], ]
m2 = t(scale(t(m2)))

ht = Heatmap(m2, col = colorRamp2(c(-1.5, 0, 1.5), c("green", "white", "red")), name = "z-score",
	show_row_dend = FALSE, show_column_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE,
	column_title = "Top 1000 genes with the highest ATC scores", use_raster = TRUE,
	heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"))
p1 = grid.grabExpr(draw(ht, heatmap_legend_side = "bottom", padding = unit(c(4, 4, 4, 4), "mm")))


set.seed(1234)
if(file.exists("HSMM_cola_consensus_partition.rds")) {
	res = readRDS("HSMM_cola_consensus_partition.rds")
} else {
	res = consensus_partition(m, top_n = 1000, cores = 4, max_k = 8)
	saveRDS(res, file = "HSMM_cola_consensus_partition.rds")
}

plot_ecdf = function(object) {
	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "consensus value (x)", ylab = "P(X <= x)", main = "eCDF curve")
	for(i in seq_along(object@k)) {
		consensus_mat = get_consensus(object, k = object@k[i])
		f = ecdf(consensus_mat[lower.tri(consensus_mat)])
		x = seq(0, 1, length = 100)
		y = f(x)
		x = c(0, x)
		y = c(0, y)
		lines(x, y, col = i)
	}
	legend("bottomright", pch = 15, legend = paste0("k = ", object@k), col = seq_along(object@k), ncol = 2)
}
p2 = ~plot_ecdf(res)

stats = get_stats(res)

p3 = ~plot(stats[, "k"], stats[, "1-PAC"], type = "b", xlab = "Number of subgroups (k)", ylab = "1-PAC")
p4 = ~plot(stats[, "k"], stats[, "mean_silhouette"], type = "b", xlab = "Number of subgroups (k)", ylab = "Mean Silhouette")
p5 = ~plot(stats[, "k"], stats[, "concordance"], type = "b", xlab = "Number of subgroups (k)", ylab = "Concordance")


make_single_plot = function(object, k) {	
	
	class_df = get_classes(object, k)
	class_ids = class_df$class

	membership_mat = get_membership(object, k)
	col_fun = colorRamp2(c(0, 1), c("white", "red"))

	membership_each = get_membership(object, k, each = TRUE)
	membership_each = t(membership_each)
	mat_col_od = cola:::column_order_by_group(factor(class_ids, levels = sort(unique(class_ids))), membership_each)


	ht = Heatmap(membership_each, name = "Class", show_row_dend = FALSE, show_column_dend = FALSE, col = cola_opt$color_set_2,
		column_order = mat_col_od, cluster_columns = FALSE,
		cluster_row_slices = FALSE,
		top_annotation = HeatmapAnnotation(Prob = membership_mat,
			Class = class_ids, col = c(list(Class = cola_opt$color_set_2), Prob = col_fun),
			show_annotation_name = FALSE, height = unit(2, "cm"), simple_anno_size_adjust = TRUE),
		show_column_names = FALSE,
		row_title = NULL, use_raster = TRUE,
		height = unit(1, "null")
	)

	consensus_mat = get_consensus(object, k = k)

	ht = ht %v% Heatmap(consensus_mat, col = colorRamp2(c(0, 1), c("white", "blue")),
		row_order = mat_col_od, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE, use_raster = TRUE,
		height = unit(1.2, "null"))
	draw(ht, main_heatmap = "Class",
		show_heatmap_legend = FALSE, show_annotation_legend = FALSE, column_title = qq("k = @{k}"))
}


p6 = grid.grabExpr(make_single_plot(res, k = 2), height = fig_height*1.6/2.6)
p7 = grid.grabExpr(make_single_plot(res, k = 3), height = fig_height*1.6/2.6)
p8 = grid.grabExpr(make_single_plot(res, k = 4), height = fig_height*1.6/2.6)
p9 = grid.grabExpr(make_single_plot(res, k = 5), height = fig_height*1.6/2.6)
p10 = grid.grabExpr(make_single_plot(res, k = 6), height = fig_height*1.6/2.6)
p11 = grid.grabExpr(make_single_plot(res, k = 7), height = fig_height*1.6/2.6)
p12 = grid.grabExpr(make_single_plot(res, k = 8), height = fig_height*1.6/2.6)


library(cowplot)

p = plot_grid( 
		plot_grid(p1, p2, p3, p4, p5, nrow = 1), 
		plot_grid(p6, p7, p8, p9, p10, p11, p12, nrow = 1), 
		ncol = 1, rel_heights = c(1, 1.6))

pdf("figure2.pdf", width = fig_width, height = fig_height)
print(p)
dev.off()
