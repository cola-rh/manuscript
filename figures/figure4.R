

fig_width = 16
fig_height = 12

library(cola)
library(GetoptLong)
library(ComplexHeatmap)
library(circlize)

library(HSMMSingleCell)
data(HSMM_expr_matrix)
data(HSMM_sample_sheet)

# `HSMM_expr_matrix` is a FPKM matrix
m = adjust_matrix(log10(HSMM_expr_matrix + 1), verbose = FALSE)
gt = readRDS("~/workspace/cola_hc/gene_type_gencode_v17.rds")
m = m[gt[rownames(m)] == "protein_coding", , drop = FALSE]


### normal consensus partitioning results and classifications
res = readRDS("HSMM_cola_consensus_partition.rds")
cola_2_groups = as.character(get_classes(res, k = 2)[, "class"])
cola_6_groups = as.character(get_classes(res, k = 6)[, "class"])

### perform hierarchical consensus partitioning
set.seed(1234)
if(file.exists("HSMM_cola_hierarchical_partition.rds")) {
	rh = readRDS("HSMM_cola_hierarchical_partition.rds")
} else {
	rh = hierarchical_partition(m, cores = 4)
	saveRDS(rh, file = "HSMM_cola_hierarchical_partition.rds")
}

cola_cl = data.frame(cola_2_groups = cola_2_groups, cola_6_groups = cola_6_groups)

### t-SNE plot
p1 = ~dimension_reduction(rh, method = "t-SNE", top_n = 1000, pc = 10)

### top rows overlap at each node
rl = rh@list
nodes = setdiff(all_nodes(rh), all_leaves(rh))
rl = rl[nodes]
all_top_value_list = lapply(rl, function(x) {
	x@row_index[order(x@top_value_list)[seq_len(max(x@top_n))]]
})
names(all_top_value_list) = paste0(names(all_top_value_list), "|", sapply(rl, function(x) x@top_value_method))

cm = make_comb_mat(all_top_value_list)
cm = cm[comb_size(cm) > 20]
p2 = grid.grabExpr(draw(UpSet(cm, column_title = "Overlap of top n features on each node")))

### signature heatmap
sig_tb = get_signatures(rh, plot = FALSE)

cl = get_classes(rh)
ht = Heatmap(t(scale(t(m[sig_tb$which_row, ]))), name = "Z-score", column_title = qq("@{length(unique(cl))} groups, @{nrow(sig_tb)} signatures under FDR < 0.05"),
	col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
	top_annotation = HeatmapAnnotation(Class = cl, col = list(Class = rh@subgroup_col)),
	bottom_annotation = HeatmapAnnotation(df = cola_cl, col = list(cola_2_groups = cola_opt$color_set_2, cola_6_groups = cola_opt$color_set_2),
		annotation_legend_param = list(cola_6_groups = list(ncol = 2))),
	show_row_names = FALSE, show_column_names = FALSE, show_row_dend = FALSE,
	cluster_columns = cola:::calc_dend(rh), 
	column_split = length(unique(cl)), row_km = 5, row_title = paste0("km", 1:5))
p3 = grid.grabExpr(ht <- draw(ht, merge_legend = TRUE, padding = unit(c(8, 2, 2, 2), "mm")))

### functional enrichment on row signatures
gene_list = lapply(row_order(ht), function(ind) rownames(m)[sig_tb$which_row][ind])
names(gene_list) = paste0("km", 1:5)

lt_GO = lapply(gene_list, function(x) functional_enrichment(x)$BP)

library(simplifyEnrichment)
p_go = grid.grabExpr(simplifyGOFromMultipleLists(lt_GO, word_cloud_grob_param = list(max_width = unit(6, "cm"))), 
	width = fig_width*0.6667, height = fig_height*0.5)


### stability of the classification
rh_list = readRDS("HSMM_cola_hierarchical_partition_rh_list.rds")
clt = lapply(rh_list, function(x) relabel_class(x, get_classes(rh), return_map = FALSE))
clm = do.call(cbind, clt)

library(clue)
p_stability = grid.grabExpr(draw(Heatmap(clm, name = "Class", row_order = do.call(order, clt), col = rh@subgroup_col,
	cluster_columns = hclust(cl_dissimilarity(cl_ensemble(list = lapply(clt, as.cl_hard_partition)))),
	column_title = "Stability of classification\nfrom 20 runs",
	show_heatmap_legend = FALSE,
	show_column_dend = FALSE), padding = unit(c(8, 2, 2, 2), "mm")))

### putting all figures together
library(cowplot)

pdf("figure4.pdf", width = fig_width, height = fig_height)
p = plot_grid( plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1.6)), 
	           plot_grid(p3, p_go, p_stability, nrow = 1, rel_widths = c(1, 1.6, 0.4)), nrow = 2, rel_heights = c(1, 1.2))
print(p)
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
