

library(ComplexHeatmap)
library(RColorBrewer)
library(cola)
library(grid)
library(circlize)
library(matrixStats)
library(cowplot)
library(GetoptLong)

fig_width = 12
fig_height = 10

set.seed(54)
mean_diff1 = rnorm(100)

m1 = do.call(rbind, lapply(1:100, function(i) {
	c(rnorm(50, mean = mean_diff1[i]), rnorm(50, mean = -mean_diff1[i]))
}))

mean_diff2 = rnorm(100)/2
mean_diff2[order(abs(mean_diff1))] = mean_diff2[order(abs(mean_diff2), decreasing = TRUE)]

m2 = do.call(rbind, lapply(1:100, function(i) {
	c(rnorm(10, mean = mean_diff2[i]), rnorm(10, mean = -mean_diff2[i]))
}))

m = cbind(m1, m2)

group = rep(c("A1", "A2", "B1", "B2"), times = c(50, 50, 10, 10))
group_col = structure(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"), names = c("A1", "A2", "B1", "B2"))


p0 = grid.grabExpr(draw(Heatmap(m, name = "mat", 
	top_annotation = HeatmapAnnotation(Group = group, col = list(Group = group_col)),
	show_row_dend = FALSE, column_title = "Heatmap of the random dataset",
	row_dend_reorder = mean_diff1, column_dend_reorder = as.numeric(factor(group))), merge_legend = TRUE, padding = unit(c(5, 5, 5, 5), "mm")),
	width = fig_width/3, height = fig_height/2)

p1 = ~dimension_reduction(m, pch = 16, col = group_col[group], method = "PCA", main = "PCA on the random matrix", scale_rows = FALSE)



hive_plot = function() {
	lt = list(
		"complete" = rowSds(m),
		"group_A" = rowSds(m[, 1:100]),
		"group_B" = rowSds(m[, 101:120])
	)

	rg = sapply(lt, range)
	n = length(lt[[1]])

	pos = list()
	pos[[1]] = list(x = lt[[1]]*sqrt(3)/2, y = -lt[[1]]/2)
	pos[[2]] = list(x = rep(0, n), y = lt[[2]])
	pos[[3]] = list(x = -lt[[3]]*sqrt(3)/2, y = -lt[[3]]/2)

	rank_list = lapply(lt, rank)

	calc_bezier = function(p0, p2, sign = 1) {
		l = sqrt( (p0[1] - p2[1])^2 + (p0[2] - p2[2])^2 )
		sin_theta = (p0[1] - p2[1])/l
		cos_theta = (p2[2] - p0[2])/l
		h = l/2
		p1 = numeric(2)
		p1[1] = (p0[1] + p2[1])/2 + h*cos_theta*sign
		p1[2] = (p0[2] + p2[2])/2 + h*sin_theta*sign
		circlize:::quadratic.bezier(p0, p1, p2)
	}

	grid.newpage()
	pushViewport(viewport(xscale = c(-2, 3), 
		                  yscale = c(-2, 3.3),
		                  width = unit(1, "snpc")*0.95, height = unit(1, "snpc")*0.95))

	all_rank_diff = c(abs(rank_list[[1]] - rank_list[[2]]),
		          abs(rank_list[[1]] - rank_list[[3]]),
		          abs(rank_list[[2]] - rank_list[[3]]))

	rank_diff_fun = colorRamp2(seq(0, quantile(all_rank_diff, 0.9), length = 11), rev(brewer.pal(11, "Spectral")), transparency = 0.6)

	for(i in 1:n) {
		p0 = c(pos[[1]]$x[i], pos[[1]]$y[i])
		p2 = c(pos[[2]]$x[i], pos[[2]]$y[i])

		seg = calc_bezier(p0, p2)
		rank_diff = abs(rank_list[[1]][i] - rank_list[[2]][i])
		grid.lines(seg[, 1], seg[, 2], default.units = "native", gp = gpar(col = rank_diff_fun(rank_diff)))
	}

	for(i in 1:n) {
		p0 = c(pos[[3]]$x[i], pos[[3]]$y[i])
		p2 = c(pos[[2]]$x[i], pos[[2]]$y[i])

		seg = calc_bezier(p0, p2, sign = -1)
		rank_diff = abs(rank_list[[2]][i] - rank_list[[3]][i])
		grid.lines(seg[, 1], seg[, 2], default.units = "native", gp = gpar(col = rank_diff_fun(rank_diff)))
	}

	for(i in 1:n) {
		p0 = c(pos[[3]]$x[i], pos[[3]]$y[i])
		p2 = c(pos[[1]]$x[i], pos[[1]]$y[i])

		seg = calc_bezier(p0, p2)
		rank_diff = abs(rank_list[[1]][i] - rank_list[[3]][i])
		grid.lines(seg[, 1], seg[, 2], default.units = "native", gp = gpar(col = rank_diff_fun(rank_diff)))
	}


	rg_x = range(lt[[1]])*sqrt(3)/2
	rg_y = -range(lt[[1]])/2
	grid.lines(c(0, rg_x[2]*1.2), c(0, rg_y[2]*1.2), default.units = "native", gp = gpar(col = "#AAAAAA"), arrow = arrow(angle = 20, length = unit(0.1, "inch")))
	# grid.lines(rg_x, rg_y, default.units = "native", gp = gpar(lwd = 4, lineend = 1, col = "#444444"))
	grid.text(round(min(lt[[1]]), 2), rg_x[1], rg_y[1], default.units = "native", just = c(1.3, 0.5), gp = gpar(fontsize = 7))
	grid.text(round(max(lt[[1]]), 2), rg_x[2], rg_y[2], default.units = "native", just = c(-0.3, 0.5), gp = gpar(fontsize = 7))
	grid.text("Complete\nmatrix", unit(2.2, "native"), unit(-1.8, "native"), just = "top")

	l = lt[[1]] >= median(lt[[1]])
	grid.points(pos[[1]]$x[l], pos[[1]]$y[l], default.units  = "native", pch = 16, gp = gpar(col = "red"), size = unit(4, "pt"))

	rg_x = c(0, 0)
	rg_y = range(lt[[2]])
	grid.lines(c(0, rg_x[2]*1.2), c(0, rg_y[2]*1.2), default.units = "native", gp = gpar(col = "#AAAAAA"), arrow = arrow(angle = 20, length = unit(0.1, "inch")))
	# grid.lines(rg_x, rg_y, default.units = "native", gp = gpar(lwd = 4, lineend = 1, col = "#444444"))
	grid.text(round(min(lt[[2]]), 2), rg_x[1], rg_y[1], default.units = "native", just = c(0.5, 1.5), gp = gpar(fontsize = 7))
	grid.text(round(max(lt[[2]]), 2), rg_x[2], rg_y[2], default.units = "native", just = c(0.5, -0.5), gp = gpar(fontsize = 7))
	grid.text("Group A1/A2", unit(0.2, "native"), unit(3.1, "native"), just = "left")

	l = lt[[1]] >= median(lt[[1]])
	grid.points(pos[[2]]$x[l], pos[[2]]$y[l], default.units  = "native", pch = 16, gp = gpar(col = "red"), size = unit(4, "pt"))

	rg_x = -range(lt[[3]])*sqrt(3)/2
	rg_y = -range(lt[[3]])/2
	grid.lines(c(0, rg_x[2]*1.2), c(0, rg_y[2]*1.2), default.units = "native", gp = gpar(col = "#AAAAAA"), arrow = arrow(angle = 20, length = unit(0.1, "inch")))
	# grid.lines(rg_x, rg_y, default.units = "native", gp = gpar(lwd = 4, lineend = 1, col = "#444444"))
	grid.text(round(min(lt[[3]]), 2), rg_x[1], rg_y[1], default.units = "native", just = c(-0.3, 0.5), gp = gpar(fontsize = 7))
	grid.text(round(max(lt[[3]]), 2), rg_x[2], rg_y[2], default.units = "native", just = c(1.3, 0.5), gp = gpar(fontsize = 7))
	grid.text("Group B1/B2", unit(-1.5, "native"), unit(-1.5, "native"), just = "top")

	l = lt[[1]] >= median(lt[[1]])
	grid.points(pos[[3]]$x[l], pos[[3]]$y[l], default.units  = "native", pch = 16, gp = gpar(col = "red"), size = unit(4, "pt"))

	rank_diff_fun = colorRamp2(seq(0, quantile(all_rank_diff, 0.9), length = 11), rev(brewer.pal(11, "Spectral")))
	lgd = Legend(title ="Rank difference", col_fun = rank_diff_fun, direction = "horizontal", title_position = "lefttop")
	draw(lgd, x = unit(0, "npc"), y = unit(-0.1, "npc"), just = c("left", "bottom"))
}


plot_ecdf = function(object) {
	plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "Consensus value (x)", ylab = "P(X <= x)", main = "eCDF curve")
	for(i in seq_along(object@k)) {
		consensus_mat = get_consensus(object, k = object@k[i])
		f = ecdf(consensus_mat[lower.tri(consensus_mat)])
		x = seq(0, 1, length = 100)
		y = f(x)
		x = c(0, x)
		y = c(0, y)
		lines(x, y, col = i)
	}
	legend("bottomright", pch = 15, legend = paste0("k = ", object@k), col = seq_along(object@k))
}


# copied from cola
membership_heatmap = function(object, k, internal = FALSE, 
	anno = object@anno, anno_col = get_anno_col(object),
	show_column_names = FALSE, column_names_gp = gpar(fontsize = 8), ...) {

	if(missing(k)) stop_wrap("k needs to be provided.")

	class_df = get_classes(object, k)
	class_ids = class_df$class

	membership_mat = get_membership(object, k)
	col_fun = colorRamp2(c(0, 1), c("white", "red"))

	membership_each = get_membership(object, k, each = TRUE)
	membership_each = t(membership_each)
	mat_col_od = cola:::column_order_by_group(factor(class_ids, levels = sort(unique(class_ids))), membership_each)

	col = cola_opt$color_set_1[1:k]

	if(is.null(anno)) {
		bottom_anno = NULL
	} else {
		if(is.atomic(anno)) {
			anno_nm = deparse(substitute(anno))
			anno = data.frame(anno)
			colnames(anno) = anno_nm
			if(!is.null(anno_col)) {
				anno_col = list(anno_col)
				names(anno_col) = anno_nm
			}
		} else if(ncol(anno) == 1) {
			if(!is.null(anno_col)) {
				if(is.atomic(anno_col)) {
					anno_col = list(anno_col)
					names(anno_col) = colnames(anno)
				}
			}
		}

		if(is.null(anno_col)) {
			bottom_anno = HeatmapAnnotation(df = anno,
				show_annotation_name = !internal, annotation_name_side = "right")
		} else {
			bottom_anno = HeatmapAnnotation(df = anno, col = anno_col,
				show_annotation_name = !internal, annotation_name_side = "right")
		}
	}

	param = get_param(object, k, unique = FALSE)

	top_n_level = unique(param$top_n)
	suppressWarnings(n_row_col <- structure(brewer.pal(length(top_n_level), "Accent"), names = top_n_level))
		
	ht_opt$HEATMAP_LEGEND_PADDING = unit(-2, "mm")

	ht = Heatmap(membership_each, name = "Class", show_row_dend = FALSE, show_column_dend = FALSE, col = cola_opt$color_set_2,
		column_title = ifelse(internal, "", qq("Membership heatmap, k = @{k}")), 
		column_order = mat_col_od, cluster_columns = FALSE,
		row_split = factor(param$top_n, levels = sort(unique(param$top_n))),
		cluster_row_slices = FALSE,
		top_annotation = HeatmapAnnotation(Probability = membership_mat,
			Class = class_ids, col = c(list(Class = cola_opt$color_set_2), Probability = col_fun),
			show_annotation_name = !internal),
		bottom_annotation = bottom_anno,
		show_column_names = show_column_names, column_names_gp = column_names_gp,
		row_title = NULL
	)
	draw(ht, main_heatmap = "Class", 
		show_heatmap_legend = FALSE, show_annotation_legend = !internal, padding = unit(c(5, 5, 5, 5), "mm"))

}

p2 = grid.grabExpr(hive_plot())

res = consensus_partition(m, partition_method = "kmeans", top_value_method = "SD", scale_rows = FALSE,
	anno = data.frame(Group = group), anno_col = list(Group = group_col), top_n = 50)

p3 = ~plot_ecdf(res)
p4 = grid.grabExpr(membership_heatmap(res, k = 3), width = fig_width/3, height = fig_height/2)
p5 = grid.grabExpr(membership_heatmap(res, k = 4), width = fig_width/3, height = fig_height/2)

# rh = hierarchical_partition(m, partition_method = "kmeans", top_value_method = "SD", anno = group, anno_col = group_col, top_n = 50)
# collect_classes(rh)

pdf("figure1.pdf", width = fig_width, height = fig_height)
print(plot_grid(p0, p1, p3, p4, p5, p2, ncol = 3))
dev.off()

