


library(cola)
library(GetoptLong)
library(ComplexHeatmap)
library(circlize)

tumor_type = read.table("~/workspace/cola_hc/GSE90496_tumor_type.csv", sep = ";", header = TRUE, comment.char = "")
tumor_col = structure(unique(tumor_type$color), names = unique(tumor_type$tumor_type))
tumor_type = structure(tumor_type$tumor_type, names = tumor_type$subtype)

lt = readRDS("~/workspace/cola_hc/GSE90496_islands_processed.rds")
mat = as.matrix(lt$mat)
anno = lt$anno
anno$subtype2 = gsub("\\s+", " ", anno$subtype2)
anno$subtype2 = gsub("\\s+$", "", anno$subtype2)

anno = data.frame(meth_class = anno$subtype2, tumor_type = tumor_type[anno$subtype2])

anno_col = list(
    tumor_type = tumor_col
)

mat = adjust_matrix(mat)

### t-SNE plot
# s = rowSds(mat)
# dimension_reduction(mat[order(-s)[1:5000], ], control = list(theta = 0, pca = FALSE), 
#     col = col[anno$subtype2], method = "t-SNE")

# mat2 = mat[order(-s)[1:30000], ]
# s2 = ATC(mat2, cores = 4)
# dimension_reduction(mat2[order(-s2)[1:5000], ], control = list(theta = 0, pca = FALSE), 
#     col = col[anno$subtype2], method = "t-SNE")
# s2 = ATC(mat2, k_neighbours = 100, cores = 4)
# dimension_reduction(mat2[order(-s2)[1:5000], ], control = list(theta = 0, pca = FALSE), 
#     col = col[anno$subtype2], method = "t-SNE")

### perform hierarchical partitioning 
set.seed(123)
if(!file.exists("~/workspace/cola_hc/GSE90496_cola_rh.rds")) {
    rh = hierarchical_partition(mat, cores = 4, 
    	top_value_method = c("SD", "ATC"), max_k = 8, 
    	partition_method = c("kmeans", "skmeans"),
        scale_rows = FALSE, anno = anno, top_n = 1000, 
        subset = 500, group_diff = 0.25, min_n_signatures = 1000,
        filter_fun = function(mat) {
            s = rowSds(mat)
            order(-s)[1:30000]
        })

    saveRDS(rh, file = "~/workspace/cola_hc/GSE90496_cola_rh.rds")
} else {
    rh = readRDS("~/workspace/cola_hc/GSE90496_cola_rh.rds")
}

# ## do a subset
# res = list()
# for(i in 1:10) {
#     res[[i]] = consensus_partition_by_down_sampling(mat2, top_value_method = "ATC", partition_method = c("skmeans"), max_k = 10,
#         scale_rows = FALSE, anno = anno, top_n = 1000, cores = 4, subset = 500)
# }


### correspondance to original classification
tb = anno
tb$cola_class = get_classes(rh)

overlap_coefficient = function(x, y) {
    le1 = unique(x)
    le2 = unique(y)

    om = matrix(nrow = length(le1), ncol = length(le2))
    dimnames(om) = list(le1, le2)
    for(a in le1) {
        for(b in le2) {
            om[a, b] = sum(x == a & y == b)/min(sum(x == a), sum(y == b))
        }
    }
    om
}


get_column_order = function(om, row_order = order(-rowSums(om))) {
    scoreCol = function(x) {
        score = 0
        for(i in 1:length(x)) {
            if(x[i] > 0) {
                score = score + 2^(length(x)-i*1/(x[i] > 0.4))
            }
        }
        return(score)
    }
    scores = apply(om[row_order, ,drop = FALSE], 2, scoreCol)
    order(scores, decreasing=TRUE)
}

reorder_om = function(om) {
    row_order = order(-rowSums(om))
    column_order = get_column_order(om, row_order)
    om = om[row_order, column_order]

    # reorder rows
    for(i in 1:10) {
        om = t(om)
        od = get_column_order(om, 1:nrow(om))
        # reorder columns
        om = t(om[, od])
        od = get_column_order(om, 1:nrow(om))
        om = om[, od]
        print(dim(om))
    }
    om
}

make_correspondance = function(tb, which = "tumor_type", plot = "heatmap") {
    om = overlap_coefficient(tb[[which]], tb$cola_class)
    om = reorder_om(om)
    
    if(plot == "heatmap") {
        draw(Heatmap(om, name = "Overlap coeffcient", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
            show_row_names = FALSE, show_column_names = FALSE,
            column_title = "cola HCP classification", column_title_side = "bottom",
            row_title = ifelse(which == "tumor_type", "Tumor types", "Methylation classes"), row_title_side = "left",
            cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = which == "meth_class",
            heatmap_legend_param = list(direction = "horizontal", at = seq(0, 1, by = 0.2), title_position = "lefttop")),
        heatmap_legend_side = "bottom")
    } else {

        df = readRDS(system.file("extdata", "CNS_tumour_classification.rds", package = "spiralize"))
        meth_col = structure(names = unique(df$meth_class), unique(df$meth_col))
        tumor_col = structure(names = unique(as.vector(df$tumor_type)), unique(df$tumor_col))

        if(which == "tumor_type") {
            col = tumor_col
        } else {
            col = meth_col
        }

        foo = table(tb[[which]], tb$cola_class)
        foo = foo[rownames(om), colnames(om)]
        foo = as.matrix(foo)
        w1 = rowSums(foo)
        w2 = colSums(foo)

        cum1 = cumsum(w1)
        cum2 = cumsum(w2)

        grid.newpage()
        pushViewport(viewport(yscale = c(0, sum(w1)), height = unit(1, "npc") - unit(4, "mm")))
        for(i in seq_len(nrow(foo))) {
            for(j in seq_len(ncol(foo))) {
                if(foo[i, j] == 0) next
                x11 = cum1[i] - w1[i] + sum(foo[i, 1:j]) - foo[i, j]
                x12 = x11 + foo[i, j]

                x21 = cum2[j] - w2[j] + sum(foo[1:i, j]) - foo[i, j]
                x22 = x21 + foo[i, j]

                pt1 = circlize:::get_bezier_points(x11, 0.2, x21, 0.8, xlim = c(0, sum(w1)), ylim = c(0, 1))
                pt2 = circlize:::get_bezier_points(x12, 0.2, x22, 0.8, xlim = c(0, sum(w1)), ylim = c(0, 1))

                grid.polygon(
                    c(0.2, 0.2, pt2[, 2], rev(pt1[, 2]), 0.8), default.unit = "native",
                    c(x11, x12, pt2[, 1], rev(pt1[, 1]), x11),
                    gp = gpar(fill = add_transparency(col[rownames(foo)][i], 0.25), col = NA)
                )
            }
        }

        grid.rect(y = cum1, x = unit(0.2, "native") - ifelse(seq_along(cm1) %% 2 == 1, 1, 0)*unit(3, "mm"), height = w1, width = unit(3, "mm"), just = c("right", "top"), default.units = "native",
            gp = gpar(fill = col[rownames(foo)], col = NA))
        grid.text(ifelse(which == "tumor_type", "Tumor types", "Methylation classes"), x = unit(0.2, "native") - unit(10, "mm"), y = 0.5, rot = 90, just = "bottom")
    
        grid.rect(y = cum2, x = unit(0.8) + ifelse(seq_along(cm1) %% 2 == 1, 1, 0)*unit(3, "mm"), height = w2, width = unit(3, "mm"), just = c("left", "top"), default.units = "native",
            gp = gpar(fill = rh@subgroup_col[colnames(foo)], col = NA))
        grid.text("cola HCP classification", x = unit(0.8, "native") + unit(10, "mm"), y = 0.5, rot = 90, just = "top")
    }
}


tb = anno
tb$cola_class = get_classes(rh)

p1 = grid.grabExpr(make_correspondance(tb, "tumor_type", ""))
p2 = grid.grabExpr(make_correspondance(tb, "meth_class", ""))
p3 = grid.grabExpr(make_correspondance(tb, "tumor_type", "heatmap"))
p4 = grid.grabExpr(make_correspondance(tb, "meth_class", "heatmap"))

library(cowplot)

pdf("figure6.pdf", width = 12, height = 8)
plot_grid(p1, p2, plot_grid(p3, p4, ncol = 1, rel_heights = c(1, 3)), nrow = 1, rel_widths = c(1, 1, 3))
dev.off()


library(matrixStats)
overall_classification_agreement = function(x, y) {
    m = overlap_coefficient(x, y)

    tb1 = table(x)
    v = rowMaxs(m)
    names(v) = rownames(m)
    v = v[names(tb1)]
    sum(v*tb1)/sum(tb1)
}

overall_classification_agreement(tb$cola_class, tb$tumor_type)
overall_classification_agreement(tb$cola_class, tb$meth_class)


