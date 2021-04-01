
# Integrative analysis
DefaultAssay(object = cd45_combined) <- "integrated"

cd45_combined <- subset(cd45_combined, cells = setdiff(names(cd45_combined@active.ident), c(doublets[[2]],doublets[[4]])))

cd45_combined <- ScaleData(object = cd45_combined,
                          vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"),
                          verbose = T)

cd45_combined <- analyze_merged(cd45_combined, group.levels = stages, 
                                verbose = T, npcs = 50, dims = 1:30, nnei = 70, k.param = 70, min.dist = 0.6, spread = 1.5, resolution = 0.3)

# Visualization

svglite(file = "./figures/cd45_merge.svg")
plot_merge(cd45_combined)
dev.off()

plot_cluster(cd45_combined, label = F)

plot_split(cd45_combined, colors = get_colors(seq(levels(cd45_combined$seurat_clusters)), pal = "Set3"))

svglite(file = "./figures/cd45_markers.svg", width = 15, height = 8)
plot_features(cd45_combined, features = c("Cd3d","Gzma","S100a8","Ebf1","Lyz2","Siglech","Iglv1","Col1a1","Plbd1","Ptprc"), ncol = 5)
dev.off()

svglite(file = "./figures/cd45_markers_2.svg", width = 15, height = 8)
plot_features(cd45_combined, features = c("Cd4","Cd8a", "Cd8b1"), ncol = 3)
dev.off()



FeaturePlot(cd45_combined, features = c("Il17a","Il17b","Il17c","Il17d","Il25","Il17f"), split.by = "group", 
            cols = c("lightgrey", "red"))

FeaturePlot(cd45_combined, features = c("Il17ra","Il17rb","Il17rc","Il17re"), split.by = "group", 
            cols = c("lightgrey", "red"))

plot_features(cd45_combined, features = c("Cd3d","Cd4", "Cd8a", "Cd8b1"), ncol = 2)
plot_features(cd45_combined, features = c("Pdcd1","Ctla4"), ncol = 2)
plot_features(cd45_combined, features = c("Havcr2","Lag3","Ido1","Tnfrsf4","Vsir","Icos","Gzmb","Foxp3"), ncol = 4)


# Identify cell markers

cd45_markers <- FindAllMarkers(cd45_combined, only.pos = T, logfc.threshold = 0.1)

# Heatmap

plot_heatmap(cd45_combined, cd45_markers, 5, cluster_pal = c("Set3"))

# Relabeling

cd45_labels <- c("T Cell 1",
                 "T Cell 2",
                 "NK Cell", 
                 "Neutrophil", 
                 "B Cell",
                 "Macrophage",
                 "pDC",
                 "Plasma Cell",
                 "cDC",
                 "T Cell 3",
                 "T Cell 4"
                 )

cd45_levels <- cd45_labels[c(1,2,10,11,3,5,8,6,4,7,9)]

cd45_combined <- rename_cluster(cd45_combined, cd45_labels)

# New visualization

cd45_index <- c(1:4,12,10,9,8,6,7,5)

svglite(file = "./figures/cd45_cluster.svg")
plot_cluster(cd45_combined, label = F, levels = cd45_levels, self_set_color = T, self_colors = get_colors(cd45_index))
dev.off()

# Statistics

svglite(file = "./figures/cd45_group_count.svg")
plot_stat(cd45_combined, "group_count", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 3)
dev.off()

svglite(file = "./figures/cd45_cluster_count.svg")
plot_stat(cd45_combined, "cluster_count", group_levels = stages, cluster_levels = cd45_levels,
          self_set_color = T, self_colors = get_colors(cd45_index))
dev.off()

svglite(file = "./figures/cd45_prop_fill.svg")
plot_stat(cd45_combined, "prop_fill", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 3,
          self_set_color = T, self_colors = get_colors(cd45_index))
dev.off()

svglite(file = "./figures/cd45_prop_diverge.svg")
plot_stat(cd45_combined, "prop_diverge", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 0.8)
dev.off()

plot_stat(cd45_combined, "prop_multi", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 1.5)


# DE analysis

cd45_diff <- find_diff_genes(dataset = cd45_combined, clusters = cd45_levels, groups = stages, logfc = 0)

cd45_hallmark_gsea <- test_GSEA(cd45_diff, cd45_levels, pathways.hallmark)

svglite(file = "./figures/cd45_gsea.svg", width = 12, height = 12)
plot_GSEA(cd45_hallmark_gsea, p_cutoff = 0.1, levels = cd45_levels)
dev.off()

cd45_kegg_gsea <- test_GSEA(cd45_diff, cd45_levels, pathways.kegg)

tiff(file = "./figures/cd45_gsea_kegg.tiff", width = 2400, height = 14400, res = 200)
plot_GSEA(cd45_kegg_gsea, p_cutoff = 0.1, levels = cd45_levels)
dev.off()

# Violin plot

plot_measure(cd45_combined, features = c("Il17a","Il17b","Il17c","Il17d","Il25","Il17f"), plot_type = "cluster_group", show = "box")

plot_measure(cd45_combined, features = c("Il17ra","Il17rb","Il17rc","Il17re"), plot_type = "cluster_group", show = "box")

plot_measure(cd45_combined, features = c("Pdcd1","Cd274"), plot_type = "cluster_group", show = "box")

plot_measure(cd45_combined, 
             features = c("Cxcl9", "Mmp9", "Mmp3", "Crp", "Il1rn", "Ccl6", "Mpo", "Cd40", "Icam1", "Vcam1", "Il33", "Eng", "Igfbp5"), 
             plot_type = "cluster_group", show = "box")

plot_measure(cd45_combined, 
             features = c("Ccl11", "Ccl21a", "Ccl21b", "Ccl21c", "Ccl5", "Ccl22", "F3", "Postn", "Igfbp3", "Il1a", "Serpine1", "Vegfa"), 
             plot_type = "cluster_group", show = "box")

plot_measure(cd45_combined, 
             features = c("Cxcl5", "Prl2c2", "Ldlr", "Rarres2", "Ccl12", "Adipoq", "Igfbp2", "Il4", "Apcs", "Retn", "Gas6", "Serpinf1"), 
             plot_type = "cluster_group", show = "box")


# Feature plot

plot_measure_cluster(cd45_combined, features = c("Il17a","Il17b","Il17c","Il17d","Il25","Il17f"), plot_type = "group")

plot_measure_cluster(cd45_combined, features = c("Il17ra","Il17rb","Il17rc","Il17re"), plot_type = "group")


plot_measure_cluster(cd45_combined, features = c("Cd3d","Cd4", "Cd8a", "Cd8b1", "Pdcd1","Ctla4"))

plot_measure_cluster(cd45_combined, features = c("Havcr2","Lag3","Ido1","Tnfrsf4","Vsir","Icos","Gzmb","Foxp3"))

plot_measure_cluster(cd45_combined, features = c("Pdcd1","Cd274"), plot_type = "group")


svg("./figures/Il17a_f_ra_rc_cd45_merged.svg", width = 17, height = 17)
plot_measure_dim_2(cd45_combined, c("Il17a","Il17f","Il17ra","Il17rc"))
dev.off()

svg("./figures/Il17a_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Il17a"), split_by = "group")
dev.off()

svg("./figures/Il17f_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Il17f"), split_by = "group")
dev.off()

svg("./figures/Il17ra_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Il17ra"), split_by = "group")
dev.off()

svg("./figures/Il17rc_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Il17rc"), split_by = "group")
dev.off()

svg("./figures/Pdcd1_Cd274_cd45_merged.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Pdcd1","Cd274"))
dev.off()

svg("./figures/Pdcd1_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Pdcd1"), split_by = "group")
dev.off()

svg("./figures/Cd274_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Cd274"), split_by = "group")
dev.off()

svg("./figures/Cd3d_Cd4_Cd8a_Cd8b1_Pdcd1_Ctla4_cd45_merged.svg", width = 25.5, height = 17)
plot_measure_dim_2(cd45_combined, c("Cd3d","Cd4","Cd8a","Cd8b1","Pdcd1","Ctla4"))
dev.off()

svg("./figures/Cd3d_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Cd3d"), split_by = "group")
dev.off()

svg("./figures/Cd4_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Cd4"), split_by = "group")
dev.off()

svg("./figures/Cd8a_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Cd8a"), split_by = "group")
dev.off()

svg("./figures/Cd8b1_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Cd8b1"), split_by = "group")
dev.off()

svg("./figures/Ctla4_cd45_split.svg", width = 17, height = 8.5)
plot_measure_dim_2(cd45_combined, c("Ctla4"), split_by = "group")
dev.off()

# Test

test2 <- cd45_combined

test2$celltype <- Idents(test2)

DefaultAssay(test2) <- "RNA"

svglite(file = "./figures/IL17_1.svg", width = 12, height = 12)
plots <- VlnPlot(test2, features = c("S100a9", "S100a8", "Fosb"), split.by = "group", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

svglite(file = "./figures/IL17_2.svg", width = 12, height = 12)
plots <- VlnPlot(test2, features = c("Il17a", "Il17ra"), split.by = "group", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()




