
library(Seurat)
library(tidyverse)
library(infercnv)

gfp_label_1 <- as.data.frame(Idents(gfp_combined)) %>%
        rownames_to_column('Index')
gfp_label_2 <- as.data.frame(Idents(gfp_combined)) %>%
        rownames_to_column('Index')

colnames(gfp_label_1)[2] <- 'Label'
colnames(gfp_label_2)[2] <- 'Label'

gfp_label_1$Label <- as.character(gfp_label_1$Label)
gfp_label_2$Label <- as.character(gfp_label_2$Label)

gfp_label_1 <- gfp_label_1 %>%
        filter(str_detect(Index, "1$")) %>%
        filter(str_detect(Label, "(Basal|Luminal)")) %>%
        as_tibble()

gfp_label_2 <- gfp_label_2 %>%
        filter(str_detect(Index, "2$")) %>%
        filter(str_detect(Label, "(Basal|Luminal)")) %>%
        as_tibble()



# use immune

gfp_label_1 <- gfp_label_1 %>%
        filter(str_detect(Index, "1$")) %>%
        mutate(Label = ifelse(str_detect(Label, "(Basal|Luminal)"), Label, 'Normal')) %>%
        as_tibble()

gfp_label_2 <- gfp_label_2 %>%
        filter(str_detect(Index, "2$")) %>%
        mutate(Label = ifelse(str_detect(Label, "(Basal|Luminal)"), Label, 'Normal')) %>%
        as_tibble()

# output

write_tsv(gfp_label_1, "./data/gfp_label_1.tsv", col_names = F)
write_tsv(gfp_label_2, "./data/gfp_label_2.tsv", col_names = F)

gfp_object_1 <- gfp_combined@assays$RNA@counts
gfp_object_1 <- gfp_object_1[,gfp_label_1$Index]

saveRDS(gfp_object_1, "./data/gfp_object_1.rds")

gfp_object_2 <- gfp_combined@assays$RNA@counts
gfp_object_2 <- gfp_object_2[,gfp_label_2$Index]

saveRDS(gfp_object_2, "./data/gfp_object_2.rds")

# create the infercnv object
infercnv_gfp_1 <- CreateInfercnvObject(raw_counts_matrix = "./data/gfp_object_1.rds",
                                    annotations_file="./data/gfp_label_1.tsv",
                                    delim = "\t",
                                    gene_order_file="./data/mouse_gene_ordering.tsv",
                                    ref_group_names = "Normal")

infercnv_gfp_2 <- CreateInfercnvObject(raw_counts_matrix = "./data/gfp_object_2.rds",
                                       annotations_file="./data/gfp_label_2.tsv",
                                       delim = "\t",
                                       gene_order_file="./data/mouse_gene_ordering.tsv",
                                       ref_group_names = "Normal")

infercnv_gfp_1_new <- CreateInfercnvObject(raw_counts_matrix = "./data/gfp_object_1.rds",
                                       annotations_file="./data/gfp_label_1.tsv",
                                       delim = "\t",
                                       gene_order_file="./data/mouse_gene_ordering.tsv",
                                       ref_group_names = NULL)

infercnv_gfp_2_new <- CreateInfercnvObject(raw_counts_matrix = "./data/gfp_object_2.rds",
                                       annotations_file="./data/gfp_label_2.tsv",
                                       delim = "\t",
                                       gene_order_file="./data/mouse_gene_ordering.tsv",
                                       ref_group_names = NULL)

#030721

infercnv_gfp_1 <- CreateInfercnvObject(raw_counts_matrix = "./data/gfp_object_1.rds",
                                           annotations_file="./data/gfp_label_1.tsv",
                                           delim = "\t",
                                           gene_order_file="./data/mouse_gene_ordering.tsv",
                                           ref_group_names = 'Luminal 1')

infercnv_gfp_2 <- CreateInfercnvObject(raw_counts_matrix = "./data/gfp_object_2.rds",
                                           annotations_file="./data/gfp_label_2.tsv",
                                           delim = "\t",
                                           gene_order_file="./data/mouse_gene_ordering.tsv",
                                           ref_group_names = 'Luminal 1')

# perform infercnv operations to reveal cnv signal
infercnv_gfp_1_new <- infercnv::run(infercnv_gfp_1_new,
                                cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir = "./data/InferCNV/GFP_1_new",  # dir is auto-created for storing outputs
                                cluster_by_groups = T,   # cluster
                                denoise = T,
                                HMM = T)

infercnv_gfp_2_new <- infercnv::run(infercnv_gfp_2_new,
                                cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir = "./data/InferCNV/GFP_2_new",  # dir is auto-created for storing outputs
                                cluster_by_groups = T,   # cluster
                                denoise = T,
                                HMM = T)

infercnv_gfp_1 <- infercnv::run(infercnv_gfp_1,
                                    cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir = "./data/InferCNV/GFP_1",  # dir is auto-created for storing outputs
                                    cluster_by_groups = T,   # cluster
                                    denoise = T,
                                    HMM = T)

infercnv_gfp_2 <- infercnv::run(infercnv_gfp_2,
                                    cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                    out_dir = "./data/InferCNV/GFP_2",  # dir is auto-created for storing outputs
                                    cluster_by_groups = T,   # cluster
                                    denoise = T,
                                    HMM = T)


# score

cnv_1 <- read.table("./data/InferCNV/GFP_1_new/infercnv.21_denoised.observations.txt", sep = ' ')
cnv_1_score <- c()
for (i in 1:ncol(cnv_1)) {
        cnv_1_score[i] <- sum(cnv_1[,i] >= 1.15) + sum(cnv_1[,i] <= 0.85)
}

gfp_label_1 <- tibble(Index = colnames(cnv_1), Score = cnv_1_score) %>%
        left_join(gfp_label_1, by = 'Index')

ggplot(gfp_label_1) +
        geom_boxplot(aes(x = Label, y = Score))


cnv_2 <- read.table("./data/InferCNV/GFP_2_new/infercnv.21_denoised.observations.txt", sep = ' ')
cnv_2_score <- c()
for (i in 1:ncol(cnv_2)) {
        cnv_2_score[i] <- sum(cnv_2[,i] >= 1.15) + sum(cnv_2[,i] <= 0.85)
}

gfp_label_2 <- tibble(Index = colnames(cnv_2), Score = cnv_2_score) %>%
        left_join(gfp_label_2, by = 'Index')

ggplot(gfp_label_2) +
        geom_boxplot(aes(x = Label, y = Score))


# aneuploidy

infercnv_gfp_1_new@gene_order %>%
        rownames_to_column('Gene')

cnv_1_score <- as.data.frame(cnv_1) %>%
        rownames_to_column('Gene') %>%
        pivot_longer(cols = !Gene, names_to = 'Index', values_to = 'CNV') %>%
        left_join(infercnv_gfp_1_new@gene_order %>%
                          rownames_to_column('Gene'), by = 'Gene') %>%
        left_join(gfp_label_1, by = 'Index')

cnv_1_score <- cnv_1_score %>%
        mutate(Ab = CNV >= 1.1 | CNV <= 0.9) %>%
        group_by(Index, chr, Label) %>%
        summarize(Score = sum(Ab),
                  n = n()) %>%
        mutate(Pct = Score/n,
               S_Pct = Pct > 0.7) %>%
        group_by(Label, Index) %>%
        summarise(Aneuploidy = sum(S_Pct))


cnv_2_score <- as.data.frame(cnv_2) %>%
        rownames_to_column('Gene') %>%
        pivot_longer(cols = !Gene, names_to = 'Index', values_to = 'CNV') %>%
        left_join(infercnv_gfp_2_new@gene_order %>%
                          rownames_to_column('Gene'), by = 'Gene') %>%
        left_join(gfp_label_2, by = 'Index')

cnv_2_score <- cnv_2_score %>%
        mutate(Ab = CNV >= 1.1 | CNV <= 0.9) %>%
        group_by(Index, chr, Label) %>%
        summarize(Score = sum(Ab),
                  n = n()) %>%
        mutate(Pct = Score/n,
               S_Pct = Pct > 0.7) %>%
        group_by(Label, Index) %>%
        summarise(Aneuploidy = sum(S_Pct))



library(ggpubr)

library(patchwork)

ggviolin(cnv_1_score, x = "Label", y = "Aneuploidy",
         fill = "Label", palette = 'jco', add = 'jitter', add.params = list(alpha = 0.2, size = 0.5)) + 
        stat_compare_means(comparisons = list(c("Basal 1", "Basal 2"),
                                              c("Basal 1", "Luminal 1"),
                                              c("Basal 1", "Luminal 2")),
                           label.y = c(3, 3.5, 4)) +
        ggtitle('Indolent')
        

ggviolin(cnv_2_score, x = "Label", y = "Aneuploidy",
         fill = "Label", palette = 'jco',  add = 'jitter', add.params = list(alpha = 0.2, size = 0.5)) + 
        stat_compare_means(comparisons = list(c("Basal 1", "Basal 2"),
                                              c("Basal 1", "Luminal 1"),
                                              c("Basal 1", "Luminal 2")),
                           label.y = c(7, 8, 9)) +
        ggtitle('Aggressive')

cnv_1_score %>%
        mutate(Group = 'Indolent') %>%
        rbind(cnv_2_score %>%
                      mutate(Group = 'Aggressive')) %>%
        ggviolin(x = "Group", y = "Aneuploidy", fill = 'Group', palette = 'Set2',
                 add = 'jitter', add.params = list(alpha = 0.2, size = 0.5)) +
        stat_compare_means()

write_tsv(as.data.frame(rbind(
cnv_1_score %>%
        group_by(Label, Aneuploidy) %>%
        summarise(n = n()) %>%
        mutate(Group = 'Indolent'),
cnv_2_score %>%
        group_by(Label, Aneuploidy) %>%
        summarise(n = n()) %>%
        mutate(Group = 'Aggressive'))), file = './tables/Aneuploidy_Stat.tsv')


# Wilcox

wilcox.test(cnv_1_score %>% filter(Label == 'Basal 1') %>% pull(Aneuploidy),
            cnv_1_score %>% filter(Label == 'Basal 2') %>% pull(Aneuploidy))

wilcox.test(cnv_1_score %>% filter(Label == 'Basal 1') %>% pull(Aneuploidy),
            cnv_1_score %>% filter(Label == 'Luminal 1') %>% pull(Aneuploidy))

wilcox.test(cnv_1_score %>% filter(Label == 'Basal 1') %>% pull(Aneuploidy),
            cnv_1_score %>% filter(Label == 'Luminal 2') %>% pull(Aneuploidy))


wilcox.test(cnv_2_score %>% filter(Label == 'Basal 1') %>% pull(Aneuploidy),
            cnv_2_score %>% filter(Label == 'Basal 2') %>% pull(Aneuploidy))

wilcox.test(cnv_2_score %>% filter(Label == 'Basal 1') %>% pull(Aneuploidy),
            cnv_2_score %>% filter(Label == 'Luminal 1') %>% pull(Aneuploidy))

wilcox.test(cnv_2_score %>% filter(Label == 'Basal 1') %>% pull(Aneuploidy),
            cnv_2_score %>% filter(Label == 'Luminal 2') %>% pull(Aneuploidy))





        



