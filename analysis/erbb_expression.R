
library(Seurat)
library(tidyverse)

kegg_mm$`ErbB signaling pathway`

gfp_combined@assays$RNA@data@Dimnames[[1]]

kegg_mm$`ErbB signaling pathway`[kegg_mm$`ErbB signaling pathway` %in% gfp_combined@assays$integrated@var.features]
erbb_genes <- kegg_mm$`ErbB signaling pathway`[kegg_mm$`ErbB signaling pathway` %in% gfp_combined@assays$RNA@data@Dimnames[[1]]]

erbb_exp <- t(as.matrix(gfp_combined@assays$RNA@data)[erbb_genes,])

erbb_exp <- as.data.frame(erbb_exp)

Idents(gfp_combined)
gfp_combined$orig.ident

erbb_exp <- erbb_exp %>%
        rownames_to_column(var = 'Index') %>%
        add_column(Group = gfp_combined$orig.ident, .before = 2) %>%
        add_column(Cluster = Idents(gfp_combined), .before = 3)

gene_wilcox <- function(df, cluster, gene) {
        
        if (cluster != 'All') {
                v1 <- df %>%
                        filter(Group == 'Early' & Cluster == cluster) %>%
                        `[[`(gene)
                
                v2 <- df %>%
                        filter(Group == 'Advanced' & Cluster == cluster) %>%
                        `[[`(gene)
        } else {
                v1 <- df %>%
                        filter(Group == 'Early') %>%
                        `[[`(gene)
                
                v2 <- df %>%
                        filter(Group == 'Advanced') %>%
                        `[[`(gene)
        }
        
        diff <- mean(v2) - mean(v1)
        pval <- wilcox.test(v1, v2)$p.value
        
        return(c(diff, pval))
        
}

gene_wilcox(erbb_exp, 'All', 'Erbb2')

erbb_res <- list()
erbb_res$Gene <- c()
erbb_res$Cluster <- c()
erbb_res$Difference <- c()
erbb_res$p.value <- c()

for (i in erbb_genes) {
        for (j in c('All', 'Basal 1', 'Basal 2', 'Luminal 1', 'Luminal 2')) {
                erbb_res$Gene <- c(erbb_res$Gene, i)
                erbb_res$Cluster <- c(erbb_res$Cluster, j)
                erbb_res$Difference <- c(erbb_res$Difference, gene_wilcox(erbb_exp, j ,i)[1])
                erbb_res$p.value <- c(erbb_res$p.value, gene_wilcox(erbb_exp, j ,i)[2])
        }
}

erbb_res <- cbind(erbb_res[[1]], erbb_res[[2]], erbb_res[[3]], erbb_res[[4]])
erbb_res <- as_tibble(erbb_res)
colnames(erbb_res) <- c('Gene','Cluster','Difference','p.value')
erbb_res$Difference <- as.numeric(erbb_res$Difference)
erbb_res$p.value <- as.numeric(erbb_res$p.value)
erbb_res$p.adjust <- p.adjust(erbb_res$p.value, method = 'BH')

write_csv(erbb_res, "./tables/erbb_expression.csv")

