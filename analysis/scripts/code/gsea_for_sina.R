library(BUSpaRse)
library(tidyverse)
library(DropletUtils)
library(Matrix)
library(Seurat)
library(loomR)
methodDE <- function(seu1, seu2, ...) {
  # Make cell names unique
  seu1 <- RenameCells(seu1, add.cell.id = "a")
  seu2 <- RenameCells(seu2, add.cell.id = "b")
  seu <- merge(seu1, seu2)
  Idents(seu) <- "orig.ident"
  markers <- FindMarkers(seu, ident.1 = unique(as.character(seu1$orig.ident)), ...) %>% 
    clean_markers()
}

# need to transpose matrix
my_read_loom <- function(loom_path, cells = "obs_names", features = "var_names") {
  loom_path <- normalizePath(loom_path)
  x <- connect(loom_path)
  counts <- x$get.sparse(
    dataset = 'matrix',
    feature.names = features,
    cell.names = cells
  )
  seu <- CreateSeuratObject(t(counts)) %>% 
    NormalizeData() %>% 
    ScaleData()
}
Af <- my_read_loom("./velocity_gsea_for_lambda/Af.loom")
Bf <- my_read_loom("./velocity_gsea_for_lambda/Bf.loom")
Af$orig.ident <- "Af"
Bf$orig.ident <- "Bf"

gns <- read_csv("./data/Ensembl/gns_Homo_sapiens.csv")
kegg_kgsets <- readRDS("./data/kegg/kegg_Homo_sapiens.rds")
kegg_df <- read_csv("./data/kegg/kegg_df_Homo_sapiens.csv")

genes_use <- rownames(GetAssayData(Af, slot = "data"))
cells_use <- Cells(Af)

entrez_inter <- gns$entrezgene[gns$gene %in% genes_use]
entrez_inter <- entrez_inter[!is.na(entrez_inter)]

method_markers <- methodDE(Af, Bf, test.use = "LR", 
                           latent.vars = "nCount_RNA") %>% 
  filter(p_val_adj < 0.05)

method_markers <- method_markers %>% 
  left_join(gns, by = "gene")
method_gsea <- egsea_ora(method_markers$gene, "Homo sapiens", gns, entrez_inter, kegg_kgsets,
                         logFC = method_markers$avg_logFC) %>% 
  egsea_results_df()
out <- method_markers %>% 
  left_join(kegg_df, by = "entrezgene") %>% 
  inner_join(method_gsea, by = "gene_set") %>% 
  filter(p.adj < 0.05) %>% 
  mutate(change = cut(avg_logFC, -20:20))

write_csv(out, "./velocity_gsea_for_lambda/out.csv")

p <- ggplot(out, aes(fct_reorder(gene_set, gene_set, length, .desc = TRUE), 
                fill = fct_rev(change))) +
  geom_bar() +
  scale_fill_viridis_d(option = "E", name = "logFC", direction = -1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(y = "Number of genes", x = "gene set")
ggsave("./velocity_gsea_for_lambda/out.png", 
       p, device = "png", width = 9, height = 6)
