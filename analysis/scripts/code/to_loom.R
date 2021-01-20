# Convert to loom file to try scvi
library(Seurat)
library(loomR)

seu <- readRDS("./output/BronnerLab_Sha/seu_concat.rds")
seu_loom <- as.loom(seu, assay = "RNA", filename = "./output/BronnerLab_Sha/seu_concat.loom")
seu_loom <- connect("./output/BronnerLab_Sha/seu_concat.loom", "r+")
seu_loom$add.col.attribute(list(BatchID = as.integer(factor(seu$orig.ident, 
                                                            levels = c("premigratory.crest", 
                                                                       "migratory.crest", 
                                                                       "circumpharyngeal.crest",
                                                                       "outflowtract.crest"))) - 1L))
seu_loom$close_all()
