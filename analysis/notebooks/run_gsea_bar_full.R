# Make the gsea bar plot for all datasets
source("./gsea_bar_full.R")

dataset_meta <- tribble(~name_regex,	~species, ~mapping,
                        "mouse-EMTAB7320_v2",	"Mus musculus", "org.Mm.eg.db",
                        "mouse-heart1k_v2",	"Mus musculus", "org.Mm.eg.db",
                        "mouse-heart1k_v3",	"Mus musculus", "org.Mm.eg.db",
                        "human_mouse-hgmm10k_v3",	c("Homo sapiens", "Mus musculus"), c("org.Hs.eg.db", "org.Mm.eg.db"),
                        "human_mouse-hgmm1k_v2",	c("Homo sapiens", "Mus musculus"), c("org.Hs.eg.db", "org.Mm.eg.db"),
                        "human_mouse-hgmm1k_v3",	c("Homo sapiens", "Mus musculus"), c("org.Hs.eg.db", "org.Mm.eg.db"),
                        "mouse-neuron10k_v3",	"Mus musculus", "org.Mm.eg.db",
                        "human-pbmc10k_v3",	"Homo sapiens", "org.Hs.eg.db",
                        "human-pbmc1k_v3",	"Homo sapiens", "org.Hs.eg.db",
                        "zebrafish-SRR6956073_v2", "Danio rerio", "org.Dr.eg.db",
                        "mouse-SRR6998058_v2", "Mus musculus", "org.Mm.eg.db",
                        "rat-SRR7299563_v2", "Rattus norvegicus", "org.Rn.eg.db",
                        #"SRR7692543_v2", "Macaca fascicularis",
                        "mouse-SRR8206317_v2", "Mus musculus", "org.Mm.eg.db",
                        "arabidopsis-SRR8257100_v2", "Arabidopsis thaliana", "org.At.tair.db",
                        "human-SRR8327928_v2", "Homo sapiens", "org.Hs.eg.db",
                        "fly-SRR8513910_v2", "Drosophila melanogaster", "org.Dm.eg.db",
                        "human-SRR8524760_v2", "Homo sapiens", "org.Hs.eg.db",
                        "mouse-SRR8599150_v2", "Mus musculus", "org.Mm.eg.db",
                        "worm-SRR8611943_v2", "Caenorhabditis elegans", "org.Ce.eg.db",
                        "mouse-SRR8639063_v2", "Mus musculus", "org.Mm.eg.db")
dataset_meta <- dataset_meta %>% 
  mutate(name_use = str_remove(name_regex, "_\\?"))

# Now begins the run
dataset_meta <- dataset_meta %>% 
  mutate(plot = pmap(list(species, name_regex, mapping, name_use), ~gsea_bar_full(..1, ..2, ..3, ..4)))

# # Change all labels for qq plots into full, not truncated names
walk2(dataset_meta$name_use, dataset_meta$species,
      function(x, y) {
        qq_df <- read_csv(paste0("../../data/gsea_qq/", x, ".csv"))
        go_def <- map_dfr(y, function(s) {
          read_csv(paste0("../../reference/GO/go_def", str_replace(s, " ", "_"), ".csv"))
        })
        if (any(!is.na(qq_df$label))) {
          qq_df <- qq_df %>% 
            left_join(go_def) %>% 
            mutate(label = case_when(
              !is.na(label) & !is.na(GO_full_name) & label != GO_full_name ~ GO_full_name,
              TRUE ~ label
            ))
          write_csv(qq_df, paste0("../../data/gsea_qq/", x, ".csv"))
        }
      })

# # Change all labels for qq plots into full, not truncated names
walk2(dataset_meta$name_use, dataset_meta$species,
      function(x, y) {
        bar_df <- read_csv(paste0("../../data/gsea_bar/", x, ".csv"))
        go_def <- map_dfr(y, function(s) {
          read_csv(paste0("../../reference/GO/go_def", str_replace(s, " ", "_"), ".csv"))
        })
        bar_df <- bar_df %>% 
          left_join(go_def) %>% 
          mutate(label = case_when(
            (!is.na(label) & !is.na(GO_full_name) & label != GO_full_name) ~ GO_full_name,
            TRUE ~ label
          )) %>% 
          dplyr::select(p_val:description, change, gene_set = label)
        write_csv(bar_df, paste0("../../data/gsea_bar_name_fix/", x, ".csv"))
      })
