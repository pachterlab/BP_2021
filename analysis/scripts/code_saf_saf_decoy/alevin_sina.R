# Convert alevin output to sparse matrix for Sina
library(tidyverse)
source("./code/read_counnt_output.R")

alevin_outs <- list.dirs(path = "./to_convert",recursive = TRUE)
alevin_outs <- alevin_outs[str_detect(alevin_outs, "alevin$")]
save_paths <- str_replace(list.dirs(path = "./to_convert", recursive = FALSE), 
                          "to_convert", "alevin_converted")

mats <- map2(alevin_outs, save_paths, ~ readAlevin_sparse(.x, save_mtx = TRUE, save_path = .y))
