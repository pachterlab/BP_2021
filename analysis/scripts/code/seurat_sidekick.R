library(lmtest)
library(pbapply)
library(future.apply)
library(future)
library(igraph)
library(leidenbase)

GroupSingletons2 <- function(ids, SNN, min_cluster_size, group.singletons = TRUE, 
                             verbose = TRUE) {
  # identify singletons
  singletons <- names(x = which(x = table(ids) < min_cluster_size))
  singletons <- intersect(x = unique(x = ids), singletons)
  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }
  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to
  cluster_names <- as.character(x = unique(x = ids))
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  new.ids <- ids
  for (i in singletons) {
    i.cells <- names(which(ids == i))
    for (j in cluster_names) {
      j.cells <- names(which(ids == j))
      subSNN <- SNN[i.cells, j.cells]
      set.seed(1) # to match previous behavior, random seed being set in WhichCells
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
      } else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
    ids[i.cells] <- closest_cluster
  }
  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }
  return(ids)
}

#' My version of FindClusters
#' 
#' Leiden clustering using the R package `leidenbase`. I have trouble using
#' anything that depends on `reticulate` on servers. 
#' @inheritParams Seurat::FindClusters
FindClusters2 <- function(object, graph.name = NULL, resolution = 0.8, 
                          n.iter = 10, partition.type = "RBConfigurationVertexPartition",
                          min_cluster_size = 5, group.singletons = TRUE,
                          random.seed = 1, verbose = TRUE) {
  graph.name <- graph.name %||% paste0(DefaultAssay(object = object), "_snn")
  if (!graph.name %in% names(x = object)) {
    stop("Provided graph.name not present in Seurat object")
  }
  if (!is(object = object[[graph.name]], class2 = "Graph")) {
    stop("Provided graph.name does not correspond to a graph object.")
  }
  g <- graph_from_adjacency_matrix(ceiling(object[[graph.name]]))
  part <- leiden_find_partition(g, resolution_parameter = resolution,
                                seed = random.seed, verbose = verbose,
                                partition_type = partition.type,
                                num_iter = n.iter)
  clusts <- part$membership
  names(clusts) <- Cells(object)
  clusts <- GroupSingletons2(clusts, object[[graph.name]],
                             min_cluster_size = min_cluster_size,
                             group.singletons = group.singletons,
                             verbose = verbose)
  if (group.singletons || !any(unique(clusts) == "singleton")) {
    clusts <- factor(as.character(clusts), 
                     levels = sort(unique(as.numeric(clusts))))
  } else {
    clusts <- factor(clusts, levels = c(seq_len(length(unique(clusts)) - 1), "singleton"))
  }
  object <- AddMetaData(object, clusts, 
                        col.name = paste0(graph.name, "_res.", resolution))
  object <- AddMetaData(object, clusts, 
                        col.name = "seurat_clusters")
  Idents(object) <- "seurat_clusters"
  object
}

#' Logistic regressionn DE test for spliced and unspliced
#' 
#' This function tries to find genes that effectively delineates between two
#' populations of cells on the spliced and unspliced phase portrait. This
#' function is a slight modification from LRDETest from Seurat that does LRDE
#' test one gene at a time, with only one matrix.
#' 
#' @param spliced Spliced matrix
#' @param unspliced Unspliced matrix
#' @param cells.1 Vector of cells in group 1
#' @param cells.2 Vector of cells in group 2
#' @param latent.vars Latent variables to include in model. Should be a data
#' frame or a matrix with rownames as cell barcodes.
#' @param verbose Print messages
LRDETest2 <- function(
  spliced, unspliced,
  cells.1,
  cells.2,
  latent.vars = NULL,
  verbose = TRUE
) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  spliced <- spliced[, rownames(group.info), drop = FALSE]
  unspliced <- unspliced[, rownames(group.info), drop = FALSE]
  latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
  my.sapply <- ifelse(
    test = verbose && nbrOfWorkers() == 1,
    yes = pbsapply,
    no = future_sapply
  )
  p_val <- my.sapply(
    X = 1:nrow(x = spliced),
    FUN = function(x) {
      if (is.null(x = latent.vars)) {
        model.data <- cbind(S = spliced[x, ], 
                            U = unspliced[x, ], 
                            group.info)
        fmla <- as.formula(object = "group ~ S + U")
        fmla2 <- as.formula(object = "group ~ 1")
      } else {
        model.data <- cbind(S = spliced[x, ], 
                            U = unspliced[x, ], 
                            group.info, latent.vars)
        fmla <- as.formula(object = paste(
          "group ~ S + U +",
          paste(colnames(x = latent.vars), collapse = "+")
        ))
        fmla2 <- as.formula(object = paste(
          "group ~",
          paste(colnames(x = latent.vars), collapse = "+")
        ))
      }
      model1 <- glm(formula = fmla, data = model.data, family = "binomial")
      model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
      lrtest <- lrtest(model1, model2)
      return(lrtest$Pr[2])
    }
  )
  to.return <- data.frame(p_val, row.names = rownames(spliced))
  return(to.return)
}

#' Just copied from Seurat's source code and stripped down

FindMarkers2 <- function(spliced, unspliced, ...) {
  UseMethod(generic = "FindMarkers2", object = spliced)
}

FindMarkers2.default <- function(
  spliced, unspliced,
  cells.1 = NULL,
  cells.2 = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1
) {
  features <- features %||% rownames(x = spliced)
  # error checking
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } 
  # feature selection (based on percentages)
  get_pct <- function(m, features, cells) {
    thresh.min <- 0
    round(
      x = rowSums(x = m[features, cells, drop = FALSE] > thresh.min) /
        length(x = cells),
      digits = 3
    )
  }
  thresh.min <- 0
  pct.1s <- get_pct(spliced, features, cells.1)
  pct.2s <- get_pct(spliced, features, cells.2)
  pct.1u <- get_pct(unspliced, features, cells.1)
  pct.2u <- get_pct(unspliced, features, cells.2)
  data.alpha1 <- cbind(pct.1s, pct.2s)
  data.alpha2 <- cbind(pct.1u, pct.2u)
  colnames(x = data.alpha1) <- c("pct.1s", "pct.2s")
  colnames(x = data.alpha2) <- c("pct.1u", "pct.2u")
  alpha.min1 <- apply(X = data.alpha1, MARGIN = 1, FUN = max)
  alpha.min2 <- apply(X = data.alpha2, MARGIN = 1, FUN = max)
  names(x = alpha.min1) <- names(alpha.min2) <- rownames(x = data.alpha1)
  features1 <- names(x = which(x = alpha.min1 > min.pct & alpha.min2 > min.pct))
  if (length(x = features) == 0) {
    stop("No features pass min.pct threshold")
  }
  alpha.diff1 <- alpha.min1 - apply(X = data.alpha1, MARGIN = 1, FUN = min)
  alpha.diff2 <- alpha.min2 - apply(X = data.alpha2, MARGIN = 1, FUN = min)
  features2 <- names(
    x = which(alpha.diff1 > min.diff.pct & alpha.diff2 > min.diff.pct)
  )
  features <- intersect(features1, features2)
  if (length(x = features) == 0) {
    stop("No features pass min.diff.pct threshold")
  }
  # feature selection (based on average difference)
  do_means <- function(data, features, cells, pseudocount.use) {
    apply(
      X = data[features, cells, drop = FALSE],
      MARGIN = 1,
      FUN = function(x) log(x = mean(x = expm1(x = x)) + pseudocount.use)
    )
  }
  data.1s <- do_means(spliced, features, cells.1, pseudocount.use)
  data.2s <- do_means(spliced, features, cells.2, pseudocount.use)
  data.1u <- do_means(unspliced, features, cells.1, pseudocount.use)
  data.2u <- do_means(unspliced, features, cells.2, pseudocount.use)
  total.diff1 <- (data.1s - data.2s)
  total.diff2 <- (data.1u - data.2u)
  features.diff <- if (only.pos) {
    names(x = which(x = total.diff1 > logfc.threshold & total.diff2 > logfc.threshold))
  } else {
    names(x = which(x = abs(x = total.diff1) > logfc.threshold & abs(total.diff2) > logfc.threshold))
  }
  features <- intersect(x = features, y = features.diff)
  if (length(x = features) == 0) {
    stop("No features pass logfc.threshold threshold")
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    # Should be cells.1 and cells.2?
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  # perform DE
  de.results <- LRDETest2(
    spliced = spliced[features, c(cells.1, cells.2), drop = FALSE],
    unspliced = unspliced[features, c(cells.1, cells.2), drop = FALSE],
    cells.1 = cells.1,
    cells.2 = cells.2,
    latent.vars = latent.vars,
    verbose = verbose
  )
  de.results[, "avg_logFC_s"] <- total.diff1[rownames(x = de.results)]
  de.results[, "avg_logFC_u"] <- total.diff2[rownames(x = de.results)]
  de.results <- cbind(de.results, data.alpha1[rownames(x = de.results), , drop = FALSE])
  de.results <- cbind(de.results, data.alpha2[rownames(x = de.results), , drop = FALSE])
  if (only.pos) {
    de.results <- de.results[de.results[, "avg_logFC_s"] > 0 | de.results[, "avg_logFC_u"] > 0, , drop = FALSE]
  }
  de.results <- de.results[order(de.results$p_val, 
                                 -de.results[, "avg_logFC_s"]
                                 -de.results[, "avg_logFC_u"]), ]
  de.results$p_val_adj = p.adjust(
    p = de.results$p_val,
    method = "bonferroni",
    n = nrow(x = spliced) # Seurat did it this way, 
    # though I think it's more conservative than it needs to be
  )
  return(de.results)
}

#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param reduction Reduction to use in differential expression testing - will test for DE on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping. Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from; note that if \code{test.use} is "negbinom", "poisson", or "DESeq2",
#' \code{slot} will be set to "counts"
#'
#' @importFrom methods is
#'
#' @rdname FindMarkers
#' @export
#' @method FindMarkers Seurat
#'
FindMarkers2.Seurat <- function(
  spliced, unspliced,
  ident.1 = NULL,
  ident.2 = NULL,
  subset.ident = NULL,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  ...
) {
  data.slot <- "data"
  assay <- assay %||% DefaultAssay(object = spliced)
  data.use1 <-  GetAssayData(object = spliced[[assay]], slot = data.slot)
  data.use2 <- GetAssayData(object = unspliced[[assay]], slot = data.slot)
  if (!identical(dimnames(data.use1), dimnames(data.use2))) {
    stop("spliced and unspliced matrices must have the same dimnames")
  }
  ident.1 <- WhichCells(object = spliced, idents = ident.1)
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(x = data.use1))) {
    bad.cells <- colnames(x = data.use1)[which(!as.character(x = ident.2) %in% colnames(x = data.use1))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = colnames(x = data.use1), y = ident.1)
    } else {
      ident.2 <- WhichCells(object = spliced, idents = ident.2)
    }
  }
  if (!is.null(x = latent.vars)) {
    latent.vars1 <- FetchData(
      object = spliced,
      vars = latent.vars,
      cells = c(ident.1, ident.2)
    )
    latent.vars2 <- FetchData(
      object = unspliced,
      vars = latent.vars,
      cells = c(ident.1, ident.2)
    )
    names(latent.vars2) <- paste0(names(latent.vars2), "_u")
    latent.vars <- cbind(latent.vars1, latent.vars2)
  }
  de.results <- FindMarkers2(
    data.use1, data.use2,
    cells.1 = ident.1,
    cells.2 = ident.2,
    features = features,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    min.diff.pct = min.diff.pct,
    verbose = verbose,
    only.pos = only.pos,
    max.cells.per.ident = max.cells.per.ident,
    random.seed = random.seed,
    latent.vars = latent.vars,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    pseudocount.use = pseudocount.use
  )
  return(de.results)
}

FindAllMarkers2 <- function(
  spliced, unspliced,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 1e-2,
  ...
) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  idents.all <- sort(x = unique(x = Idents(object = spliced)))
  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        FindMarkers2(
          spliced, unspliced,
          assay = assay,
          ident.1 = idents.all[i],
          ident.2 = NULL,
          features = features,
          logfc.threshold = logfc.threshold,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          verbose = verbose,
          only.pos = only.pos,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = random.seed,
          latent.vars = latent.vars,
          min.cells.feature = min.cells.feature,
          min.cells.group = min.cells.group,
          pseudocount.use = pseudocount.use,
          ...
        )
      },
      error = function(cond) {
        return(cond$message)
      }
    )
    if (class(x = genes.de[[i]]) == "character") {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  return(gde.all)
}

# Find DE genes between subsets of clusters

FindMarkers_clusts <- function(spliced, unspliced, 
                               clusts.1, clusts.2,
                               latent.vars = NULL, gns = NULL, ...) {
  spliced$neu_cat1 <- case_when(
    Idents(spliced) %in% clusts.2 ~ "two",
    Idents(spliced) %in% clusts.1 ~ "one",
    TRUE ~ "other"
  )
  Idents(spliced) <- "neu_cat1"
  out <- FindMarkers2(spliced, unspliced, ident.1 = "one", ident.2 = "two", 
                      latent.vars = latent.vars, ...)
  if (!is.null(gns)) {
    out <- out %>% 
      mutate(., gene = rownames(.)) %>% 
      left_join(gns, by = "gene")
  }
  return(out)
}
