## This file is provided as an addition to the Jupyter Notebook code supplementing our manuscript and shared on GitHub at https://github.com/InsermRightLab/LILRB2-Varin-et-al
## It is used to load all necessary packages as well as build all functions necessary for the Jupyter Notebook

# Loading packages

invisible(suppressPackageStartupMessages(sapply(c("patchwork","venn","ggborderline","Seurat","slingshot","tradeSeq","batchtools",
                                                  "BiocParallel","RightOmicsTools","ggplot2"),library,character.only = T)))








# Custom function to visualize colored lineages on UMAP coordinates and show pseudotime

Right_LineagePlot <- function(sds,
                              models = NULL,
                              lineages = NULL,
                              clusters = NULL,
                              knots = NULL,
                              pt.size = 1,
                              clusters.colors = "lightgrey",
                              lineages.colors = "black") {
  
  if (isTRUE("SingleCellExperiment" %in% class(models))) {
    if (is.numeric(knots)) {
      knots = unname(metadata(models)$tradeSeq$knots)[knots]
    }
    else {
      knots = unname(metadata(models)$tradeSeq$knots)
    }
  }
  
  df <- data.frame(slingReducedDim(sds))
  df <- cbind(df, data.frame(slingPseudotime(sds)))
  df$clusters <- "none"
  if (is.character(clusters) | is.factor(clusters)) {
    if (clusters == "pseudotime") {
      df$clusters = apply(df[ , grep("Lineage", colnames(df))], 1, max, na.rm = TRUE)
      ggcolors <- scale_color_viridis_c(limits = c(0, ceiling(max(df$clusters, na.rm = TRUE))))
    }
    else if (length(clusters) == 1) {
      df$clusters <- colData(sds)[ , clusters]
      ggcolors <- scale_color_manual(values = clusters.colors)
    }
    else if (length(clusters) == ncol(sds)) {
      df$clusters <- clusters
      ggcolors <- scale_color_manual(values = clusters.colors)
    }
    else error()
  }
  else {
    ggcolors <- scale_color_manual(values = clusters.colors) 
  }
  
  curve <- slingCurves(sds)
  if (is.numeric(lineages)) {
    curve <- curve[lineages]
  }
  else {
    lineages <- seq_along(curve)
  }
  
  gg <- ggplot(df, aes(x = .data[[colnames(df)[1]]], y = .data[[colnames(df)[2]]], col = clusters)) +
    geom_point(size = pt.size, show.legend = ifelse(is.null(clusters), FALSE, TRUE)) +
    ggcolors +
    labs(col = ifelse(clusters == "pseudotime", "Pseudotime", "Clusters")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.position = "bottom",
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  if (isTRUE("SingleCellExperiment" %in% class(models))) {
    knots = do.call(rbind, lapply(seq_along(curve), function(pts) {
      cu <- as.data.frame(curve[[pts]]$s[curve[[pts]]$ord, ])
      time <- df[ , grep(paste0("Lineage", lineages[pts]), colnames(df))]
      cu$time <- seq(min(time, na.rm = T), max(time, na.rm = T), length.out = nrow(cu))
      pts <- sapply(knots, function(pt) {
        as.numeric(ifelse(as.numeric(pt) > max(cu$time), nrow(cu),
                          ifelse(as.numeric(pt) < min(cu$time), 1,
                                 Position(function(pos) pos >= as.numeric(pt), cu$time))))
      })
      return(cu[pts, 1:2, drop = F])
    }))
    for (cu in seq_along(curve)) {
      gg <- gg + geom_path(data = as.data.frame(curve[[cu]]$s[curve[[cu]]$ord, ]),
                           col = ifelse(length(lineages.colors) == 1, lineages.colors, lineages.colors[cu]), linewidth = 2)
    }
    gg <- gg + geom_point(data = knots, col = "black", size = 4)
  }
  else {
    for (cu in seq_along(curve)) {
      gg <- gg + geom_borderpath(data = as.data.frame(curve[[cu]]$s[curve[[cu]]$ord, ]),
                                 col = ifelse(length(lineages.colors) == 1, lineages.colors, lineages.colors[cu]),
                                 linewidth = 2, arrow = arrow(angle = 30, length = unit(0.25, "inches"),
                                                              ends = "last", type = "open"))
    }
  }
  
  if (!is.null(clusters)) {
    if (clusters != "pseudotime") {
      gg <- gg + guides(colour = guide_legend(nrow = ifelse(length(unique(df$clusters)) > 3, ceiling(length(unique(df$clusters))/3), 1),
                                              title.position = "top", title.hjust = 0.5,
                                              override.aes = list(size = 4)))
    }
    else {
    gg <- gg + guides(colour = guide_colourbar(direction="vertical"))
    }
  }
  gg
}

                   
                   
                   
                   
                   
                   
                   
# Custom function to visualize Wald statistics, log fold change and rank score of each gene as scatterplot                   
                   
Right_WaldPlot = function(test.results, score.color = "red") {
  
  gg <- lapply(seq_along(test.results), function(lin) {
    df <- test.results[[lin]]
    gg <- ggplot(data=df, aes(x=.data[[grep("ogFC", colnames(df), value = T)]], y=waldStat, col=transientScore)) +
      labs(x='Median Fold Change',
           y='EarlyDETest Wald Statistic', col='Rank Score')  +
      geom_point(size=1, shape=16) +
      scale_color_gradientn(colours = c("grey40", "lightgrey", score.color),
                            labels = c("q10", "median", "q90"),
                            breaks = c(quantile(df$transientScore, probs = 0.1),
                                       median(df$transientScore),
                                       quantile(df$transientScore, probs = 0.90))) +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5, size = 14)) +
      annotate("text", x = 50, y = 500, size = 8, color = score.color,
               label = paste0("N=", nrow(df[df$transientScore >= quantile(df$transientScore, probs = 0.90), ]))) +
      ggtitle(names(test.results)[lin]) +
      scale_x_log10(expand=c(0,0), limits = c(0.01, 1000), breaks = c(0.01,0.1,1,10,100,1000),
                    labels=function(x) format(x, scientific = FALSE)) +
      scale_y_log10(limit=c(floor(min(df$waldStat[df$waldStat != min(df$waldStat)])),10000),
                    expand=c(0,0), labels=function(x) format(x, scientific = FALSE))
  })
  return(gg)
}

                    
                    
                    
                    
                    
                    
                    
# Full code for tradeSeq's function fitGAM(), with modification on line 473 to allow batchtools to parallelize model fitting without errors 
                    
.assignCells <- function(cellWeights) {
  if (is.null(dim(cellWeights))) {
    if (any(cellWeights == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      return(matrix(1, nrow = length(cellWeights), ncol = 1))
    }
  } else {
    if (any(rowSums(cellWeights) == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      # normalize weights
      normWeights <- sweep(cellWeights, 1,
                           FUN = "/",
                           STATS = apply(cellWeights, 1, sum)
      )
      # sample weights
      wSamp <- apply(normWeights, 1, function(prob) {
        stats::rmultinom(n = 1, prob = prob, size = 1)
      })
      # If there is only one lineage, wSamp is a vector so we need to adjust for that
      if (is.null(dim(wSamp))) {
        wSamp <- matrix(wSamp, ncol = 1)
      } else {
        wSamp <- t(wSamp)
      }
      return(wSamp)
    }
  }
}


.checks <- function(pseudotime, cellWeights, U, counts, conditions, family) {

  # counts must only have positive integer values
  if(family == "nb"){
    if (any(counts < 0)) {
      stop("All values of the count matrix should be non-negative")
    }
  }
  

  # check if pseudotime and weights have same dimensions.
  if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
    if (!identical(dim(pseudotime), dim(cellWeights))) {
      stop("pseudotime and cellWeights must have identical dimensions.")
    }
  }

  # check if dimensions of U and counts agree
  if (!is.null(U)) {
    if (!(nrow(U) == ncol(counts))) {
      stop("The dimensions of U do not match those of counts.")
    }
  }

  # check if dimensions for counts and pseudotime / cellweights agree
  if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
    if (!identical(nrow(pseudotime), ncol(counts))) {
      stop("pseudotime and count matrix must have equal number of cells.")
    }
    if (!identical(nrow(cellWeights), ncol(counts))) {
      stop("cellWeights and count matrix must have equal number of cells.")
    }
  }

  if(!is.null(conditions)){
    if(!is(conditions, "factor")){
      stop("conditions must be a vector of class factor.")
    }
  }
  
  if(any(is.na(pseudotime)[cellWeights > 0])){
    stop("Pseudotime contains NA values for non-zero weights.")
  }
  
  if(any(is.na(pseudotime))){
    warning("Pseudotime contains NA values.")
  }
}

.get_offset <- function(offset, counts) {
  if (is.null(offset)) {
    nf <- try(edgeR::calcNormFactors(counts), silent = TRUE)
    if (is(nf, "try-error")) {
      message("TMM normalization failed. Will use unnormalized library sizes",
              "as offset.\n")
      nf <- rep(1,ncol(counts))
    }
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
    if(any(libSize == 0)){
      message("Some library sizes are zero. Offsetting these to 1.\n")
      offset[libSize == 0] <- 0
    }
  }
  return(offset)
}

# TODO: make sure error messages in fitting are silent,
# but print summary at end.
# TODO: make sure warning message for knots prints after looping

.findKnots <- function(nknots, pseudotime, wSamp) {
  # Easier to recreate them all here than to pass them on
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }

  # Get the times for the knots
  tAll <- c()
  for (ii in seq_len(nrow(pseudotime))) {
    tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii,]))]
  }

  knotLocs <- stats::quantile(tAll, probs = (0:(nknots - 1)) / (nknots - 1))
  if (any(duplicated(knotLocs))) {
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest lineage
    knotLocs <- stats::quantile(t1[l1 == 1],
                                probs = (0:(nknots - 1)) / (nknots - 1))
    # if duplication still occurs, get average btw 2 points for dups.
    if (any(duplicated(knotLocs))) {
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if (max(which(dupId)) == length(knotLocs)) {
        dupId <- duplicated(knotLocs, fromLast = TRUE)
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      } else {
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      }
    }
    # if this doesn't fix it, get evenly spaced knots with warning
    if (any(duplicated(knotLocs))) {
      knotLocs <- seq(min(tAll), max(tAll), length = nknots)
    }
  }

  maxT <- max(pseudotime[,1])
  if (ncol(pseudotime) > 1) {
    maxT <- c()
    # note that first lineage should correspond to the longest, hence the
    # 100% quantile end point is captured.
    for (jj in 2:ncol(pseudotime)) {
      maxT[jj - 1] <- max(get(paste0("t", jj))[get(paste0("l",jj)) == 1])
    }
  }
  # if max is already a knot we can remove that
  if (all(maxT %in% knotLocs)) {
    knots <- knotLocs
  } else {
    maxT <- maxT[!maxT %in% knotLocs]
    replaceId <- vapply(maxT, function(ll){
      which.min(abs(ll - knotLocs))
    }, FUN.VALUE = 1)
    knotLocs[replaceId] <- maxT
    if (!all(maxT %in% knotLocs)) {
      # if not all end points are at knots, return a warning, but keep
      # quantile spaced knots.
      warning(paste0("Impossible to place a knot at all endpoints.",
                     "Increase the number of knots to avoid this issue."))
    }
    knots <- knotLocs
  }

  # guarantees that first knot is 0 and last knot is maximum pseudotime.
  knots[1] <- min(tAll)
  knots[nknots] <- max(tAll)

  knotList <- lapply(seq_len(ncol(pseudotime)), function(i){
    knots
  })
  names(knotList) <- paste0("t", seq_len(ncol(pseudotime)))

  return(knotList)
}

.fitGAM <- function(counts, U = NULL, pseudotime, cellWeights,
                    conditions,
                    genes = seq_len(nrow(counts)),
                    weights = NULL, offset = NULL, nknots = 6, verbose = TRUE,
                    parallel = FALSE, BPPARAM = BiocParallel::bpparam(),
                    aic = FALSE, control = mgcv::gam.control(), sce = TRUE,
                    family = "nb", gcv = FALSE){
  
  if(!is.null(conditions)){
    message("Fitting lineages with multiple conditions. This method has ",
            "been tested on a couple of datasets, but is still in an ",
            "experimental phase.")
  }

  if (is(genes, "character")) {
    if (!all(genes %in% rownames(counts))) {
      stop("The genes ID is not present in the models object.")
    }
    if(any(duplicated(genes))){
      stop("The genes vector contains duplicates.")
    }
    id <- match(genes, rownames(counts))
  } else {
    id <- genes
  }

  if (parallel) {
    BiocParallel::register(BPPARAM)
    if (verbose) {
      # update progress bar 40 times
      BPPARAM$tasks = as.integer(40)
      # show progress bar
      BPPARAM$progressbar = TRUE
    }
  }


  # Convert pseudotime and weights to matrices if need be
  if (is.null(dim(pseudotime))) {
    pseudotime <- matrix(pseudotime, nrow = length(pseudotime))
  }
  if (is.null(dim(cellWeights))) {
    cellWeights <- matrix(cellWeights, nrow = length(cellWeights))
  }

  .checks(pseudotime, cellWeights, U, counts, conditions, family)

  wSamp <- .assignCells(cellWeights)

  # define pseudotime for each lineage
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  # get lineage indicators for cells to use in smoothers
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }

  # offset
  offset <- .get_offset(offset, counts)

  # fit model
  ## fixed effect design matrix
  if (is.null(U)) {
    U <- matrix(rep(1, nrow(pseudotime)), ncol = 1)
  }

  ## Get the knots
  knotList <- .findKnots(nknots, pseudotime, wSamp)

  ## fit NB GAM
  ### Actually fit the model ----
  teller <- 0
  converged <- rep(TRUE, length(genes))
  counts_to_Gam <- function(y) {
    teller <<- teller + 1
    # define formula (only works if defined within apply loop.)
    nknots <- nknots
    if (!is.null(weights)) weights <- weights[teller,]
    if (!is.null(dim(offset))) offset <- offset[teller,]
    if(is.null(conditions)){
      smoothForm <- stats::as.formula(
        paste0("y ~ -1 + U + ",
               paste(vapply(seq_len(ncol(pseudotime)), function(ii){
                 paste0("s(t", ii, ", by=l", ii, ", bs='cr', id=1, k=nknots)")
               }, FUN.VALUE = "formula"),
               collapse = "+"), " + offset(offset)")
      )
    } else {
      for(jj in seq_len(ncol(pseudotime))){
        for(kk in seq_len(nlevels(conditions))){
          # three levels doesn't work. split it up and loop over both conditions and pseudotime
          # to get a condition-and-lineage-specific smoother. Also in formula.
          lCurrent <- get(paste0("l", jj))
          id1 <- which(lCurrent == 1)
          lCurrent[id1] <- ifelse(conditions[id1] == levels(conditions)[kk], 1, 0)
          assign(paste0("l", jj, "_", kk), lCurrent)
        }
      }
      smoothForm <- stats::as.formula(
        paste0("y ~ -1 + U + ",
               paste(vapply(seq_len(ncol(pseudotime)), function(ii){
                 paste(vapply(seq_len(nlevels(conditions)), function(kk){
                   paste0("s(t", ii, ", by=l", ii, "_", kk,
                          ", bs='cr', id=1, k=nknots)")
                 }, FUN.VALUE = "formula"),
                 collapse = "+")
               }, FUN.VALUE = "formula"),
               collapse="+")
               , " + offset(offset)")
      )
    }
    # fit smoother, catch errors and warnings
    s <- mgcv::s
    m <- suppressWarnings(try(withCallingHandlers({
      mgcv::gam(smoothForm, family = family, knots = knotList, weights = weights,
                control = control)},
      error = function(e){ #if errors: return try-error class
        converged[teller] <<- FALSE
        return(structure("Fitting errored",
                         class = c("try-error", "character")))
      },
      warning = function(w){ #if warning: set converged to FALSE
        converged[teller] <<- FALSE
      }), silent=TRUE))
    return(m)
  }

  #### fit models
  if (parallel) {
    gamList <- BiocParallel::bplapply(
      as.list(as.data.frame(t(as.matrix(counts)[id, ]))), # We added as.list() to fix a batchtools error, no other change in the function
      counts_to_Gam, BPPARAM = BPPARAM
    )
  } else {
    if (verbose) {
      gamList <- pbapply::pblapply(
        as.data.frame(t(as.matrix(counts)[id, ])),
        counts_to_Gam
      )
    } else {
      gamList <- lapply(
        as.data.frame(t(as.matrix(counts)[id, ])),
        counts_to_Gam
      )
    }
  }

  ### output
  if (aic) { # only return AIC
    # return(unlist(lapply(gamList, function(x){
    #   if (class(x)[1] == "try-error") return(NA)
    #   x$aic
    # })))
    aicVals <- unlist(lapply(gamList, function(x){
      if (class(x)[1] == "try-error") return(NA)
      x$aic
    }))
    if (gcv) {
      gcvVals <- unlist(lapply(gamList, function(x){
        if (class(x)[1] == "try-error") return(NA)
        x$gcv.ubre
      }))
      return(list(aicVals, gcvVals))
    } else return(aicVals)
  }

  if (sce) { #tidy output: also return X
    # tidy smoother regression coefficients
    betaAll <- lapply(gamList, function(m) {
      if (is(m, "try-error")) {
        beta <- NA
      } else {
        beta <- matrix(stats::coef(m), ncol = 1)
        rownames(beta) <- names(stats::coef(m))
      }
      return(beta)
    })
    betaAllDf <- data.frame(t(do.call(cbind,betaAll)))
    rownames(betaAllDf) <- rownames(counts)[id]

    # list of variance covariance matrices
    SigmaAll <- lapply(gamList, function(m) {
      if (is(m, "try-error")) {
        Sigma <- NA
      } else {
        Sigma <- m$Vp
      }
      return(Sigma)
    })

    # Get X, dm and knotPoints
    element <- min(which(!is.na(SigmaAll)))
    m <- gamList[[element]]
    X <- stats::predict(m, type = "lpmatrix")
    dm <- m$model[, -1]
    knotPoints <- m$smooth[[1]]$xp

    # return output
    return(list(beta = betaAllDf,
                Sigma = SigmaAll,
                X = X,
                dm = dm,
                knotPoints = knotPoints,
                converged = converged)
           )
  } else {
    return(gamList)
  }
}

fitGAM = function(counts,
                  sds = NULL,
                  pseudotime = NULL,
                  cellWeights = NULL,
                  conditions = NULL,
                  U = NULL,
                  genes = seq_len(nrow(counts)),
                  weights = NULL,
                  offset = NULL,
                  nknots = 6,
                  verbose = TRUE,
                  parallel = FALSE,
                  BPPARAM = BiocParallel::bpparam(),
                  control = mgcv::gam.control(),
                  sce = TRUE,
                  family = "nb",
                  gcv = FALSE){
    
            if (is.null(counts)) stop("Provide expression counts using counts",
                                      " argument.")

            ## either pseudotime or slingshot object should be provided
            if (is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))) {
              stop("Either provide the slingshot object using the sds ",
                   "argument, or provide pseudotime and cell-level weights ",
                   "manually using pseudotime and cellWeights arguments.")
            }
            if (!is.null(sds)) {
              # check if input is slingshotdataset or pseudotimeordering
              if (is(sds, "SlingshotDataSet") | is(sds, "PseudotimeOrdering")) {
                if (!sce) {
                  warning(paste0(
                    "If an sds argument is provided, the sce argument is ",
                    "forced to TRUE "))
                  sce <- TRUE
                }
              } else stop("sds argument must be a SlingshotDataSet or ",
                          "PseudotimeOrdering object.")
              # extract variables from slingshotdataset
              pseudotime <- slingPseudotime(sds, na = FALSE)
              cellWeights <- slingCurveWeights(sds)
            }

            if(any(is.na(pseudotime))){
              stop("The pseudotimes contain NA values, and these cannot be used",
              " for GAM fitting.")
            }
            if(any(is.na(cellWeights))){
              stop("The cellWeights contain NA values, and these cannot be used",
                   " for GAM fitting.")
            }

            if(!is.null(conditions)){
              if(class(conditions) != "factor") stop("conditions must be a factor vector.")
              if(length(conditions) != ncol(counts)){
                stop("conditions vector must have same length as number of cells.")
              }
              if(nlevels(conditions) == 1) {
                message("Only one condition was provided. Will run fitGAM without conditions")
                conditions <- NULL
              }
              if (!sce) {
                warning(paste0("If conditions, tradeSeq will return",
                               " a SingleCellExperiment object"))
              }
            }

            gamOutput <- .fitGAM(counts = counts,
                                 U = U,
                                 pseudotime = pseudotime,
                                 cellWeights = cellWeights,
                                 conditions = conditions,
                                 genes = genes,
                                 weights = weights,
                                 offset = offset,
                                 nknots = nknots,
                                 verbose = verbose,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 control = control,
                                 sce = sce,
                                 family = family,
                                 gcv = gcv)

            # old behaviour: return list
            if (!sce) {
              return(gamOutput)
            }

            # return SingleCellExperiment object
            sc <- SingleCellExperiment(assays = list(counts = counts[genes,]))
            # slingshot info
            SummarizedExperiment::colData(sc)$crv <- S4Vectors::DataFrame(
              pseudotime = pseudotime,
              cellWeights = cellWeights)
            # tradeSeq gene-level info
            df <- tibble::enframe(gamOutput$Sigma, value = "Sigma")
            df$beta <- tibble::tibble(beta = gamOutput$beta)
            df$converged <- gamOutput$converged
            suppressWarnings(rownames(df) <- rownames(counts[genes,]))
            SummarizedExperiment::rowData(sc)$tradeSeq <- df
            # tradeSeq cell-level info
            if(is.null(conditions)){
              SummarizedExperiment::colData(sc)$tradeSeq <-
                tibble::tibble(X = gamOutput$X, dm = gamOutput$dm)
            } else {
              SummarizedExperiment::colData(sc)$tradeSeq <-
                tibble::tibble(X = gamOutput$X, dm = gamOutput$dm,
                               conditions = conditions)
            }
            # metadata: tradeSeq knots
            S4Vectors::metadata(sc)$tradeSeq <-
              list(knots = gamOutput$knotPoints)
            return(sc)
          }

                    
                    
                    
                    
                    
                
                    
# Full code for tradeSeq's function evaluateK(), no modification, uses .fitGAM() internally and needs to use the modified version for batchtools parallelization
                    
.evaluateK <- function(counts, U = NULL, pseudotime, cellWeights, plot = TRUE,
                       nGenes = 500, k = 3:10, weights = NULL, offset = NULL,
                       aicDiff = 2, verbose = TRUE, conditions, parallel,
                       BPPARAM, family = "nb", gcv = FALSE, ...) {

  if (any(k < 3)) stop("Cannot fit with fewer than 3 knots, please increase k.")
  if (length(k) == 1) stop("There should be more than one k value")
  ## calculate offset on full matrix
  if (is.null(offset)) {
    nf <- try(edgeR::calcNormFactors(counts))
    if (is(nf, "try-error")) {
      message("TMM normalization failed. Will use unnormalized library sizes",
              "as offset.")
      nf <- rep(1,ncol(counts))
    }
    libSize <- colSums(as.matrix(counts)) * nf
    offset <- log(libSize)
  }

  ## AIC over knots
  geneSub <- sample(seq_len(nrow(counts)), nGenes)
  countSub <- counts[geneSub,]
  weightSub <- weights[geneSub,]
  kList <- list()
  for (ii in seq_len(length(k))) kList[[ii]] <- k[ii]
  #gamLists <- BiocParallel::bplapply(kList, function(currK){
  aicVals <- lapply(kList, function(currK){
    gamAIC <- .fitGAM(counts = countSub, U = U, pseudotime = pseudotime,
                      cellWeights = cellWeights, conditions = conditions,
                      nknots = currK, verbose = verbose, sce = FALSE, 
                      parallel = parallel, BPPARAM = BPPARAM,
                      weights = weightSub, offset = offset, aic = TRUE, 
                      family = family, gcv = gcv, ...)
  })
  #, BPPARAM = MulticoreParam(ncores))

  if (gcv) {
    aicMat <- do.call(cbind, lapply(aicVals, "[[", 1))
    colnames(aicMat) <- paste("k:", k)
    gcvMat <- do.call(cbind, lapply(aicVals, "[[", 2))
    colnames(gcvMat) <- paste("k:", k)
  } else {
    aicMat <- do.call(cbind,aicVals)
    colnames(aicMat) <- paste("k:", k)
  }


  if (plot) {
    plot_evalutateK_results(aicMat = aicMat, k = k, aicDiff = aicDiff) 
  }

  if(gcv){
    return(list(aic = aicMat,
                gcv = gcvMat))
  } else {
    return(aicMat)
  }
}

evaluateK = function(counts,
                     k = 3:10,
                     nGenes = 500,
                     sds = NULL,
                     pseudotime = NULL,
                     cellWeights = NULL,
                     U = NULL,
                     conditions = NULL,
                     plot = TRUE,
                     weights = NULL,
                     offset = NULL,
                     aicDiff = 2,
                     verbose = TRUE,
                     parallel = FALSE,
                     BPPARAM = BiocParallel::bpparam(),
                     control = mgcv::gam.control(),
                     family = "nb",
                     gcv = FALSE,
                     ...){

            ## either pseudotime or slingshot object should be provided
            if (is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))) {
              stop("Either provide the slingshot object using the sds ",
                   "argument, or provide pseudotime and cell-level weights ",
                   "manually using pseudotime and cellWeights arguments.")
            }

            if (!is.null(sds)) {
              # check if input is slingshotdataset or pseudotimeordering
              if (is(sds, "SlingshotDataSet") | is(sds, "PseudotimeOrdering")) {
                # extract variables from slingshotdataset
                pseudotime <- slingPseudotime(sds, na = FALSE)
                cellWeights <- slingCurveWeights(sds)
              }
              else stop("sds argument must be a SlingshotDataSet or ",
                        "PseudotimeOrdering object.")
            }

            if (is.null(counts)) stop("Provide expression counts using counts",
                                      " argument.")

            aicOut <- .evaluateK(counts = counts,
                                 k = k,
                                 U = U,
                                 pseudotime = pseudotime,
                                 cellWeights = cellWeights,
                                 plot = plot,
                                 nGenes = nGenes,
                                 weights = weights,
                                 offset = offset,
                                 verbose = verbose,
                                 conditions = conditions,
                                 parallel = parallel,
                                 BPPARAM = BPPARAM,
                                 control = control,
                                 family = family,
                                 gcv = gcv,
                                 ...)

            return(aicOut)

          }
