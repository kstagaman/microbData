#' @title Transform Feature Counts
#' @description Functions for transforming feature counts in a \code{microbData} object.
#' @aliases transform_counts
#' @aliases rarefy
#' @aliases center.log.ratio
#' @aliases variance.stabilize
#' @aliases relative.abundance
#' @param mD required; the \code{microbData} object with feature counts to be transformed.
#' @param update.mD logical; should this function return a new \code{microbData} object with the transformed abundances replacing the original abundances (TRUE) or just the resulting transformed abundance table (FALSE)? Default is TRUE.
#' @param quiet logical; should informational output (not warnings or errors) be suppressed? Default is FALSE.
#' @param f required for \code{transform_counts}; a function that can be applied to the samples (rows) of the Abundance table.
#' @param iters integer; for \code{rarefy}, the number of times to subsample the abundance table. Default is 999.
#' @param replace.with character; for \code{rarefy}, whether to return the first subsampling as the new abundances table in the microbData object, or to return a table of average abundances (may cause issues with other functions). Default is "first"
#' @param keep.tables logical; for \code{rarefy}, keep all the subsampled abundances tables?
#' @param min.abund integer; for \code{rarefy}, if not NULL, samples will be rarefied to the same depth as the lowest sample abundance equal to or greater than this integer. E.g., if \code{min.abund = 10000} and there four samples with 8000, 11000, 12000, and 13000 total reads, respectively, the first sample will be dropped, and last three samples will be rarefied to 11000 reads. For \code{center.log.ratio}, the interger to pass to the \code{min.reads} argument in \code{\link[CoDaSeq]{codaSeq.filter}}. Default is 1e4.
#' @param exactly.to integer; for \code{rarefy}, if not NULL, samples will be rarefied to exactly this integer. This argument supersedes \code{min.abund}. Samples with total reads lower than this number will be dropped. Default is NULL.
#' @param alpha.metrics character; for \code{rarefy}, which alpha-diversity metrics, if any should be applied to the iterations of abundance subsamplings and the results averaged (mean). Default is NULL.
#' @param beta.metrics character; for \code{rarefy}, which beta-diversity metrics, if any should be applied to the iterations of abundance subsamplings and the results averaged (mean). Default is NULL.
#' @param trim.features logical; for \code{rarefy}, should features that are no longer present in any samples after rarefaction be dropped from the \code{microbData} object? Default is TRUE.
#' @param user.seed integer; for \code{rarefy}, a user-supplied random seed to make the random subsampling process repeatable. If NULL, will just use the default \code{\link{.Random.seed}}, which will be reported at the end. Default is NULL.
#' @param min.prop numeric; for \code{center.log.ratio}, the minimum proportional abundance of a read in any sample. (See \code{\link[CoDaSeq]{codaSeq.filter}}). Default is 0.001.
#' @param min.occur numeric; for \code{center.log.ratio}, the minimum fraction of non-0 reads for each variable in all samples. (See \code{\link[CoDaSeq]{codaSeq.filter}}). Default is 0.
#' @param smpls.by.row logical; for \code{center.log.ratio}, TRUE if rows contain samples, FALSE if rows contain variables. (See \code{\link[CoDaSeq]{codaSeq.filter}}). Default is TRUE.
#' @param method character; for \code{center.log.ratio}. Choose one of: count zero multiplicative ("CZM"); geometric Bayesian multiplicative ("GBM"); square root BM ("SQ"); Bayes-Laplace BM ("BL"); user-specified hyper-parameters ("user"). (See \code{\link[zCompositions]{cmultRepl}}). Default is "CZM"
#' @param lab numeric or character; for \code{center.log.ratio}, a unique label used to denote count zeros in X. (See \code{\link[zCompositions]{cmultRepl}}). Default is 0.
#' @param design formula or matrix; required for \code{variance.stabilize}, from \code{\link[DESeq2]{DESeqDataSetFromMatrix}}: the formula expresses how the counts for each gene depend on the variables in colData. Many R formula are valid, including designs with multiple variables, e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment. See results for a variety of designs and how to extract results tables. By default, the functions in this package will use the last variable in the formula for building results tables and plotting. ~ 1 can be used for no design, although users need to remember to switch to another design for differential testing.
#' @seealso \code{\link[CoDaSeq]{codaSeq.filter}}, \code{\link[zCompositions]{cmultRepl}}, \code{\link[DESeq2]{getVarianceStabilizedData}}
#' @examples
#' data("mD_raw")
#' # serial execution
#' mD.rar <- rarefy(mD.raw, keep.tables = T, alpha.metrics = "InvSimpson", beta.metrics = "Canberra")
#'
#' # parallel execution
#' cl <- parallel::makeCluster(4L)
#' doParallel::registerDoParallel(cl)
#' mD.rar <- rarefy(mD.raw, keep.tables = T, alpha.metrics = "InvSimpson", beta.metrics = "Canberra")
#' parallel::stopCluster()
#' doParallel::registerDoSEQ()
#'
#' @export

####################################
#' @name transform_counts
#' @title Transform Feature Counts
#' @description Supply a custom function with which to transform feature counts in a \code{microbData} object.
#' @rdname transform_counts
#' @export

transform_counts <- function(mD, f, update.mD = TRUE) {
  mD@Abundances <- apply(X = mD@Abundances, MARGIN = 1, FUN = f) %>%
    t() %>%
    .[, mD@Feature.names]
  if (update.mD) {
    mD <- add.other.data(x = f, name = "Abundances transformed", mD = mD)
    return(mD)
  } else {
    return(mD@Abundances)
  }
}

####################################
#' @name rarefy
#' @title Rarefy Feature Counts
#' @description Rarefy feature counts in a \code{microbData} object to the same depth over multiple iterations and, if wanted calculate average alpha- and beta-diversity statistics.
#' @rdname transform_counts
#' @export

rarefy <- function(
    mD,
    iters = 999,
    replace.with = c("first", "average"),
    keep.tables = FALSE,
    min.abund = 10000,
    exactly.to = NULL,
    alpha.metrics = NULL,
    beta.metrics = NULL,
    trim.features = TRUE,
    user.seed = NULL,
    quiet = FALSE,
    debug = FALSE
) {
  if (is.null(user.seed)) {
    if (!quiet) { rlang::inform(paste("Random seed:", .Random.seed[1])) }
  } else {
    set.seed(user.seed)
    if (!quiet) { rlang::inform(paste("Random seed:", user.seed)) }
  }
  replace.with <- rlang::arg_match(replace.with, values = c("first", "average"))
  if (!is.null(alpha.metrics)) { alpha.metrics %<>% purrr::set_names() }
  if (!is.null(beta.metrics)) { beta.metrics %<>% purrr::set_names() }
  rarefy.to <- ifelse(
    is.null(exactly.to),
    min(sample.sums(mD)[sample.sums(mD) >= min.abund]),
    exactly.to
  )
  if (!quiet) { rlang::inform(paste("Rarefying to:", rarefy.to)) }

  mD1 <- drop.samples(
    mD,
    samples = names(sample.sums(mD)[sample.sums(mD) < rarefy.to])
  )
  zero.mat <- copy(mD1@Abundances)
  zero.mat[,] <- 0

  if (iters > nsamples(mD1)) {
    mat.list <- foreach::foreach(
      i = 1:iters,
      .verbose = debug,
      .export = c("mD1", "rarefy.to")
    ) %dopar% {
      sub.list <- apply(X = mD1@Abundances, MARGIN = 1, simplify = F, FUN = function(x) {
        sample(names(x), size = rarefy.to, replace = T, prob = x) %>%
          table() %>%
          return()
      })
      mat.i <- zero.mat
      for (smpl in names(sub.list)) {
        mat.i[smpl, names(sub.list[[smpl]])] <- sub.list[[smpl]]
      }
      return(mat.i)
    }
  } else {
    splits <- split(mD1@Abundances, row(mD1@Abundances)) %>%
      setNames(rownames(mD1@Abundances)) %>%
      lapply(setNames, colnames(mD1@Abundances))
    mat.list <- lapply(1:iters, function(i) {
      sub.list <- foreach::foreach(
        x = splits,
        .final = function(x) setNames(x, names(splits)),
        .verbose = debug,
        .export = c("rarefy.to")
      ) %dopar% {
        sample(names(x), size = rarefy.to, replace = T, prob = x) %>%
          table() %>%
          return()
      }
      mat.i <- zero.mat
      for (smpl in names(sub.list)) {
        mat.i[smpl, names(sub.list[[smpl]])] <- sub.list[[smpl]]
      }
      return(mat.i)
    })
  }

  if (!is.null(alpha.metrics)) {
    if (!quiet) {
      rlang::inform(paste("Calculating alpha-diversities:", paste(alpha.metrics, collapse = ", ")))
    }
    alpha.res <- foreach::foreach(
      mat.i = mat.list,
      .final = rbindlist,
      .inorder = FALSE,
      .verbose = debug,
      .export = c("mD1", "alpha.metrics")
    ) %dopar% {
      mD.i <- replace.abundances(mD1, mat.i)
      metrics.i <- microbData::alpha.diversity(mD.i, metrics = alpha.metrics, update.mD = F)
    }
    new.meta <- alpha.res[, lapply(.SD, mean), by = c(mD1@Sample.col), .SDcols = alpha.metrics] %>%
      merge(mD1@Metadata, by = mD1@Sample.col)
    mD1 %<>% replace.metadata(new.tbl = new.meta)
  }
  if (!is.null(beta.metrics)) {
    if (!quiet) {
      rlang::inform(paste("Calculating beta-diversities:", paste(beta.metrics, collapse = ", ")))
    }
    dist.mats <- lapply(beta.metrics, function(beta) {
      iter.mats <- foreach::foreach(
        mat.i = mat.list,
        .verbose = debug,
        .export = c("mD1", "beta")
      ) %dopar% {
        mD.i <- replace.abundances(mD1, mat.i)
        microbData::beta.diversity(mD.i, metrics = beta, update.mD = F)[[1]] %>%
          as.matrix() %>%
          return()
      }
      { Reduce("+", iter.mats)/length(iter.mats) } %>%
        as.dist() %>%
        return()
    })
    mD1@Other.data[["Beta.metrics"]] <- beta.metrics
    mD1@Distance.matrices <- dist.mats
  }
  if (replace.with == "first") {
    mD1 %<>% replace.abundances(new.tbl = mat.list[[1]])
  } else {
    mD1 %<>% replace.abundances(new.tbl = { Reduce("+", mat.list)/length(mat.list) })
  }
  if (!quiet) { rlang::inform(paste("Number of samples dropped:", sum(sample.sums(mD) < rarefy.to))) }
  if (trim.features) {
    if (is.null(beta.metrics)) {
      rlang::inform(
        "Argument `beta.metrics' is NULL, so feature trimming is skipped as it could interfere with later beta-diversity estimation"
      )
    } else {
      if (!quiet) { rlang::inform(paste("Number of features dropped:", sum(feature.sums(mD1) == 0))) }
      mD1 %<>% drop.features(features = feature.sums(mD1) == 0)
    }
  }

  mD1 %<>% add.other.data(
    x = rarefy.to,
    name = "Abundances rarefied"
  )
  if (keep.tables) {
    mD1 %<>% add.other.data(
      x = mat.list,
      name = "Rarefied.tables"
    )
  }
  return(mD1)
}

####################################
#' @name subsample.features
#' @title Subsample Feature Counts
#' @description Subsample feature counts in a \code{microbData} object to the same depth.
#' @rdname transform_counts
#' @export

subsample.features <- function(
    mD,
    min.abund = 10000,
    exactly.to = NULL,
    trim.features = TRUE,
    user.seed = NULL,
    update.mD = TRUE,
    threads = 1,
    quiet = FALSE
) {
  dqrng.installed <- require(dqrng, quietly = TRUE)
  subsmpl.to <- ifelse(
    is.null(exactly.to),
    min(sample.sums(mD)[sample.sums(mD) >= min.abund]),
    exactly.to
  )
  if (!quiet) { rlang::inform(paste("Rarefying to:", subsmpl.to)) }
  if (subsmpl.to >= 1e5 && !dqrng.installed) {
    rlang::inform(
      "Your target subsampling depth is quite high (>=100,000 reads/sample). Consider installing the `dqrng' package to speed up subsampling"
    )
  }
  mD1 <- drop.samples(
    mD,
    samples = names(sample.sums(mD)[sample.sums(mD) < subsmpl.to])
  )
  if (is.null(user.seed)) {
    if (!quiet) {
      rlang::inform(paste("Random seed:", .Random.seed[1]))
    }
  } else {
    set.seed(user.seed)
  }
  if (threads == 1) {
    res <- apply(X = mD1@Abundances, MARGIN = 1, simplify = F, FUN = function(x) {
      if (dqrng.installed) {
        dqrng::dqsample(names(x), size = subsmpl.to, replace = T, prob = x) %>%
          table() %>%
          return()
      } else {
        sample(names(x), size = subsmpl.to, replace = T, prob = x) %>%
          table() %>%
          return()
      }
    })
  } else {
    cl <- parallel::makeCluster(threads, type = "FORK")
    res <- parallel::parApply(cl = cl, X = mD1@Abundances, MARGIN = 1, simplify = F, FUN = function(x) {
      if (dqrng.installed) {
        dqrng::dqsample(names(x), size = subsmpl.to, replace = T, prob = x) %>%
          table() %>%
          return()
      } else {
        sample(names(x), size = subsmpl.to, replace = T, prob = x) %>%
          table() %>%
          return()
      }
    })
    stopCluster(cl)
  }
  mD1@Abundances[mD1@Abundances >= 0] <- 0
  for (smpl in names(res)) {
    mD1@Abundances[smpl, names(res[[smpl]])] <- res[[smpl]]
  }
  if (!quiet) {
    rlang::inform(
      paste("Number of samples dropped:", sum(sample.sums(mD) < subsmpl.to))
    )
  }
  if (trim.features) {
    if (!quiet) {
      rlang::inform(
        paste("Number of features dropped:", sum(feature.sums(mD1) == 0))
      )
    }
    mD.res <- drop.features(mD1, features = feature.sums(mD1) == 0)
  } else {
    mD.res <- mD1
  }
  if (update.mD) {
    mD.res <- add.other.data(
      x = subsmpl.to,
      name = "Abundances rarefied",
      mD = mD1
    )
    return(mD.res)
  } else {
    return(mD.res@Abundances)
  }
}

####################################
#' @name center.log.ratio
#' @title Center-log Ratio Transform Feature Counts
#' @description Transform feature counts using the center log-ratio method. This function requires packages [CoDaSeq](https://github.com/ggloor/CoDaSeq) and \code{\link{zCompositions}}.
#' @rdname transform_counts
#' @export

center.log.ratio <- function(
    mD,
    min.abund = 10000,
    min.prop = 0.001,
    min.occur = 0,
    smpls.by.row = TRUE,
    method = "CZM",
    lab = 0,
    update.mD = TRUE,
    quiet = FALSE
) {

  require(CoDaSeq)
  require(zCompositions)
  filt.mat <- codaSeq.filter(
    mD@Abundances,
    min.reads = min.abund,

    min.prop = min.prop,
    min.occurrence = min.occur,
    samples.by.row = smpls.by.row
  ) %>%
    t()
  rep.mat <- try(cmultRepl(filt.mat, method = method, label = lab, z.warning = 1), silent = T)
  if ("try-class" %in% class(rep.mat)) {
    clr.mat <- codaSeq.clr(filt.mat)
  } else {
    clr.mat <- codaSeq.clr(rep.mat)
  }

  if (!quiet) {
    rlang::inform(paste("Number of samples dropped:", nsamples(mD) - nrow(clr.mat)))
    rlang::inform(paste("Number of features dropped:", nfeatures(mD) - ncol(clr.mat)))
  }
  if (update.mD) {
    mD.res <- keep.features(mD, features = colnames(clr.mat)) %>%
      keep.samples(samples = rownames(clr.mat))
    mD.res <- replace.abundances(mD = mD.res, new.tbl = clr.mat)
    mD.res <- add.other.data(x = "CLR", name = "Abundances transformed", mD = mD.res)
    return(mD.res)
  } else {
    return(clr.mat)
  }
}

####################################
#' @name variance.stabilize
#' @title Variance Stabilize Feature Counts
#' @description Transform feature counts using the variance-stabilized method. This function requires package \code{\link{DESeq2}}.
#' @rdname transform_counts
#' @export

variance.stabilize <- function(
    mD,
    design,
    update.mD = TRUE
) {
  vs.mat <- microbData2DEseq(mD, expt.design = design) %>%
    DESeq2::estimateSizeFactors() %>%
    DESeq2::estimateDispersions() %>%
    DESeq2::getVarianceStabilizedData() %>%
    t()
  if (update.mD) {
    mD@Abundances <- vs.mat
    mD <- add.other.data(
      x = "variance-stabilized",
      name = "Abundances transformed",
      mD = mD
    )
    return(mD)
  } else {
    return(vs.mat)
  }
}

####################################
#' @name relative.abundance
#' @title Relative Abundance of Feature Counts
#' @description Transform feature counts to relative abundance. This function is the equivalent of \code{transform_counts(mD = mD, f = function(x) x / sum(x))}, but because it's a popular method of transformation, it has its own wrapper here.
#' @rdname transform_counts
#' @export

relative.abundance <- function(mD, update.mD = TRUE) {
  mD@Abundances <- apply(
    X = mD@Abundances,
    MARGIN = 1,
    FUN = function(x) x / sum(x)
  ) %>%
    t()
  if (update.mD) {
    mD <- add.other.data(
      x = "relative",
      name = "Abundances transformed",
      mD = mD
    )
    return(mD)
  } else {
    return(mD@Abundances)
  }
}

