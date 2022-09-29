#' @name transform.counts
#' @title Transform Feature Counts
#' @description Functions for transforming feature counts in a microbData object.
#' @aliases transform.counts
#' @aliases rarefy
#' @aliases center.log.ratio
#' @aliases variance.stabilize
#' @aliases relative.abundance
#' @param mD required; the microbData object with feature counts to be transformed.
#' @param update.mD logical; should this function return a new microbData object with the transformed abundances replacing the original abundances (TRUE) or just the resulting transformed abundance table (FALSE)? Default is TRUE.
#' @param quiet logical; should informational output (not warnings or errors) be suppressed? Default is FALSE.
#' @param f required for `transform.counts`; a function that can be applied to the samples (rows) of the Abundance table.
#' @param min.abund integer; for `rarefy`, if not NULL, samples will be rarefied to the same depth as the lowest sample abundance equal to or greater than this integer. E.g., if `min.abund = 10000` and there four samples with 8000, 11000, 12000, and 13000 total reads, respectively, the first sample will be dropped, and last three samples will be rarefied to 11000 reads. For `center.log.ratio`, the interger to pass to the `min.reads` argument in \code{\link[CoDaSeq]{codaSeq.filter}}. Default is 1e4.
#' @param exactly.to integer; for `rarefy`, if not NULL, samples will be rarefied to exactly this integer. This argument supersedes `min.abund`. Samples with total reads lower than this number will be dropped. Default is NULL.
#' @param trim.features logical; for `rarefy`, should features that are no longer present in any samples after rarefaction be dropped from the microbData object? Default is TRUE.
#' @param user.seed integer; for `rarefy`, a user-supplied random seed to make the random subsampling process repeatable. If NULL, will just use the default \code{\link{.Random.seed}}, which will be reported at the end. Default is NULL
#' @param min.prop numeric; for `center.log.ratio`, the minimum proportional abundance of a read in any sample. (See \code{\link[CoDaSeq]{codaSeq.filter}}). Default is 0.001.
#' @param min.occur numeric; for `center.log.ratio`, the minimum fraction of non-0 reads for each variable in all samples. (See \code{\link[CoDaSeq]{codaSeq.filter}}). Default is 0.
#' @param smpls.by.row logical; for `center.log.ratio`, TRUE if rows contain samples, FALSE if rows contain variables. (See \code{\link[CoDaSeq]{codaSeq.filter}}). Default is TRUE.
#' @param method character; for `center.log.ratio`. Choose one of: count zero multiplicative ("CZM"); geometric Bayesian multiplicative ("GBM"); square root BM ("SQ"); Bayes-Laplace BM ("BL"); user-specified hyper-parameters ("user"). (See \code{\link[zCompositions]{cmultRepl}}). Default is "CZM"
#' @param lab numeric or character; for `center.log.ratio`, a unique label used to denote count zeros in X. (See \code{\link[zCompositions]{cmultRepl}}). Default is 0.
#' @param design formula or matrix; required for `variance.stabilize`, from \code{\link[DESeq2]{DESeqDataSetFromMatrix}}: the formula expresses how the counts for each gene depend on the variables in colData. Many R formula are valid, including designs with multiple variables, e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment. See results for a variety of designs and how to extract results tables. By default, the functions in this package will use the last variable in the formula for building results tables and plotting. ~ 1 can be used for no design, although users need to remember to switch to another design for differential testing.
#' @seealso \code{\link{}}
#' @export

####################################
#' @name transform.counts
#' @title Transform Feature Counts
#' @description Supply a custom function with which to transform feature counts in a microbData object.
#' @rdname transform.counts
#' @export

transform.counts <- function(mD, f, update.mD = TRUE) {
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
#' @description Rarefy feature counts in a microbData object to the same depth.
#' @rdname transform.counts
#' @export

rarefy <- function(
    mD,
    min.abund = 1e4,
    exactly.to = NULL,
    trim.features = TRUE,
    user.seed = NULL,
    update.mD = TRUE,
    quiet = FALSE
) {
  if (is.null(exactly.to)) {
    rarefy.to <- min(sample.sums(mD)[sample.sums(mD) >= min.abund])
    if (!quiet) { rlang::inform(paste("Rarefying to:", rarefy.to)) }
  } else {
    rarefy.to <- exactly.to
  }
  mD.rar <- drop.samples(
    mD,
    samples = names(sample.sums(mD)[sample.sums(mD) < rarefy.to])
  )

  if (is.null(user.seed)) {
    if (!quiet) {
      rlang::inform(paste("Random seed:", .Random.seed[1]))
    }
  } else {
    set.seed(user.seed)
  }
  res <- apply(
    X = mD.rar@Abundances,
    MARGIN = 1,
    FUN = function(x) {
      sample(names(x), size = rarefy.to, replace = T, prob = x) %>%
        table()
    }
  )
  mD.rar@Abundances[mD.rar@Abundances >= 0] <- 0
  for (smpl in names(res)) {
    mD.rar@Abundances[smpl, names(res[[smpl]])] <- res[[smpl]]
  }
  if (!quiet) {
    rlang::inform(
      paste("Number of samples dropped:", sum(sample.sums(mD) < rarefy.to))
    )
  }
  if (trim.features) {
    if (!quiet) {
      rlang::inform(
        paste("Number of features dropped:", sum(feature.sums(mD.rar) == 0))
      )
    }
    mD.res <- drop.features(mD.rar, features = feature.sums(mD.rar) == 0)
  } else {
    mD.res <- mD.rar
  }
  if (update.mD) {
    mD.res <- add.other.data(x = rarefy.to, name = "Abundances rarefied", mD = mD)
    return(mD.res)
  } else {
    retunr(mD.res@Abundances)
  }
}

####################################
#' @name center.log.ratio
#' @title Center-log Ratio Transform Feature Counts
#' @description Transform feature counts using the center log-ratio method. This function requires packages [CoDaSeq](https://github.com/ggloor/CoDaSeq) and \code{\link{zCompositions}}.
#' @rdname transform.counts
#' @export

center.log.ratio <- function(
    mD,
    min.abund = 1e4,
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

  clr.mat <- codaSeq.filter(
    mD@Abundances,
    min.reads = min.abund,
    min.prop = min.prop,
    min.occurrence = min.occur,
    samples.by.row = smpls.by.row
  ) %>%
    t() %>%
    cmultRepl(
      method = method,
      label = lab,
      z.warning = 1
    ) %>% # replace 0 values with an estimate
    codaSeq.clr()
  if (!quiet) {
    rlang::inform(
      paste("Number of samples dropped:", nsamples(mD) - nrow(clr.mat))
    )
    rlang::inform(
      paste("Number of features dropped:", nfeatures(mD) - ncol(clr.mat))
    )
  }
  if (update.mD) {
    mD.res <- keep.features(mD, features = colnames(clr.mat)) %>%
      keep.samples(samples = rownames(clr.mat))
    mD.res <- add.other.data(x = "CLR", name = "Abundances transformed", mD = mD)
    return(mD.res)
  } else {
    return(clr.mat)
  }
}

####################################
#' @name variance.stabilize
#' @title Variance Stabilize Feature Counts
#' @description Transform feature counts using the variance-stabilized method. This function requires package \code{\link{DESeq2}}.
#' @rdname transform.counts
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
#' @description Transform feature counts to relative abundance. This function is the equivalent of `transform.counts(mD = mD, f = function(x) x / sum(x))`, but because it's a popular method of transformation, it has its own wrapper here.
#' @rdname transform.counts
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
