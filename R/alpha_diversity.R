#' @name alpha.diversity
#' @title Estimate Alpha-diversity
#' @description This function take a \code{microbData} object and estimates alpha-diversity (see `metrics`) from the abundance table.
#' @param mD required; \code{microbData} object from which data will be read.
#' @param metrics character; a vector of alpha-diveristy metrics to be estimated from the Abundances table in the \code{microbData} object. Supported metrics include "Richness" (a.k.a observed feature counts),"Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Phylogenetic" (requires a Phylogeny to be present in the \code{microbData} object), and "Fisher".
#' @param update.mD logical; should this function return a new \code{microbData} object with the alpha-diversity results added to the Metadata table and add \code{metrics} to Other.data (TRUE) or just the results of the alpha-diversity estimation (FALSE)? Default is TRUE.
#' @details This function utilizes \code{\link[vegan]{diversity}}, \code{\link[vegan]{estimateR}}, and/or \code{\link[picante]{pd}} to estimate common alpha-diversity metrics for microbial data.
#' @returns A \code{microbData} object with estimated alpha-diversity scores added to the Metadata table, OR a table of the estimated alpha-diversity scores for each sample.
#' @seealso \code{\link[vegan]{diversity}}, \code{\link[vegan]{estimateR}}, \code{\link[picante]{pd}}
#' @examples
#' data("mD_rar")
#' mD.new <- alpha.diversity(mD.rar)
#' print(mD.new)
#' @export

alpha.diversity <- function(
    mD,
    metrics = c("Chao1", "Shannon", "Simpson", "Phylogenetic"),
    update.mD = TRUE
) {
  diversity.metrics <- c("Shannon", "Simpson", "InvSimpson")
  estimateR.metrics <- c("ACE", "Chao1", "Richness")
  pd.metric <- c("Phylogenetic")
  fisherfit.metric <- c("Fisher")
  all.metrics <- c(diversity.metrics, estimateR.metrics, pd.metric, fisherfit.metric)
  if (!all(metrics %in% all.metrics)) {
    rlang::abort(
      paste(
        "Metrics supplied to `metrics' must be of the following:",
        paste(all.metrics, collapse = ", ")
      )
    )
  }
  results <- NULL
  if (any(diversity.metrics %in% metrics)) {
    div.mets <- metrics[metrics %in% diversity.metrics]
    res.dt <- NULL
    for (metric in div.mets) {
      dt <- vegan::diversity(mD@Abundances, index = tolower(metric)) %>%
        as.data.table(keep.rownames = mD@Sample.col)
      names(dt)[2] <- metric
      if (is.null(res.dt)) {
        res.dt <- dt
      } else {
        res.dt <- merge(res.dt, dt, by = mD@Sample.col)
      }
    }
    if (is.null(results)) {
      results <- res.dt
    } else {
      results <- merge(results, res.dt, by = mD@Sample.col)
    }
  }
  if (any(estimateR.metrics %in% metrics)) {
    esr.mets <- metrics[metrics %in% estimateR.metrics]
    res.dt <- vegan::estimateR(mD@Abundances) %>%
      t() %>%
      as.data.table(keep.rownames = mD@Sample.col)
    res.dt[, `:=`(se.chao1 = NULL, se.ACE = NULL)]
    names(res.dt)[2:4] <- c("Richness", "Chao1", "ACE")
    keep.cols <- c(mD@Sample.col, esr.mets)
    res.dt <- res.dt[, ..keep.cols]
    if (is.null(results)) {
      results <- res.dt
    } else {
      results <- merge(results, res.dt, by = mD@Sample.col)
    }
  }
  if (pd.metric %in% metrics) {
    res.dt <- picante::pd(samp = mD@Abundances, tree = mD@Phylogeny) %>%
      as.data.table(keep.rownames = mD@Sample.col)
    res.dt[, SR := NULL]
    names(res.dt)[2] <- pd.metric
    if (is.null(results)) {
      results <- res.dt
    } else {
      results <- merge(results, res.dt, by = mD@Sample.col)
    }
  }
  if (fisherfit.metric %in% metrics) {
    res.dt <- vegan::fisher.alpha(mD@Abundances) %>%
      as.data.table(keep.rownames = mD@Sample.col)
    names(res.dt)[2] <- fisherfit.metric
    if (is.null(results)) {
      results <- res.dt
    } else {
      results <- merge(results, res.dt, by = mD@Sample.col)
    }
  }
  if (update.mD) {
    mD@Metadata <- merge(mD@Metadata, results, by = mD@Sample.col)
    mD@Other.data[["Alpha.metrics"]] <- metrics
    return(mD)
  } else {
    return(results)
  }
}
