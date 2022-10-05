#' @name beta.dispersion
#' @title Estimate Beta-dispersion
#' @description Utilizes \code{\link[vegan]{betadisper}} to calculate beta-dispersion (distance to centroid) for each sample by a group.
#' @param mD required; the `microbData` object, with distance matrices (see \code{\link[microbData]{beta.diversity}}) from which to estimate beta-dispersion.
#' @param group required; a character string indicating by which column in the Metadata table the samples should be grouped.
#' @param metrics character; if NULL, beta-dispersion will be calculated for each distance matrix present in the `microbData` object. If metrics are provided, beta-dispersion will only be calculated for those metrics. Default is NULL.
#' @param type character; the type of analysis to perform. Use the spatial median ("median") or the group centroid ("centroid")? Default is "median". (Directly from \code{\link[vegan]{betadisper}}).
#' @param bias.adjust logical: adjust for small sample bias in beta diversity estimates? Default is FALSE.
#' @param sqrt.dist logical; take square root of dissimilarities. This often euclidifies dissimilarities. Default is FALSE. (Directly from \code{\link[vegan]{betadisper}}).
#' @param add logical; add a constant to the non-diagonal dissimilarities such that all eigenvalues are non-negative in the underlying Principal Co-ordinates Analysis (see \code{\linkl[vegan]{wcmdscale}} for details). Choice "lingoes" (or TRUE) use the recommended method of Legendre & Anderson (1999: “method 1”) and "cailliez" uses their “method 2”. Default is FALSE. (Directly from \code{\link[vegan]{betadisper}}).
#' @param update.mD logical; should this function return a new microbData object with the beta-dispersion results added to the Metadata table and add `metrics` to Other.data (TRUE) or just the results of the beta-dispersion estimation (FALSE)? Default is TRUE.
#' @details This function utilizes \code{\link[vegan]{betadisper}} to estimate beta-dispersion (i.e., distance from spatial median or group centroid in ordination space) by categorical group. See `?vegan::betadisper` for more details about the procedure that function employs.
#' @returns If `update.mD` is TRUE, a `microbData` object with estimated beta-dispersion scores added to the Metadata table and a list of the "betadisper" objects added to the Other.data slot.
#' @returns If `update.mD` is FALSE, a list containing:
#'
#' \item{Distances}{a table of the estimated beta-dispersion scores for each sample}
#' \item{Betadispers}{a list of the "betadisper" objects}
#'
#' @seealso \code{\link[vegan]{betadisper}}
#' @export

beta.dispersion <- function(
    mD,
    group,
    metrics = NULL,
    type = c("median", "centroid"),
    bias.adjust = FALSE,
    sqrt.dist = FALSE,
    add = FALSE,
    update.mD = TRUE
) {
  if (is.null(mD@Distance.matrices)) {
    rlang::abort("The Distance.matrices slot of the supplied microbData object must not by empty. See ?microbData::beta.diversity to generate distance matrices")
  }
  if (is.null(metrics)) {
    metrics <- mD@Other.data$Beta.metrics
  }
  type <- rlang::arg_match(type, values = c("median", "centroid"))
  bad.logical <- function(arg.name) {
    rlang::abort(paste0("Argument `", arg.name, "' must be logical or coercible to logical."))
  }
  bias.adjust <- as.logical(bias.adjust)
  if (is.na(bias.adjust)) { bad.logical("bias.adjust") }
  sqrt.dist <- as.logical(sqrt.dist)
  if (is.na(sqrt.dist)) { bad.logical("sqrt.dist") }
  add <- as.logical(add)
  if (is.na(add)) { bad.logical("add") }
  update.mD <- as.logical(update.mD)
  if (is.na(update.mD)) { bad.logical("update.mD") }
  res.dt <- mD@Metadata
  bd.list <- NULL
  for (beta in metrics) {
    bd <- vegan::betadisper(
      d = mD@Distance.matrices[[beta]],
      group = res.dt[[group]],
      type = type,
      bias.adjust = bias.adjust,
      sqrt.dist = sqrt.dist,
      add = add
    )
    res.dt <- data.table(Smpl = names(bd$distances), Dtc = bd$distances) %>%
      set_names(c(mD@Sample.col, paste0(beta, ".dispersion"))) %>%
      setkeyv(mD@Sample.col) %>%
      merge(x = res.dt, y = .)
    bd.list[[paste(beta, group, type, sep = "_")]] <- bd
  }
  if (update.mD) {
    replace.metadata(mD, res.dt) %>%
      add.other.data(name = "Beta.dispersion", x = bd.list) %>%
      return()
  } else {
    return(list(Distances = res.dt, Betadispers = bd.list))
  }
}
