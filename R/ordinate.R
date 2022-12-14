#' @name ordinate
#' @title Create an Ordination
#' @description Create an ordination from a \code{microbData} object. The \code{microbData} must have Distance Matrices (can be created with \code{\link[microbData]{beta.diversity}}), and an ordination will be created for each matrix unless otherwise specified.
#' @param mD required; \code{microbData} object from which data will be read.
#' @param method character; the name of the ordination method to be used. Currently supported methods are "PCoA" (principal coordinate analysis; \code{\link[ape]{pcoa}}), "NMDS" (nonmetric multidimensional scaling; \code{\link[vegan]{metaMDS}}), and "dbRDA" (distance-based redundancy analysis; \code{\link[vegan]{capscale}}). The last requires a formula, as it is a constrained ordination method. Default is "dbRDA".
#' @param formula formula; in the form of e.g., \code{~ Var1 + Var2 * Var3}. This is only required for and used by the dbRDA method. Default is NULL.
#' @param optimize logical; should \code{\link[vegan]{ordistep}} be run on the dbRDA model and the results returned in place of the full model? Default is FALSE.
#' @param only.for character or integer; a vector of the names or the indices of the matrices contained in the Distance.matrices slot of the microbData. If NULL, will make an ordination for all matrices. Default is NULL.
#' @param include.feature.scores logical; some ordination methods, like "NMDS" and "dbRDA" allow feature (often referred to a 'species') scores as well as site scores. If true, these will be calculated. Default is FALSE.
#' @param update.mD logical; should this function return a new \code{microbData} object with the ordinations list (called "Ordinations") added to Other.data (TRUE) or just a list of the results (FALSE)? Default is TRUE.
#' @param ... other arguments passed to ordination methods. See their documentation for appropriate arguments.
#' @seealso \code{\link[ape]{pcoa}}, \code{\link[vegan]{metaMDS}}, \code{\link[vegan]{capscale}}, \code{\link[ape]{pcoa}},  \code{\link[vegan]{sppscores}}
#' @export

ordinate <- function(
    mD,
    method = "dbRDA",
    formula = NULL,
    optimize = FALSE,
    only.for = NULL,
    include.feature.scores = FALSE,
    update.mD = TRUE,
    ...
) {
  if (!{"microbData" %in% class(mD)}) {
    rlang::abort("Object supplied to `mD' must be of class 'microbData'")
  }
  method <- rlang::arg_match(method, values = c("PCoA", "NMDS", "dbRDA"))
  if (method == "dbRDA" & is.null(formula)) {
    rlang::abort("The argument supplied to `method' is 'dbRDA', but `formula' is NULL. Please supply a formula")
  }
  if (!is.null(formula)) {
    if (!purrr::is_formula(formula)) {
      rlang::abort("The argument supplied to `formula' must be of class 'formula'")
    }
  }
  if (!is.logical(optimize)) {
    rlang::abort("The argument supplied to `optimize' must be of class 'logical'")
  }
  if (method != "dbRDA" & optimize) {
    rlang::warn(
      "The argument `optimize' is set to TRUE, but `method' is not 'dbRDA', this option will be ignored"
      )
  }

  if (is.null(only.for)) {
    dist.mat.list <- mD@Distance.matrices
  } else {
    dist.mat.list <- mD@Distance.matrices[only.for]
  }
  ord.list <- lapply(dist.mat.list, function(dist.mat) {
    if (method == "PCoA") {
      ord <- ape::pcoa(D = dist.mat, ...)
    } else if (method == "NMDS") {
      ord <- vegan::metaMDS(comm = dist.mat, ...)
      if (include.feature.scores) {
        vegan::sppscores(ord) <- mD@Abundances
      }
    } else {
      mod.frm <- update.formula(formula, dist.mat ~ .)
      environment(mod.frm) <- environment()
      if (optimize) {
        mod.data <- copy(mD@Metadata)[, .SD, .SDcols = c(all.vars(formula), mD@Sample.col)] %>%
          .[complete.cases(.)]
        dist.mat <- usedist::dist_subset(d = dist.mat, idx = mod.data[[mD@Sample.col]])
      } else {
        mod.data <- copy(mD@Metadata)
      }
      if (include.feature.scores) {
        ord <- vegan::capscale(
          formula = mod.frm,
          data = mod.data,
          comm = mD@Abundances,
          ...
        )
      } else {
        ord <- vegan::capscale(formula = mod.frm, data = mod.data, ...)
      }
      if (optimize) {
        ord <- vegan::ordistep(object = ord, direction = "both")
      }
    }
    return(ord)
  })
  names(ord.list) <- paste0(names(ord.list), "_", method)
  if (update.mD) {
    add.other.data(
      x = ord.list, name = ifelse(optimize, "Optimized.ordinations", "Ordinations"), mD = mD
      ) %>%
      return()
  } else {
    return(ord.list)
  }
}
