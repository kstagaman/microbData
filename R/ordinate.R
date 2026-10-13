####################################
# Ordination method registry: a dictionary of ordination methods, each describing
# what inputs it needs and how to compute it. Add new ordination methods here to
# make them available via `ordinate()`.
####################################

ordination.registry <- list(
  PCoA = list(
    requires.distance = TRUE,
    requires.formula = FALSE,
    supports.optimize = FALSE,
    fn = function(mD, dist.mat, formula, optimize, include.feature.scores, ...) {
      ape::pcoa(D = dist.mat, ...)
    }
  ),
  NMDS = list(
    requires.distance = TRUE,
    requires.formula = FALSE,
    supports.optimize = FALSE,
    fn = function(mD, dist.mat, formula, optimize, include.feature.scores, ...) {
      ord <- vegan::metaMDS(comm = dist.mat, ...)
      if (include.feature.scores) {
        vegan::sppscores(ord) <- mD@Abundances
      }
      ord
    }
  ),
  dbRDA = list(
    requires.distance = TRUE,
    requires.formula = TRUE,
    supports.optimize = TRUE,
    fn = function(mD, dist.mat, formula, optimize, include.feature.scores, ...) {
      mod.frm <- update.formula(formula, dist.mat ~ .)
      environment(mod.frm) <- environment()
      if (optimize) {
        meta.dt <- get.metadata(mD, return.as = "data.table")
        mod.data <- copy(meta.dt)[, .SD, .SDcols = c(all.vars(formula), mD@Sample.col)] %>%
          .[complete.cases(.)]
        dist.mat <- usedist::dist_subset(d = dist.mat, idx = mod.data[[mD@Sample.col]])
      } else {
        mod.data <- copy(get.metadata(mD, return.as = "data.table"))
      }
      if (include.feature.scores) {
        ord <- vegan::capscale(
          formula = mod.frm,
          data = mod.data,
          comm = get.abundances(mD),
          ...
        )
      } else {
        ord <- vegan::capscale(formula = mod.frm, data = mod.data, ...)
      }
      if (optimize) {
        ord <- vegan::ordistep(object = ord, direction = "both")
      }
      ord
    }
  ),
  PCA = list(
    requires.distance = FALSE,
    requires.formula = FALSE,
    supports.optimize = FALSE,
    fn = function(mD, dist.mat, formula, optimize, include.feature.scores, ...) {
      stats::prcomp(get.abundances(mD), ...)
    }
  )
)

#' @name ordinate
#' @title Create an Ordination
#' @description Create an ordination from a \code{microbData} object. Methods that operate on distance matrices (PCoA, NMDS, dbRDA) require the \code{microbData} to have Distance Matrices (can be created with \code{\link[microbData]{beta.diversity}}), and an ordination will be created for each matrix unless otherwise specified. Methods that operate directly on the abundance table (PCA) do not require distance matrices.
#' @param mD required; \code{microbData} object from which data will be read.
#' @param method character; the name of the ordination method to be used. Currently supported methods are "PCoA" (principal coordinate analysis; \code{\link[ape]{pcoa}}), "NMDS" (nonmetric multidimensional scaling; \code{\link[vegan]{metaMDS}}), "dbRDA" (distance-based redundancy analysis; \code{\link[vegan]{capscale}}), and "PCA" (principal components analysis; \code{\link[stats]{prcomp}}). "dbRDA" requires a formula, as it is a constrained ordination method. "PCA" operates on the abundance table directly rather than a distance matrix. Default is "dbRDA".
#' @param formula formula; in the form of e.g., \code{~ Var1 + Var2 * Var3}. This is only required for and used by the dbRDA method. Default is NULL.
#' @param optimize logical; should \code{\link[vegan]{ordistep}} be run on the dbRDA model and the results returned in place of the full model? Default is FALSE.
#' @param only.for character or integer; a vector of the names or the indices of the matrices contained in the Distance.matrices slot of the microbData. If NULL, will make an ordination for all matrices. Ignored for methods that do not operate on distance matrices (e.g., PCA). Default is NULL.
#' @param include.feature.scores logical; some ordination methods, like "NMDS" and "dbRDA" allow feature (often referred to a 'species') scores as well as site scores. If true, these will be calculated. Default is FALSE.
#' @param update.mD logical; should this function return a new \code{microbData} object with the new ordinations appended (unnamed) to the \code{Ordinations} slot and a tracking table appended to \code{Other.data$Ordination.log} (TRUE), or just a \code{list} with elements \code{Ordinations} (the unnamed list of results) and \code{Log} (the tracking table for just this call) (FALSE)? Default is TRUE.
#' @param ... other arguments passed to ordination methods. See their documentation for appropriate arguments.
#' @details Ordinations are stored as an \strong{unnamed} list in the \code{Ordinations} slot, since, e.g., multiple dbRDA ordinations with different formulas can be run on the same distance matrix, making a naming scheme based on method and distance metric alone ambiguous/colliding. Instead, a tracking table is kept in \code{Other.data$Ordination.log} with one row per ordination (in the same order as the \code{Ordinations} list) recording the \code{Index} into the \code{Ordinations} list, the \code{Method} used, the \code{Distance.metric} it was computed from (\code{NA} for methods like PCA that do not use one), the \code{Formula} used (\code{NA} if not applicable), and whether it was \code{Optimized} via \code{\link[vegan]{ordistep}}.
#' @seealso \code{\link[ape]{pcoa}}, \code{\link[vegan]{metaMDS}}, \code{\link[vegan]{capscale}}, \code{\link[stats]{prcomp}}, \code{\link[vegan]{sppscores}}
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
    rlang::abort("Object supplied to `mD` must be of class 'microbData`")
  }
  method <- rlang::arg_match(method, values = names(ordination.registry))
  ord.spec <- ordination.registry[[method]]

  if (ord.spec$requires.formula && is.null(formula)) {
    rlang::abort(
      sprintf("The argument supplied to `method` is '%s', but `formula` is NULL. Please supply a formula", method)
    )
  }
  if (!is.null(formula)) {
    if (!purrr::is_formula(formula)) {
      rlang::abort("The argument supplied to `formula` must be of class 'formula'")
    }
  }
  if (!is.logical(optimize)) {
    rlang::abort("The argument supplied to `optimize` must be of class 'logical'")
  }
  if (!ord.spec$supports.optimize & optimize) {
    rlang::warn(
      sprintf("The argument `optimize` is set to TRUE, but `method` ('%s') does not support it, this option will be ignored", method)
    )
  }
  formula.str <- if (!is.null(formula)) paste(deparse(formula), collapse = " ") else NA_character_

  if (ord.spec$requires.distance) {
    dist.mat.list <- get.distance.matrices(mD)
    if (is.null(dist.mat.list) || length(dist.mat.list) == 0) {
      rlang::abort(
        sprintf("The method '%s' operates on distance matrices, but `mD` has none. Use `beta.diversity()` to create them first.", method)
      )
    }
    if (!is.null(only.for)) {
      dist.mat.list <- dist.mat.list[only.for]
    }
    ord.list <- lapply(dist.mat.list, function(dist.mat) {
      ord.spec$fn(
        mD = mD,
        dist.mat = dist.mat,
        formula = formula,
        optimize = optimize,
        include.feature.scores = include.feature.scores,
        ...
      )
    })
    log.dt <- data.table(
      Method = method,
      Distance.metric = names(dist.mat.list),
      Formula = formula.str,
      Optimized = optimize
    )
  } else {
    if (!is.null(only.for)) {
      rlang::warn(
        sprintf("The argument `only.for` is ignored for method '%s' since it does not operate on distance matrices", method)
      )
    }
    ord.list <- list(
      ord.spec$fn(
        mD = mD,
        dist.mat = NULL,
        formula = formula,
        optimize = optimize,
        include.feature.scores = include.feature.scores,
        ...
      )
    )
    log.dt <- data.table(
      Method = method,
      Distance.metric = NA_character_,
      Formula = formula.str,
      Optimized = optimize
    )
  }
  ord.list <- unname(ord.list)

  if (update.mD) {
    start.idx <- length(mD@Ordinations)
    mD@Ordinations <- c(mD@Ordinations, ord.list)
    log.dt[, Index := start.idx + seq_len(.N)]
    setcolorder(log.dt, "Index")
    mD@Other.data[["Ordination.log"]] <- rbindlist(
      list(mD@Other.data[["Ordination.log"]], log.dt),
      use.names = TRUE
    )
    return(mD)
  } else {
    return(list(Ordinations = ord.list, Log = log.dt))
  }
}
