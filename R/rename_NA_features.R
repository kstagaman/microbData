#' @name rename.NA.features
#' @title
#' @description
#' @param mD required; microbData object to read from and write to.
#' @param force.split logical; should the lowest feature assignment (e.g. ASV) be appended to replacement assigment? If TRUE, this will prevent features that have NAs at the same assigment level from receiving identical assigments (keeping separate during future glomming). If FALSE, NAs will be replaced by a string including the previous level assigment and the current assignment level. Default is FALSE.
#' @param level.order character; a vector specifying the hierarchical ordering of the assignment levels. Must be identical to, or a subset of the column names in the Feature Assignments table. If NULL, the order will be inferred as Feature column lowest and the rest in order as the table columns. Default is NULL.
#' @seealso \code{\link{}}
#' @export
#' @examples
#'


rename.NA.features <- function(mD, force.split = FALSE, level.order = NULL) {
  feature.tbl <- copy(mD@Features)
  if (!is.null(level.order)) {
    if (!all(level.order %in% names(feature.tbl))) {
      rlang::abort(
        "Levels supplied to `level.order' must be the same as or a subset of the column names in the Feature Table of the microbData object."
      )
    } else {
      assign.lvls <- level.order
    }
  } else {
    assign.lvls <- c(names(feature.tbl)[2:ncol(feature.tbl)], mD@Feature.col)
  }
  feature.tbl <- copy(feature.tbl)[, ..assign.lvls]
  for (col in 1:{ncol(feature.tbl) - 1}) {
    taxlevel <- names(feature.tbl)[col]
    curr.col <- feature.tbl[[col]]
    to.replace <- is.na(curr.col) |
      grepl("Unknown_Family", curr.col) |
      grepl("Incertae_Sedis", curr.col)
    if (col == 1) {
      if (force.split) {
        replacing.feature.names <- mD@Feature.names[to.replace]
        curr.col[to.replace] <- paste0("NA_", taxlevel, "_", replacing.feature.names)
      } else {
        curr.col[to.replace] <- paste0("NA_", taxlevel)
      }
    } else {
      prev.col <- feature.tbl[[col - 1]]
      if (force.split) {
        replacing.feature.names <- mD@Feature.names[to.replace]
        curr.col[to.replace] <- do.call(
          paste,
          c(
            cbind(
              str_remove_all(prev.col[to.replace], paste(replacing.feature.names, collapse = "|")),
              replacing.feature.names,
            ),
            sep = "_"
          )
        )
      } else {
        curr.col[to.replace] <- paste0(prev.col[to.replace], "_", taxlevel)
      }
    }
    feature.tbl[[col]] <- curr.col
  }
  mD@Features <- feature.tbl

}
