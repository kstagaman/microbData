#' @name rename.NA.assignments
#' @title Rename NA Feature Assignments
#' @description This function searches through the Assignments table of a \code{microbData} object and replaces any NA assignments with a more descriptive term based on the higher order assignment provided. This prevents unrelated features lacking labels from getting grouped together simply because they're labelled "NA".
#' @param mD required; \code{microbData} object to read from and write to.
#' @param force.split logical; should the lowest feature assignment (e.g. ASV) be appended to replacement assignment? If TRUE, this will prevent features that have NAs at the same assignment level from receiving identical assignments (keeping separate during future glomming). If FALSE, NAs will be replaced by a string including the previous level assignment and the current assignment level. Default is FALSE.
#' @param level.order character; a vector specifying the hierarchical ordering of the assignment levels. Must be identical to, or a subset of the column names in the Feature Assignments table. If NULL, the order will be inferred as Feature column lowest and the rest in order as the table columns. Default is NULL.
#' @export

rename.NA.assignments <- function(mD, force.split = FALSE, level.order = NULL) {
  assignment.tbl <- copy(mD@Assignments)
  if (!is.null(level.order)) {
    if (!all(level.order %in% names(assignment.tbl))) {
      rlang::abort(
        "Levels supplied to `level.order' must be the same as or a subset of the column names in the Feature Table of the microbData object."
      )
    } else {
      assign.lvls <- level.order
    }
  } else {
    assign.lvls <- c(names(assignment.tbl)[2:ncol(assignment.tbl)], mD@Feature.col)
  }
  assignment.tbl <- copy(assignment.tbl)[, ..assign.lvls]
  for (col in 1:{ncol(assignment.tbl) - 1}) {
    taxlevel <- names(assignment.tbl)[col]
    curr.col <- assignment.tbl[[col]]
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
      prev.col <- assignment.tbl[[col - 1]]
      if (force.split) {
        replacing.feature.names <- mD@Feature.names[to.replace]
        curr.col[to.replace] <- do.call(
          paste,
          c(
            cbind(
              str_remove_all(
                prev.col[to.replace],
                paste(replacing.feature.names, collapse = "|")
              ),
              replacing.feature.names,
            ),
            sep = "_"
          )
        )
      } else {
        curr.col[to.replace] <- paste0(prev.col[to.replace], "_", taxlevel)
      }
    }
    assignment.tbl[[col]] <- curr.col
  }
  mD@Assignments <- assignment.tbl
  return(mD)
}
