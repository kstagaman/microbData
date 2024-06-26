#' @name update_feature_names
#' @title Update Feature Names
#' @description A simple way to replace the feature names in a \code{microbData} object.
#' @param mD required; \code{microbData} object to read from and write to.
#' @param new.names required; a character vector of new feature names. This vector must be the same length as the original Feature Names and in the same order.
#' @export

update_feature_names <- function(mD, new.names) {
  if (!is.null(mD@Feature.names)) {
    if (length(mD@Feature.names) != length(new.names)) {
      rlang::abort(
        "The length of the vector `new.names' is not equal to the length of the original Feature Names."
      )
    }
    old.names <- copy(mD@Feature.names)
    names.dt <- data.table(
      Old = old.names,
      New = new.names
    ) %>% setkey(Old)
    colnames(mD@Abundances) <- names.dt[colnames(mD@Abundances)]$New
    mD@Abundances <- mD@Abundances[, sort(colnames(mD@Abundances))]
    mD@Assignments[[mD@Feature.col]] <- names.dt[mD@Assignments[[mD@Feature.col]]]$New
    setkeyv(mD@Assignments, mD@Feature.col)
    mD@Assignments[new.names]
    mD@Phylogeny$tip.label <- names.dt[mD@Phylogeny$tip.label]$New
    mD@Feature.names <- names.dt[mD@Feature.names]$New
  } else {
    mD@Feature.names <- new.names
  }
  return(mD)
}
