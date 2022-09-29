#' @name subset.features
#' @title Keep, Drop, or Filter Features in a microbData object
#' @description Functions to specify specific features to keep/drop by name or to filter features based on variable(s) in the Metadata
#' @aliases keep.features
#' @aliases drop.features
#' @aliases filter.features
#' @param mD required; the microbData object to be updated.
#' @param features character or logical; a vector of feature names to keep or drop. Also can be a logical vector that is either named or in the exact same order as those in the Feature Names slot.
#' @param ... logical expression(s) by which to filter features based on the Features table, e.g. `Kingdom == "Bacteria" & Class != "Mitochondria"`.
#' @seealso \code{\link[microbData]{microbData}}


#' @name keep.features
#' @title Keep Features
#' @description Create a new microbData objects with just the specified features
#' @rdname subset.features
#' @export

keep.features <- function(mD, features) {
  if (is.logical(features)) {
    if (is.null(names(features))) {
      if (length(features) != nfeatures(mD)) {
        rlang::abort(
          "The length of the supplied logical vector does not match the number of features in the supplied microbData object"
        )
      } else {
        to.keep <- mD@Feature.names[features]
      }
    } else {
      to.keep <- names(features)[features]
    }
  } else {
    to.keep <- features
  }
  if (!all(to.keep %in% mD@Feature.names)) {
    rlang::abort("One or more of the features supplied to `features' is not in the microbData Feature Names")
  }
  mD@Abundances <- mD@Abundances[, to.keep]
  mD@Features <- copy(mD@Features)[to.keep]
  if (!is.null(mD@Phylogeny)) {
    to.drop <- to.keep[!{mD@Feature.names %in% to.keep}]
    mD@Phylogeny <- ape::drop.tip(mD@Phylogeny, tip = to.drop)
  }
  mD@Feature.names <- to.keep
  return(mD)
}

#' @name drop.features
#' @title Drop Features
#' @description Create a new microbData objects without the specified features
#' @rdname subset.features
#' @export

drop.features <- function(mD, features) {
  if (is.logical(features)) {
    if (is.null(names(features))) {
      if (length(features) != nfeatures(mD)) {
        rlang::abort(
          "The length of the supplied logical vector does not match the number of features in the supplied microbData object"
        )
      } else {
        to.keep <- mD@Feature.names[!features]
        to.drop <- mD@Feature.names[features]
      }
    } else {
      to.keep <- names(features)[!features]
      to.drop <- names(features)[features]
    }
  } else {
    to.keep <- mD@Feature.names[!{mD@Feature.names %in% features}]
    to.drop <- mD@Feature.names[mD@Feature.names %in% features]
  }
  if (!all(to.keep %in% mD@Feature.names)) {
    rlang::abort("One or more of the features to be kept is not in the microbData Feature Names")
  }
  mD@Abundances <- mD@Abundances[, to.keep]
  mD@Features <- copy(mD@Features)[to.keep]
  if (!is.null(mD@Phylogeny)) {
    mD@Phylogeny <- ape::drop.tip(mD@Phylogeny, tip = to.drop)
  }
  mD@Feature.names <- to.keep
  return(mD)
}

#' @name filter.features
#' @title Filter Features
#' @description Create a new microbData objects with features that match the filter
#' @rdname subset.features
#' @export

filter.features <- function(mD, ...) {
  copy(mD@Features)[...][[mD@Feature.col]] %>%
    keep.features(mD = mD, features = .) %>%
    return()
}

