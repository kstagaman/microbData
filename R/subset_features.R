#' @title Keep, Drop, or Filter Features in a \code{microbData} object
#' @description Functions to specify specific features to keep/drop by name or to filter features based on variable(s) in the Metadata
#' @aliases keep.features
#' @aliases drop.features
#' @aliases filter.features
#' @aliases remove.eukarya
#' @param mD required; the \code{microbData} object to be updated.
#' @param features character or logical; a vector of feature names to keep or drop. Also can be a logical vector that is either named or in the exact same order as those in the Feature Names slot.
#' #' @param track logical; should the names of the dropped samples be added to Other.data? Default is TRUE.
#' @param ... logical expression(s) by which to filter features based on the Features table, e.g. \code{Kingdom == "Bacteria" & Family != "Mitochondria"}.
#' @seealso \code{\link[microbData]{microbData}}, \code{\link[ape]{drop.tip}}

####################################
#' @name keep.features
#' @title Keep Features
#' @description Create a new \code{microbData} objects with just the specified features
#' @rdname subset.features
#' @export

keep.features <- function(mD, features, track = TRUE) {
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
  if (track) {
    mD <- add.other.data(
      mD = mD,
      name = "Dropped.features",
      x = mD@Feature.names[!{mD@Feature.names %in% to.keep}]
    )
  }
  keep.order <- order(colSums(mD@Abundances[, to.keep]), decreasing = T)
  mD@Abundances <- mD@Abundances[, to.keep[keep.order]]
  mD@Features <- copy(mD@Features)[to.keep[keep.order]]
  if (!is.null(mD@Phylogeny)) {
    to.drop <- mD@Feature.names[!{mD@Feature.names %in% to.keep}]
    mD@Phylogeny <- ape::drop.tip(mD@Phylogeny, tip = to.drop)
  }
  mD@Feature.names <- to.keep[keep.order]
  return(mD)
}

####################################
#' @name drop.features
#' @title Drop Features
#' @description Create a new \code{microbData} objects without the specified features
#' @rdname subset.features
#' @export

drop.features <- function(mD, features, track = TRUE) {
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
  if (track) {
    mD <- add.other.data(mD = mD, name = "Dropped.features", x = features)
  }
  mD@Abundances <- mD@Abundances[, to.keep]
  mD@Features <- copy(mD@Features)[to.keep]
  if (!is.null(mD@Phylogeny)) {
    mD@Phylogeny <- ape::drop.tip(mD@Phylogeny, tip = to.drop)
  }
  mD@Feature.names <- to.keep
  return(mD)
}

####################################
#' @name filter.features
#' @title Filter Features
#' @description Create a new \code{microbData} objects with features that match the filter
#' @rdname subset.features
#' @export

filter.features <- function(mD, ..., track = TRUE) {
  mD.return <- copy(mD@Features)[...][[mD@Feature.col]] %>%
    keep.features(mD = mD, features = ., track = FALSE)
  if (track) {
    mD.return <- add.other.data(mD = mD.return, name = "Feature.filter", x = deparse(substitute(...)))
  }
  return(mD.return)
}

####################################
#' @name remove.eukarya
#' @title Remove ASVs with eukaryotic taxonomic assignments
#' @description Create a new \code{microbData} objects without ASVs assignment to eukaryotic taxonomy
#' @rdname subset.features
#' @export

remove.eukarya <- function(mD, track = TRUE) {
  euks <- "(Kingdom == 'Bacteria' | Kingdom == 'Archaea') & Order != 'Chloroplast' & Family != 'Mitochondria'"
  copy(mD@Features)[eval(parse(text = euks))][[mD@Feature.col]] %>%
    keep.features(mD = mD, features = .) %>%
    add.other.data(name = "Feature.filter", x = "removed Eukarya") %>%
    return()
}
