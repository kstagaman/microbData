#' @title Keep, Drop, or Filter Samples in a \code{microbData} object
#' @description Functions to specify specific samples to keep/drop by name or to filter samples based on variable(s) in the Metadata
#' @aliases keep.samples
#' @aliases drop.samples
#' @aliases filter.samples
#' @param mD required; the \code{microbData} object to be updated.
#' @param samples character or logical; a vector of sample names to keep or drop. Also can be a logical vector that is either named or in the exact same order as those in the Sample Names slot.
#' @param track logical; should the names of the dropped samples be added to Other.data? Default is TRUE.
#' @param ... logical expression(s) by which to filter samples based on Metadata, e.g. \code{Weight > 10 & Genotype != "WT"}.
#' @seealso \code{\link[microbData]{microbData}}, \code{\link[ape]{drop.tip}}, \code{\link[usedist]{dist_subset}}


#' @name keep.samples
#' @title Keep Samples
#' @description Create a new \code{microbData} objects with just the specified samples
#' @rdname subset.samples
#' @export

keep.samples <- function(mD, samples, track = TRUE) {
  if (is.logical(samples)) {
    if (is.null(names(samples))) {
      if (length(samples) != nsamples(mD)) {
        rlang::abort(
          "The length of the supplied logical vector does not match the number of samples in the supplied microbData object"
        )
      } else {
        to.keep <- mD@Sample.names[samples]
      }
    } else {
      to.keep <- names(samples)[samples]
    }
  } else {
    to.keep <- samples
  }
  if (!all(to.keep %in% mD@Sample.names)) {
    rlang::abort("One or more of the samples supplied to `samples' is not in the microbData Sample Names")
  }
  if (track) {
    mD <- add.other.data(
      mD = mD,
      name = "Dropped.samples",
      x = mD@Sample.names[!{mD@Sample.names %in% to.keep}]
    )
  }
  mD@Metadata <- copy(mD@Metadata)[to.keep]
  mD@Abundances <- mD@Abundances[to.keep, ]
  if (!is.null(mD@Distance.matrices)) {
    if (class(mD@Distance.matrices) == "dist") {
      mD@Distance.matrices <- usedist::dist_subset(mD@Distance.matrices, to.keep)
    } else {
      mD@Distance.matrices <- lapply(mD@Distance.matrices, function(dist.mat) {
        usedist::dist_subset(dist.mat, to.keep) %>%
          return()
      })
    }
  }
  mD@Sample.names <- to.keep
  return(mD)
}

#' @name drop.samples
#' @title Drop Samples
#' @description Create a new \code{microbData} objects without the specified samples
#' @rdname subset.samples
#' @export

drop.samples <- function(mD, samples, track = TRUE) {
  if (is.logical(samples)) {
    if (is.null(names(samples))) {
      if (length(samples) != nsamples(mD)) {
        rlang::abort(
          "The length of the supplied logical vector does not match the number of samples in the supplied microbData object"
        )
      } else {
        to.keep <- mD@Sample.names[!samples]
      }
    } else {
      to.keep <- names(samples)[!samples]
    }
  } else {
    to.keep <- mD@Sample.names[!{mD@Sample.names %in% samples}]
  }
  if (!all(to.keep %in% mD@Sample.names)) {
    rlang::abort("One or more of the samples to be kept is not in the microbData Sample Names")
  }
  if (track) {
    mD <- add.other.data(mD = mD, name = "Dropped.samples", x = samples)
  }
  mD@Metadata <- copy(mD@Metadata)[to.keep]
  mD@Abundances <- mD@Abundances[to.keep, ]
  if (!is.null(mD@Distance.matrices)) {
    if (class(mD@Distance.matrices) == "dist") {
      mD@Distance.matrices <- usedist::dist_subset(mD@Distance.matrices, to.keep)
    } else {
      mD@Distance.matrices <- lapply(mD@Distance.matrices, function(dist.mat) {
        usedist::dist_subset(dist.mat, to.keep) %>%
          return()
      })
    }
  }
  mD@Sample.names <- to.keep
  return(mD)
}

#' @name filter.samples
#' @title Filter Samples
#' @description Create a new \code{microbData} objects with samples that match the filter
#' @rdname subset.samples
#' @export

filter.samples <- function(mD, ..., track = TRUE) {
  mD.return <- copy(mD@Metadata)[...][[mD@Sample.col]] %>%
    keep.samples(mD = mD, samples = ., track = FALSE)
  if (track) {
    mD.return <- add.other.data(mD = mD.return, name = "Sample.filter", x = deparse(substitute(...)))
  }
  return(mD.return)
}
