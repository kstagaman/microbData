#' @title Keep, Drop, or Filter Samples in a microbData object
#' @description Functions to specify specific samples to keep/drop by name or to filter samples based on variable(s) in the Metadata
#' @aliases keep.samples
#' @aliases drop.samples
#' @aliases filter.samples
#' @param mD required; the microbData object to be updated.
#' @param samples character or logical; a vector of sample names to keep or drop. Also can be a logical vector that is either named or in the exact same order as those in the Sample Names slot.
#' @param ... logical expression(s) by which to filter samples based on Metadata, e.g. `Weight > 10 & Genotype != "WT"`.
#' @seealso \code{\link[microbData]{microbData}}


#' @name keep.samples
#' @title Keep Samples
#' @description Create a new microbData objects with just the specified samples
#' @rdname subset.samples
#' @export

keep.samples <- function(mD, samples) {
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
#' @description Create a new microbData objects without the specified samples
#' @rdname subset.samples
#' @export

drop.samples <- function(mD, samples) {
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
#' @description Create a new microbData objects with samples that match the filter
#' @rdname subset.samples
#' @export

filter.samples <- function(mD, ...) {
  copy(mD@Metadata)[...][[mD@Sample.col]] %>%
    keep.samples(mD = mD, samples = .) %>%
    return()
}
