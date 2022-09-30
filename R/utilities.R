#' @title Utility Functions
#' @description Functions to do quick, helpful things
#' @aliases nsamples
#' @aliases nfeatures
#' @aliases sample.sums
#' @aliases feature.sums
#' @param mD required; the microbData object to be updated.
#' @param sort character; either "increasing", "decreasing", or NULL. Increasing will set `sort(decreasing = FALSE)` and decreasing will set `sort(decreasing = TRUE)`. NULL will not sort output. Default is NULL.
#' @seealso \code{\link[microbData]{microbData}}

#' @name nsamples
#' @title Number of Samples
#' @description Get the number of samples in a microbData object.
#' @rdname utilities
#' @export

nsamples <- function(mD) {
  length(mD@Sample.names) %>%
    return()
}

#' @name nfeatures
#' @title Number of Features
#' @description Get the number of features in a microbData object.
#' @rdname utilities
#' @export

nfeatures <- function(mD) {
  length(mD@Feature.names) %>%
    return()
}

#' @name sample.sums
#' @title Total Count of Features in Each Sample
#' @description Get the total abundance of all features for each sample in a microbData object.
#' @rdname utilities
#' @export

sample.sums <- function(mD, sort = NULL) {
  if (is.null(sort)) {
    rowSums(mD@Abundances) %>%
      return()
  } else {
    sort = tolower(sort)
    rowSums(mD@Abundances) %>%
      sort(
        decreasing = ifelse(
          sort == "decreasing", T,
          ifelse(
            sort == "increasing", F,
            rlang::abort("Argument `sort' must be 'increasing', 'decreasing', or NULL")
          )
        )
      ) %>%
      return()
  }
}

#' @name feature.sums
#' @title Total Count of Features Across All Samples
#' @description Get the total abundance of each features across all samples in a microbData object.
#' @rdname utilities
#' @export

feature.sums <- function(mD, sort = NULL) {
  if (is.null(sort)) {
    colSums(mD@Abundances) %>%
      return()
  } else {
    sort = tolower(sort)
    colSums(mD@Abundances) %>%
      sort(
        decreasing = ifelse(
          sort == "decreasing", T,
          ifelse(
            sort == "increasing", F,
            rlang::abort("Argument `sort' must be 'increasing', 'decreasing', or NULL")
          )
        )
      ) %>%
      return()
  }
}
