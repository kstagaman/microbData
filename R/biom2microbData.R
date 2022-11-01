#' @name biom2microbData
#' @title BIOM to microbData
#' @description Convert a biom object to a \code{microbData} object
#' @param biom required; the biom object to be converted to a \code{microbData} object.
#' @param feature.prefix character; If NULL, feature IDs will remain unchanged. If a character, will be the prefix followed by integers to re-ID features. E.g., if \code{feature.prefix = "ASV"}, features will be renamed "ASV001", "ASV002", ...with \code{\link[numbered.features]{numbered.features}}. The number of zeroes in the IDs will be determined by the total number of features. Default is NULL.
#' @param rename.NA.assignments logical; if FALSE, NAs in the Feature Assignments table will be ignored. If TRUE, NAs in the Feature Assignments table will be renamed with \code{\link[microbData]{rename.NA.taxa}}. Default is FALSE
#' @param ... additional arguments to pass to \code{\link[microbData]{numbered.features}} and/or \code{\link[microbData]{rename.NA.taxa}}.
#' @seealso \code{\link[biom]{biom}}, \code{\link[microbData]{numbered.features}}, \code{\link[microbData]{rename.NA.taxa}}
#' @export

biom2microbData <- function(
    biom,
    feature.prefix = NULL,
    rename.NA.assignments = FALSE,
    ...
) {
  require(biomformat)
  vargs <- rlang::list2(...)
  if (!{"biom" %in% class(biom)}) {
    rlang::abort("Object supplied to argument `biom' must be of class 'biom'.")
  }
  smpl.dt <- as(biomformat::sample_metadata(biom), "data.frame") %>%
    as.data.table(keep.rownames = "Sample") %>%
    setkey(Sample)
  abund.mat <- as(biomformat::biom_data(biom), "matrix") # Sort abundance table if it's not sequences
  if (!identical(row.names(as(sample_metadata(biom), "data.frame")), rownames(abund.mat))) {
    abund.mat <- t(abund.mat)
  }
  abund.mat <- abund.mat[, names(sort(colSums(abund.mat), decreasing = T))]
  mD <- microbData(metadata = smpl.dt, abundances = abund.mat)
  if (!is.null(biomformat::observation_metadata(biom))) {
    assign.dt <- as(observation_metadata(biom), "matrix") %>%
      as.data.table(keep.rownames = "Feature") %>%
      setkey(Feature)
    mD <- add.assignments(assignments = assign.dt, mD = mD)
  }
  if (!is.null(feature.prefix)) {
    mD <- numbered.features(
      mD,
      prefix = feature.prefix,
      old.IDs.file = vargs$old.IDs.file
    )
  }
  if (rename.NA.assignments) {
    if (is.null(vargs$force.split)) { vargs$force.split <- FALSE }
    mD <- rename.NA.assignments(
      mD,
      force.split = vargs$force.split,
      level.order = vargs$level.order
    )
  }
  return(mD)
}
