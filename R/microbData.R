#' @name microbData
#' @title Create a `microbData` object
#' @description Creates an object with associated microbiome data of class `microbData`
#' @param metadata required; must be a data.table, data.frame, or matrix with sample metadata. If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Sample" (which will be used to infer sample names). If a keyed data.table (see \code{\link[data.table]{setkey}}), `sample.names` will be set from the keyed column.
#' @param abundances required; must be a matrix of feature (taxon/function) abundance values, either counts or relative. This table will be coerced to have sample names as rownames since that is more often the required orientation for other tools.
#' @param features data.table, data.frame, or matrix; the table containing taxonomic or functional assignments for each feature (taxon/function). If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Feature" (which will be used to infer feature names). If a keyed data.table (see \code{\link[data.table]{setkey}}), `feature.names` will be set from the keyed column. Default is NULL.
#' @param phylogeny phylo; a phylogenetic tree of features (taxonomic). Default is NULL.
#' @param distance.matrices dist or list; a single distance matrix of class "dist" or a list of distance matrices. Default is NULL.
#' @param sample.names character; a vector of the sample names. If not provided directly here, will be inferred from `metadata`. Default is NULL.
#' @param feature.names character; a vector of feature (taxon/function) names. If not provided directly here, will be inferred from `abundances`. Default is NULL.
#' @param other.data list; a named list of other data to associate with the microbiome data. Can be things like covariate categories of interest or alpha- and beta-diveristy metrics to be inlcuded in the analysis. This is for your reference only and will not be implicitly used by any of the associated functions in the microbData package.
#' @seealso \code{\link{head}}
#' @export
#' @examples
#'

microbData <- function(
    metadata,
    abundances,
    features = NULL,
    phylogeny = NULL,
    distance.matrices = NULL,
    sample.names = NULL,
    feature.names = NULL,
    other.data = NULL
) {
  table.classes <- c("matrix", "data.frame", "data.table")
  feat.col <- NULL
  if (!any(class(metadata) %in% table.classes)) {
    rlang::abort(
      "The table supplied to `metadata' must be of class `data.table', `data.frame', or `matrix'"
    )
  }
  if (!{"matrix" %in% class(abundances)}) {
    rlang::abort(
      "The table supplied to `abundances' must be of class `matrix'"
    )
  }
  if (!is.null(features)) {
    if (!any(class(features) %in% table.classes)) {
      rlang::abort(
        "The table supplied to `features' must be of class `data.table', `data.frame', or `matrix'"
      )
    }
  }
  if (!is.null(phylogeny)) {
    if (!{"phylo" %in% class(phylogeny)}) {
      rlang::abort(
        "The tree supplied to `phylogeny' must be of class `phylo'"
      )
    }
  }
  if (!is.null(distance.matrices)) {
    if (class(distance.matrices) == "list") {
      for (i in seq_along(distance.matrices)) {
        if (class(distance.matrices[[i]]) != "dist") {
          rlang::abort(
            paste(
            "The elements in the list supplied to `distance.matrices' must all be of class `dist'.",
            "Element at index", i, "is of class", class(distance.matrices[[i]])
            )
          )
        }
      }
    } else if (class(distance.matrices) != "dist") {
      rlang::abort(
        "The object supplied to `distance.matrices' must be of class `dist' or `list'"
      )
    }
  }
  if (!is.null(other.data)) {
    if (class(other.data) != "list") {
      rlang::abort(
        "The data supplied to `other.data' must be of class `list'"
      )
    }
    if (is.null(names(other.data))) {
      rlang::abort(
        "The list supplied to `other.data' must have names"
      )
    }
  }

  if (!{"data.table" %in% class(metadata)}) {
    metadata <- as.data.table(metadata, keep.rownames = "Sample")
    setkey(metadata, Sample)
    smpl.col <- "Sample"
  }
  if (!is.null(features)) {
    if (!{"data.table" %in% class(features)}) {
      features <- as.data.table(features, keep.rownames = "Feature")
      setkey(features, Feature)
    }
  }

  if (!{"sorted" %in% names(attributes(metadata))} & is.null(sample.names)) {
    rlang::abort(
      "The data.table supplied to `metadata' is not sorted by sample names and `sample.names' is also NULL, please supply sample names by either using `data.table::setkey' on the data.table or providing a character vector of sample names."
    )
  } else {
    smpl.col <- attributes(metadata)$sorted
  }
  if (!is.null(features)) {
    if (!{"sorted" %in% names(attributes(features))} & is.null(feature.names)) {
      rlang::abort(
        "The data.table supplied to `features' is not sorted by feature names and `feature.names' is also NULL, please supply feature names by either using `data.table::setkey' on the data.table or providing a character vector of feature names."
      )
    } else {
      feat.col <- attributes(features)$sorted
    }
  }

  if (is.null(sample.names)) {
    sample.names <- as.character(metadata[[attributes(metadata)$sorted]])
  }
  if (!identical(sort(sample.names), sort(rownames(abundances)))) {
    if (identical(sort(sample.names), sort(colnames(abundances)))) {
      abundances <- t(abundances)
    } else {
      rlang::abort(
        "The sample names supplied in either `sample.names' or the key column of `metadata' do no match the sample names in the `abundances' matrix"
      )
    }
  }
  if (is.null(feature.names)) {
    feature.names <- sort(colnames(abundances))
  }
  return(
    new(
      "microbData",
      Metadata = metadata,
      Abundances = abundances,
      Features = features,
      Phylogeny = phylogeny,
      Distance.matrices = distance.matrices,
      Sample.names = sample.names,
      Feature.names = feature.names,
      Sample.col = smpl.col,
      Feature.col = feat.col,
      Other.data = other.data
    )
  )
}
