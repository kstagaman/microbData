#' @name microbData
#' @title Create a \code{microbData} object
#' @description Creates an object with associated microbiome data of class \code{microbData}
#' @param metadata required; must be a data.table, data.frame, or matrix with sample metadata. If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Sample" (which will be used to infer sample names). If a keyed data.table (see \code{\link[data.table]{setkey}}), \code{sample.names} will be set from the keyed column.
#' @param abundances required; must be a matrix of feature (taxon/function) abundance values, either counts or relative. This table will be coerced to have sample names as row names since that is more often the required orientation for other tools.
#' @param features data.table, data.frame, or matrix; the table containing taxonomic or functional assignments for each feature (taxon/function). If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Feature" (which will be used to infer feature names). If a keyed data.table (see \code{\link[data.table]{setkey}}), \code{feature.names} will be set from the keyed column. Default is NULL.
#' @param phylogeny phylo; a phylogenetic tree of features (taxonomic). Default is NULL.
#' @param distance.matrices dist or list; a single distance matrix of class "dist" or a list of distance matrices. Default is NULL.
#' @param sample.names character; a vector of the sample names. If not provided directly here, will be inferred from \code{metadata}. Default is NULL.
#' @param feature.names character; a vector of feature (taxon/function) names. If not provided directly here, will be inferred from \code{abundances}. Default is NULL.
#' @param other.data list; a named list of other data to associate with the microbiome data. Can be things like covariate categories of interest or alpha- and beta-diveristy metrics to be included in the analysis. This is for your reference only and will not be implicitly used by any of the associated functions in the \code{microbData} package.
#' @details This function creates a \code{microbData} object that is designed to make it simple for the user to perform basic microbiome analyses such as alpha- and beta-diversity estimation, sample and feature filtering, ordination, etc. The \code{microbData} object stores information like what metrics have been estimated from the data as well as the option to store distance matrices and ordinations in association with the underlying data to make it simple for the user to keep track of analysis steps as well as sharing data and results with colleagues. Furthermore, \code{microbData} objects utilize the power of \code{data.table}s for fast filtering, summarizing and merging.
#' @returns A \code{microbData} object.
#' @slot Metadata A \code{data.table} containing values for the covariates associated with each sample.
#' @slot Abundances A \code{matrix} with abundance counts for each feature, e.g., ASVs, KOs, etc.
#' @slot Features A \code{data.table} containing higher order assignments for each feature, e.g., taxonomy for ASVs or modules and pathways for KOS.
#' @slot Phylogeny A phylogenetic tree for features.
#' @slot Distance.matrices A distance matrix or a list of distance matrices for beta-diversity between each sample.
#' @slot Sample.names A vector of the sample IDs (if not supplied directly, taken from the key column in the Metadata table).
#' @slot Feature.names A vector of the feature IDs (if not supplied directly, taken from the column names in the Abundances table).
#' @slot Sample.col A character string that identifies the column in the Metadata table that contains the sample names (primarily use internally, for consistency).
#' @slot Feature.col A character string that identifies the column in the Features table that contains the feature names (primarily use internally, for consistency).
#' @slot Other.data A list of any other elements you want to associate with this data. Many of the functions in this package add information to this slot for tracking steps and results.
#' @examples
#' ## load data
#' data("metadata_dt")  # loads example metadata.dt
#' setkey(metadata.dt, Sample) # this will tell \code{microbData} that this column contains our sample names
#' data("asv_mat")      # loads example asv.mat
#' data("taxonomy_dt")  # loads example taxonomy.dt
#' setkey(taxonomy.dt, ASV) # this will tell \code{microbData} that this column contains our feature names
#' data("phylogeny")    # loads example phylogeny
#'
#' ## create microbData object
#'
#' mD1 <- microbData(
#'   metadata = metadata.dt,
#'   abundances = asv.mat,
#'   features = taxonomy.dt,
#'   phylogeny = phylogeny
#' )
#' print(mD1)
#' @export

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
      "The data.table supplied to `metadata' is not keyed by sample names and `sample.names' is also NULL, please supply sample names by either using `data.table::setkey' on the data.table or providing a character vector of sample names."
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
    feature.names <- colnames(abundances)[order(colSums(abundances), decreasing = T)]
  }
  return(
    new(
      "microbData",
      Metadata = metadata,
      Abundances = abundances[, order(colSums(abundances), decreasing = T)],
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
