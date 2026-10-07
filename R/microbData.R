#' @name microbData
#' @title Create a \code{microbData} object
#' @description Creates an object with associated microbiome data of class \code{microbData}
#' @param metadata required; must be a data.table, data.frame, or matrix with sample metadata. If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Sample" (which will be used to infer sample names). If a keyed data.table (see \code{\link[data.table]{setkey}}), \code{sample.names} will be set from the keyed column.
#' @param abundances required; must be a matrix of feature (taxon/function) abundance values, either counts or relative. This table will be coerced to have sample names as row names since that is more often the required orientation for other tools.
#' @param assignments data.table, data.frame, or matrix; the table containing taxonomic or functional assignments for each feature (taxon/function). If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Feature" (which will be used to infer feature names). If a keyed data.table (see \code{\link[data.table]{setkey}}), \code{feature.names} will be set from the keyed column. Default is NULL.
#' @param phylogeny phylo; a phylogenetic tree of features (taxonomic). Default is NULL.
#' @param distance.matrices dist or list; a single distance matrix of class "dist" or a list of distance matrices. Default is NULL.
#' @param sample.names character; a vector of the sample names. If not provided directly here, will be inferred from \code{metadata}. Default is NULL.
#' @param feature.names character; a vector of feature (taxon/function) names. If not provided directly here, will be inferred from \code{abundances}. Default is NULL.
#' @param other.data list; a named list of other data to associate with the microbiome data. Can be things like covariate categories of interest or alpha- and beta-diveristy metrics to be included in the analysis. This is for your reference only and will not be implicitly used by any of the associated functions in the \code{microbData} package.
#' @details This function creates a \code{microbData} object that is designed to make it simple for the user to perform basic microbiome analyses such as alpha- and beta-diversity estimation, sample and feature filtering, ordination, etc. The \code{microbData} object stores information like what metrics have been estimated from the data as well as the option to store distance matrices and ordinations in association with the underlying data to make it simple for the user to keep track of analysis steps as well as sharing data and results with colleagues. Furthermore, \code{microbData} objects utilize the power of \code{data.table}s for fast filtering, summarizing and merging.
#' @returns A \code{microbData} object.
#' @slot Metadata A \code{data.table} containing values for the covariates associated with each sample.
#' @slot Abundances A \code{matrix} with abundance counts for each feature, e.g., ASVs, KOs, etc.
#' @slot Assignments A \code{data.table} containing higher order assignments for each feature, e.g., taxonomy for ASVs or modules and pathways for KOs.
#' @slot Phylogeny A phylogenetic tree for features.
#' @slot Distance.matrices A distance matrix or a list of distance matrices for beta-diversity between each sample.
#' @slot Sample.names A vector of the sample IDs (if not supplied directly, taken from the key column in the Metadata table).
#' @slot Feature.names A vector of the feature IDs (if not supplied directly, taken from the column names in the Abundances table).
#' @slot Sample.col A character string that identifies the column in the Metadata table that contains the sample names (primarily use internally, for consistency).
#' @slot Feature.col A character string that identifies the column in the Assignments table that contains the feature names (primarily use internally, for consistency).
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
#'   assignments = taxonomy.dt,
#'   phylogeny = phylogeny
#' )
#' print(mD1)
#' @export

microbData <- function(
    metadata,
    abundances,
    sample.names.col = NULL,
    sample.names = NULL,
    assignments = NULL,
    feature.names.col = NULL,
    feature.names = NULL,
    phylogeny = NULL,
    distance.matrices = NULL,
    other.data = NULL
) {
  table.classes <- c("matrix", "data.frame", "data.table", "tbl_df")
  feat.col <- NULL
  if (!any(class(metadata) %in% table.classes)) {
    rlang::abort(
      "The table supplied to `metadata' must be of class `data.table`, `data.frame`, `tbl_df`, or `matrix`"
    )
  }
  if (!{"matrix" %in% class(abundances)}) {
    rlang::abort(
      "The table supplied to `abundances` must be of class `matrix`"
    )
  }
  if (!is.null(assignments)) {
    if (!any(class(assignments) %in% table.classes)) {
      rlang::abort(
        "The table supplied to `assignments' must be of class `data.table`, `data.frame`, `tbl_df`, or `matrix`"
      )
    }
  }
  if (!is.null(phylogeny)) {
    if (!{"phylo" %in% class(phylogeny)}) {
      rlang::abort(
        "The tree supplied to `phylogeny` must be of class `phylo`"
      )
    }
  }
  if (!is.null(distance.matrices)) {
    if (class(distance.matrices) == "list") {
      for (i in seq_along(distance.matrices)) {
        if (class(distance.matrices[[i]]) != "dist") {
          rlang::abort(
            paste(
              "The elements in the list supplied to `distance.matrices` must all be of class `dist`.",
              "Element at index", i, "is of class", class(distance.matrices[[i]])
            )
          )
        }
      }
    } else if (class(distance.matrices) != "dist") {
      rlang::abort(
        "The object supplied to `distance.matrices` must be of class `dist` or `list`"
      )
    }
  }
  if (!is.null(other.data)) {
    if (class(other.data) != "list") {
      rlang::abort(
        "The data supplied to `other.data` must be of class `list`"
      )
    }
    if (is.null(names(other.data))) {
      rlang::abort(
        "The list supplied to `other.data` must have names"
      )
    }
  }
  metadata.checker <- table.checkers[[class(metadata)[1]]]
  metadata.check <- metadata.checker(metadata, sample.names.col, sample.names, slot = "metadata")

  if (!is.null(assignments)) {
    assignments.checker <- table.checkers[[class(assignments)[1]]]
    assignments.check <- assignments.checker(metadata, sample.names.col, sample.names, slot = "assignments")
  } else {
    assignments.check <- list(TBL = NULL, NC = NULL, NS = NULL)
  }

  if (!identical(sort(metadata.check$ns), sort(rownames(abundances)))) {
    if (identical(sort(metadata.check$ns), sort(colnames(abundances)))) {
      abundances <- t(abundances)
    } else {
      rlang::abort(
        "The sample names supplied in either `sample.names` or the key column of `metadata` do no match the sample names in the `abundances' matrix"
      )
    }
  }

  return(
    new(
      "microbData",
      Metadata = metadata.check$TBL,
      Abundances = abundances[, order(colSums(abundances), decreasing = T)],
      Assignments = assignments.check$TBL,
      Phylogeny = phylogeny,
      Distance.matrices = distance.matrices,
      Sample.names = metadata.check$SN,
      Feature.names = assignments.check$SN,
      Sample.col = metadata.check$NC,
      Feature.col = assignments.check$NC,
      Other.data = other.data
    )
  )
}


####################################
#' @title Table checkers
#' @description A dictionary of functions to check that certain table are in the correct mD format
#' @noRd

table.checkers <- list(
  data.table = function(tbl, nc = NULL, ns = NULL, slot = c("metadata", "assignments")) {
    if ("sorted" %in% names(attributes(tbl))) {
      smpl.col <- attributes(tbl)$sorted
    } else if (!is.null(nc)) {
      smpl.col <- nc
      setkeyv(tbl, nc)
    } else if (!is.null(ns)) {
      smpl.col <- names(tbl)[vapply(tbl, function(col) setequal(col, ns), logical(1))]
      setkeyv(tbl, smpl.col)
    } else {
      type <- ifelse(slot == "metadata", "sample", "feature")
      rlang::abort(
        sprintf(
          "The data.table supplied to `%1$s` is not keyed by %2$s names and both `%2$s.names.col` & `%2$s.names` are NULL, please use `data.table::setkey` on the data.table or provide the name of the appropriate column or a character vector of names.",
          slot,
          type
        )
      )
    }
    smpl.names <- as.character(tbl[[smpl.col]])
    return(list(TBL = tbl, NC = smpl.col, NS = smpl.names))
  },
  tbl_df = function(tbl, nc = NULL, ns = NULL, slot = c("metadata", "assignments")) {
    if (!is.null(nc)) {
      smpl.col <- nc
    } else if (!is.null(ns)) {
      smpl.col <- names(tbl)[vapply(tbl, function(col) setequal(col, ns), logical(1))]
    } else {
      type <- ifelse(slot == "metadata", "sample", "feature")
      rlang::abort(
        sprintf(
          "Both `%1$s.names.col` & `%1$s.names` are NULL, please supply sample names by either providing the name of the appropriate column or a character vector of names.",
          type
        )
      )
    }
    smpl.names <- as.character(tbl[[smpl.col]])
    return(list(TBL = tbl, NC = smpl.col, NS = smpl.names))
  },
  data.frame = function(tbl, nc = NULL, ns = NULL) {
    if (!is.null(nc)) {
      if (nc == 0 | nc == "rn") {
        tbl %<>% as.data.table(keep.rownames = "Sample")
        setkey(tbl, Sample)
        smpl.col <- "Sample"
      } else {
        tbl %<>% as.data.table()
        setkeyv(tbl, nc)
        smpl.col <- nc
      }
    } else if (!is.null(ns)) {
      if (setequal(row.names(tbl), ns)) {
        smpl.col <- "Sample"
        tbl %<>% as.data.table(keep.rownames = smpl.col)
      } else {
        smpl.col <- names(tbl)[vapply(tbl, function(col) setequal(col, ns), logical(1))]
        tbl %<>% as.data.table()
      }
    } else {
      type <- ifelse(slot == "metadata", "sample", "feature")
      rlang::abort(
        sprintf(
          "Both `%1$s.names.col` & `%1$s.names` are NULL, please supply sample names by either providing the name of the appropriate column or a character vector of names.",
          type
        )
      )
    }
    setkeyv(tbl, smpl.col)
    smpl.names <- as.character(tbl[[smpl.col]])
    return(list(TBL = tbl, NC = smpl.col, NS = smpl.names))
  },
  matrix = function(tbl, nc = NULL, ns = NULL) {
    if (!is.null(nc)) {
      if (nc == 0 | nc == "rn") {
        smpl.col <- "Sample"
        tbl %<>% as.data.table(keep.rownames = smpl.col)
        setkeyv(tbl, smpl.col)
      } else {
        tbl %<>% as.data.table()
        setkeyv(tbl, nc)
        smpl.col <- nc
      }
    } else if (!is.null(ns)) {
      if (setequal(rownames(tbl), ns)) {
        smpl.col <- "Sample"
        tbl %<>% as.data.table(keep.rownames = smpl.col)
      } else {
        smpl.col <- colnames(tbl)[vapply(tbl, function(col) setequal(col, ns), logical(1))]
        tbl %<>% as.data.table()
      }
    } else {
      type <- ifelse(slot == "metadata", "sample", "feature")
      rlang::abort(
        sprintf(
          "Both `%1$s.names.col` & `%1$s.names` are NULL, please supply sample names by either providing the name of the appropriate column or a character vector of names.",
          type
        )
      )
    }
    setkeyv(tbl, smpl.col)
    smpl.names <- as.character(tbl[[smpl.col]])
    return(list(TBL = tbl, NC = smpl.col, NS = smpl.names))
  }
)
