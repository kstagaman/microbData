#' @title Add or replace data in appropriate slot of a microbData object
#' @description Functions for adding/replacing metadata, abundances, features, phylogeny, distance matrices, and other data to an established microbData object.
#' @aliases replace.metadata
#' @aliases replace.abundances
#' @aliases add.distance.matrices
#' @aliases add.features
#' @aliases add.other.data
#' @aliases add.phylogeny
#' @param mD required; the microbData object to be updated.
#' @param new.tbl required for `replace.metadata` and `replace.abundances`; a data.table or matrix to replace the current Metadata or Abundances slot. Must have the same samples and/or features names as table it replaces.
#' @param distance.matrices required for `add.distance.matrices`; a single distance matrix  of class "dist" or a list of distance matrices of class "dist".
#' @param features required for `add.features`; data.table, data.frame, or matrix table containing taxonomic or functional assignments for each feature (taxon/function). If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Feature" (which will be used to infer feature names). If a keyed data.table (see \code{\link[data.table]{setkey}}), `feature.names` will be set from the keyed column.
#' @param feature.names character; a vector of feature (taxon/function) names. If not provided directly here, will be inferred from `features`. Default is NULL.
#' @param phylo required for `add.phylogeny`; phylo, a phylogenetic tree object.
#' @param x required; a list, vector, or array (data.frame, data.table, matrix, ...) to add to microbData object in the appropriate (usually Other.data) slot.
#' @param name required for `add.other.data`; the name to give the object in the Other.data list (e.g., "Ordinations" or "Beta.metrics").
#' @param replace logical; for `add.other.data`. If TRUE, the new data will replace any data currently in that element of the Other.data list. If FALSE, the new data will be appended to any data currently in that element of the Other.data list. Default is FALSE.
#' @seealso \code{\link[microbData]{microbData}}, \code{\link{slot}}

####################################
#' @name replace.metadata
#' @title Replace Metadata
#' @description Replace the data.table in the Metadata slot of a microbData object with a new one
#' @rdname update.microbData
#' @export

replace.metadata <- function(mD, new.tbl) {
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  if (!{"data.table" %in% class(new.tbl)}) {
    new.tbl <- as.data.table(metadata, keep.rownames = mD@Sample.col) %>%
      setkeyv(mD@Sample.col)
  } else if (!{"sorted" %in% names(attributes(new.tbl))}) {
    rlang::abort(
      "The data.table supplied to `new.tbl' is not keyed by sample names, please supply sample names by using `data.table::setkey'."
    )
  } else {
    new.smpl.col <- attributes(new.tbl)$sorted
  }
  if (!identical(sort(mD@Sample.names), sort(new.tbl[[new.smpl.col]]))) {
    rlang::abort(
      "The data.table supplied to `new.tbl' does not have the same samples as the Metadata table it is replacing."
      )
  } else {
    mD@Metadata <- new.tbl
    if (mD@Sample.col != new.smpl.col) { mD@Sample.col <- new.smpl.col }
    return(mD)
  }
}

####################################
#' @name replace.abundances
#' @title Replace Abundances Table
#' @description Replace the matrix in the Abundances slot of a microbData object with a new one
#' @rdname update.microbData
#' @export

replace.abundances <- function(mD, new.tbl) {
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  if (!{"matrix" %in% class(new.tbl)}) {
    rlang::abort(
      "The table supplied to `new.tbl' must be a matrix."
    )
  } else if (!identical(sort(rownames(mD@Abundances)), sort(rownames(new.tbl)))) {
    rlang::abort(
      "The matrix supplied to `new.tbl' does not have the same samples as the Abundances table it is replacing."
    )
  } else if (!identical(sort(colnames(mD@Abundances)), sort(colnames(new.tbl)))) {
    rlang::abort(
      "The matrix supplied to `new.tbl' does not have the same features as the Abundances table it is replacing."
    )
  } else {
    mD@Abundances <- new.tbl[rownames(mD@Abundances), colnames(mD@Abundances)]
    return(mD)
  }
}

####################################
#' @name add.distance.matrices
#' @title Add Distance Matrices
#' @description Add one or a list of distance matrices to an already created microbData object
#' @rdname update.microbData
#' @export

add.distance.matrices <- function(mD, distance.matrices) {
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
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
  } else {
    mD@Distance.matrices <- distance.matrices
    return(mD)
  }
}

####################################
#' @name add.features
#' @title Add Features Table
#' @description Add a features table to an already created microbData object
#' @rdname update.microbData
#' @export

add.features <- function(mD, features, feature.names = NULL) {
  table.classes <- c("matrix", "data.frame", "data.table")
  if (!any(class(features) %in% table.classes)) {
    rlang::abort(
      "The table supplied to `features' must be of class `data.table', `data.frame', or `matrix'"
    )
  }
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  if (!{"data.table" %in% class(features)}) {
    features <- as.data.table(features, keep.rownames = "Feature")
    setkey(features, Feature)
  }
  if (is.null(feature.names)) {
    if (!{"sorted" %in% names(attributes(features))}) {
      rlang::abort(
        "The data.table supplied to `features' is not sorted by feature names and `feature.names' is also NULL, please supply feature names by either using `data.table::setkey' on the data.table or providing a character vector of feature names."
      )
    } else {
      feat.col <- attributes(features)$sorted
    }
  } else {
    mD@Feature.names <- feature.names
  }

  mD@Features <- features
  mD@Feature.col <- feat.col
  return(mD)
}

####################################
#' @name add.phylogeny
#' @title Add a Phylogenetic Tree to microbData Object
#' @description Add a phylogenetic to an already created microbData object
#' @rdname update.microbData
#' @export

add.phylogeny <- function(mD, phylo) {
  if (!any(class(phylo) %in% "phylo")) {
    rlang::abort(
      "The tree supplied to `phylo' must be of class `phylo'"
    )
  }
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  if (identical(sort(phylo$tip.label), sort(colnames(mD@Abundances)))) {
    mD@Phylogeny <- phylo
    return(mD)
  } else {
    rlang::abort(
      "The tip labels in the provided `phylo' do NOT match the feature names in the `mD' abundance table"
    )
  }
}

####################################
#' @name add.other.data
#' @title Add Other.data to a microbData Object
#' @description Add a vector, list, etc... to an already created microbData object in the Other.data slot
#' @rdname update.microbData
#' @export

add.other.data <- function(mD, name, x, replace = FALSE) {
  if (class(name) != "character") {
    rlang::abort(
      "The argument supplied to `name' must be of class `character'"
    )
  }
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  if (is.null(mD@Other.data[[name]]) | replace) {
    mD@Other.data[[name]] <- x
  } else {
    mD@Other.data[[name]] <- c(mD@Other.data[[name]], x)
  }
  return(mD)
}

