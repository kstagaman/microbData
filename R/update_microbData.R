#' @title Add or replace data in appropriate slot of a \code{microbData} object
#' @description Functions for adding/replacing metadata, abundances, features, phylogeny, distance matrices, and other data to an established \code{microbData} object.
#' @aliases replace.metadata
#' @aliases replace.abundances
#' @aliases add.distance.matrices
#' @aliases add.assignments
#' @aliases replace.assignments
#' @aliases add.other.data
#' @aliases add.phylogeny
#' @aliases replace.phylogeny
#' @param mD required; the \code{microbData} object to be updated.
#' @param new.tbl required for \code{replace.metadata} and \code{replace.abundances}; a data.table or matrix to replace the current Metadata or Abundances slot. Must have the same samples and/or features names as table it replaces.
#' @param distance.matrices required for \code{add.distance.matrices}; a single distance matrix  of class "dist" or a list of distance matrices of class "dist".
#' @param features required for \code{add.assignments}; data.table, data.frame, or matrix table containing taxonomic or functional assignments for each feature (taxon/function). If not already a data.table, the data will be converted to a data.table and the row names will be saved under a column called "Feature" (which will be used to infer feature names). If a keyed data.table (see \code{\link[data.table]{setkey}}), \code{feature.names} will be set from the keyed column.
#' @param feature.names character; a vector of feature (taxon/function) names. If not provided directly here, will be inferred from \code{features}. Default is NULL.
#' @param phylo required for \code{add.phylogeny}; phylo, a phylogenetic tree object.
#' @param x required; a list, vector, or array (data.frame, data.table, matrix, ...) to add to \code{microbData} object in the appropriate (usually Other.data) slot.
#' @param name required for \code{add.other.data}; the name to give the object in the Other.data list (e.g., "Ordinations" or "Beta.metrics").
#' @param replace logical; for \code{add.other.data}. If TRUE, the new data will replace any data currently in that element of the Other.data list. If FALSE, the new data will be appended to any data currently in that element of the Other.data list. Default is FALSE.
#' @details These function can replace data in or add data to an existing \code{microbData} object. The \code{@} operator, or the \code{\link{slot}} function can also accomplish this same task. However, for replacing metadata, and abundances, these functions check to make sure the replacement tables match the samples and/or features of the original tables, helping to safeguard against user error. Additionally, these functions are more user friendly and make sure data get into the correct slot than using the base functions.
#' @returns A \code{microbData} object with the new data.
#' @seealso \code{\link[microbData]{microbData}}, \code{\link{slot}}

####################################
#' @name replace.metadata
#' @title Replace Metadata
#' @description Replace the data.table in the Metadata slot of a \code{microbData} object with a new one
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
#' @description Replace the matrix in the Abundances slot of a \code{microbData} object with a new one
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
  } else {
    mD@Abundances <- new.tbl[
      rownames(mD@Abundances),
      names(sort(colSums(new.tbl), decreasing = T))
    ]
    if (!is.null(mD@Assignments)) {
      if (!identical(sort(colnames(mD@Abundances)), sort(mD@Feature[[mD@Feature.col]]))) {
        rlang::warn(
          "The feature names (colnames) in the new Abundances table are not all identical to the names in the Assignments table."
        )
      }
    }
    if (!is.null(mD@Phylogeny)) {
      if (!identical(sort(colnames(mD@Abundances)), sort(mD@Phylogeny$tip.label))) {
        rlang::warn(
          "The feature names (colnames) in the new Abundances table are not all identical to the names in the Phylogenetic tree."
        )
      }
    }
    if (!identical(sort(colnames(mD@Abundances)), sort(mD@Feature.names))) {
      mD@Feature.names <- colnames(mD@Abundances)
    }
    return(mD)
  }
}

####################################
#' @name add.distance.matrices
#' @title Add Distance Matrices
#' @description Add one or a list of distance matrices to an already created \code{microbData} object
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
#' @name add.assignments
#' @title Add Assignments Table
#' @description Add a features table, or replace existing features table in an already created \code{microbData} object. An alias for this function is \code{replace.assignments}.
#' @rdname update.microbData
#' @export

add.assignments <- function(mD, features, feature.names = NULL) {
  table.classes <- c("matrix", "data.frame", "data.table")
  if (!any(class(features) %in% table.classes)) {
    rlang::abort(
      "The table supplied to `features' must be of class `data.table', `data.frame', or `matrix'"
    )
  }
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  feat.col <- ifelse(is.null(mD@Assignments), "Feature", mD@Feature.col)
  if (!{"data.table" %in% class(features)}) {
    features <- as.data.table(features, keep.rownames = feat.col) %>%
      setkeyv(feat.col)
  } else if (is.null(feature.names)) {
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
  if (!is.null(mD@Assignments)) {
    if (!identical(sort(mD@Feature.names), sort(features[[feat.col]]))) {
      rlang::abort(
        "The data.table supplied to `features' does not have the same feature names as the Assignments table it is replacing."
      )
    }
  }
  mD@Assignments <- features
  if (is.null(mD@Feature.col) | mD@Feature.col != feat.col) { mD@Feature.col <- feat.col }
  return(mD)
}

#' @export

replace.assignments <- add.assignments

####################################
#' @name add.phylogeny
#' @title Add a Phylogenetic Tree to \code{microbData} Object
#' @description Add a phylogenetic to an already created \code{microbData} object
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

#' @export

replace.phylogeny <- add.phylogeny

####################################
#' @name add.other.data
#' @title Add Other.data to a \code{microbData} Object
#' @description Add a vector, list, etc... to an already created \code{microbData} object in the Other.data slot
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

