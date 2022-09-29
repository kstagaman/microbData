#' @name get.microbData
#' @title Get data from appropriate slot of a microbData object
#' @description Functions grabbing features, phylogeny, distance matrices, or other data from a microbData object. The function `get.microbData` allows the user to programmatically get data from various slots as it takes the slot name as a character and returns the data in that slot. The other functions are very simple wrappers for the `@` extraction operator to make code more human-readable, or just because they may be easier to understand.
#' @aliases get.abundances
#' @aliases get.distance.matrices
#' @aliases get.features
#' @aliases get.feature.col
#' @aliases get.feature.names
#' @aliases get.metadata
#' @aliases get.other.data
#' @aliases get.phylogeny
#' @aliases get.sample.col
#' @aliases get.sample.names
#' @aliases get.covariate.names
#' @param mD required; the microbData object to get data from.
#' @param slot.name required for `get.microbData`; the name of the slot (e.g., Metadata, Abundances, ...) you wish to extract from the `mD`. To see what slot names are valid, use the function \code{\link{slotNames}} on the mD.
#' @param as.DT logical; for `get.microbData` and `get.abundances`. The Abundances table is the only one stored as a matrix in the microbData. However, if you would like to retrieve it and convert it to a data.table, you can set this to TRUE. Default is FALSE.
#' @param location required for `get.other.data`; the index or name of the list element to get from the Other.data slot.
#' @seealso \code{\link{slot}}, \code{\link{\@}}

####################################
#' @title Get Abundances
#' @description Get the abundance table from a microbData object.
#' @rdname get.microbData
#' @export

get.microbData <- function(mD, slot.name, as.DT = FALSE) {
  tbl <- slot(mD, slot.name)
  if (as.DT & !{c("data.table", "phylo") %in% class(tbl)} ) {
    return(as.data.table(tbl, keep.rownames = mD@Sample.col))
  } else {
    return(tbl)
  }
}

####################################
#' @name get.abundances
#' @title Get Abundances
#' @description Get the abundance table from a microbData object.
#' @rdname get.microbData
#' @export

get.abundances <- function(mD, as.DT = FALSE) {
  if (as.DT) {
    return(as.data.table(mD@Abundances, keep.rownames = mD@Sample.col))
  } else {
    return(mD@Abundances)
  }
}

####################################
#' @name get.distance.matrices
#' @title Get Distance Matrices
#' @description Get the distance matrices from a microbData object.
#' @rdname get.microbData
#' @export

get.distance.matrices <- function(mD) {
  return(mD@Distance.matrices)
}

####################################
#' @name get.features
#' @title Get Features
#' @description Get the features table from a microbData object.
#' @rdname get.microbData
#' @export

get.features <- function(mD) {
  return(mD@Features)
}

####################################
#' @name get.feature.col
#' @title Get Feature Column Name
#' @description Get the name of the column containing specific feature IDs from the features table in a microbData object.
#' @rdname get.microbData
#' @export

get.feature.col <- function(mD) {
  return(mD@Feature.col)
}

####################################
#' @name get.feature.names
#' @title Get Feature Names
#' @description Get the specific IDs of features in a microbData object.
#' @rdname get.microbData
#' @export

get.feature.names <- function(mD) {
  return(mD@Feature.names)
}

####################################
#' @name get.metadata
#' @title Get Metadata
#' @description Get the metadata table from a microbData object.
#' @rdname get.microbData
#' @export

get.metadata <- function(mD) {
  return(mD@Metadata)
}

####################################
#' @name get.other.data
#' @title Get Something from Other.data
#' @description Get an element from the Other.data list of a microbData object.
#' @rdname get.microbData
#' @export

get.other.data <- function(mD, location) {
  return(mD@Other.data[[location]])
}

####################################
#' @name get.phylogeny
#' @title Get Phylogeny
#' @description Get the phylogenetic tree from a microbData object.
#' @rdname get.microbData
#' @export

get.phylogeny <- function(mD) {
  return(mD@Phylogeny)
}

####################################
#' @name get.sample.col
#' @title Get Sample Column Name
#' @description Get the name of the column containing specific sample IDs from the samples table in a microbData object.
#' @rdname get.microbData
#' @export

get.sample.col <- function(mD) {
  return(mD@Sample.col)
}

####################################
#' @name get.sample.names
#' @title Get Sample Names
#' @description Get the specific IDs of samples in a microbData object.
#' @rdname get.microbData
#' @export

get.sample.names <- function(mD) {
  return(mD@Sample.names)
}

####################################
#' @name get.covariate.names
#' @title Get Covariate Names
#' @description Get the column names from the Metadata (minus sample column name) in a microbData object.
#' @rdname get.microbData
#' @export

get.covariate.names <- function(mD) {
  names(mD@Metadata) %>%
    str_subset(pattern = paste0("^", mD@Sample.col, "$"), negate = T) %>%
    return()
}

