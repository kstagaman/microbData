#' @name merge.microbData
#' @title Merge Two MicrobDatas
#' @description Merge two microbData objects.
#' @param mD1 required; a microbData object you want to merge.
#' @param mD2 required; a microbData object you want to merge.
#' @seealso \code{\link[microbData]{microbData}}
#' @export

merge.microbData <- function(mD1, mD2) {
  if (mD1@Sample.col != mD2@Sample.col) {
    rlang::abort(
      "The Sample Column names for these two microbData object do not match, and therefore cannot be merged"
    )
  }
  merged.abunds.dt <- rbindlist(
    list(
      as.data.table(mD1@Abundances, keep.rownames = mD1@Sample.col),
      as.data.table(mD2@Abundances, keep.rownames = mD1@Sample.col)
    ),
    fill = TRUE
  )
  sample.names <- merged.abunds.dt[[mD1@Sample.col]]
  merged.abunds.mat <- as.matrix(merged.abunds.dt[, (mD1@Sample.col) := NULL])
  rownames(merged.abunds.mat) <- sample.names
  merged.abunds.mat[is.na(merged.abunds.mat)] <- 0
  new.mD <- microbData(
    metadata = rbindlist(list(mD1@Metadata, mD2@Metadata), fill = TRUE),
    abundances = merged.abunds.mat
  )
  if (!is.null(mD1@Features) & !is.null(mD2@Features)) {
    merged.features <- merge(mD1@Features, mD2@Features, by = mD1@Feature.col, all = TRUE)
    new.mD <- add.features(features = merged.features, mD = new.mD)
  } else if (!is.null(mD1@Features) | !is.null(mD2@Features)) {
    rlang::inform(
      "One of the supplied microbData objects had a Features table, but the other did not, so none has been included in the merged microbData object."
    )
  }
  if (!is.null(mD1@Phylogeny) & !is.null(mD2@Phylogeny)) {
    if (identical(mD1@Phylogeny, mD2@Phylogeny)) {
      new.mD <- add.phylogeny(mD1@Phylogeny)
    } else {
      rlang::inform(
        "The phylogenies of the two microbData objects are not identical, so none has been included in the merged microbData object. Please, re-infer a phylogeny for the merged samples together to add to the microbData object."
      )
    }
  } else {
    rlang::inform(
      "One of the supplied microbData objects had a Phylogeny, but the other did not, so none has been included in the merged microbData object."
    )
  }
  return(new.mD)
}
