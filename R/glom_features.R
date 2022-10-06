#' @name glom.features
#' @title Agglomerate Feature Abundances by Level
#' @description Sum feature counts by the specified assignment level per sample. E.g., sum ASVs by their Genus assignment.
#' @param mD required; the \code{microbData} object to be glommed. It must contain a Features table.
#' @param level required; a character matching a column name in the Features table by which to glom (sum) abundances.
#' @param decreasing.assignments logical or NULL; in the Features table, are assignments decreasing in hierarchy (e.g. Phylum down to Genus), or increasing (e.g. KO up to Pathway). If you do not want the returning Feature table trimmed of lower level assignments, set to NULL. Default is TRUE.
#' @details This function will sum feature abundances for any assignment level that is present in the Features table. Note the \code{decreasing.assignments} argument. It is assumed that the Features table will be organized in some sort of hierarchy, and will only keep the columns of that feature table that are the same level as an higher than the column passed to \code{level}.
#' @details Also note that this function will drop any phylogenetic tree that was present as it cannot accurately reflect the phylogeny of the glommed features. Also note that things like count transformation and alpha-, beta-diversity, and beta-dispersion results already added to Metadata, or notes about them in Other.data will remain there, even if no longer accurate. Therefore it is important to run this function on a \code{microbData} object with raw counts and no other analyses conducted on it, beyond, perhaps, filtering.
#' @returns a \code{microbData} object with a new Abundances and Feature table that reflect the glomming process.
#' @export

glom.features <- function(mD, level, decreasing.assignments = TRUE) {
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  if (is.null(mD@Features)) {
    rlang::abort("There is no Features table in this microbData object. Please add a Features table and retry.")
  }
  if (!{level %in% names(mD@Features)}) {
    rlang::abort(
      paste0(
        "The argument supplied to `level' must match a column name in the Features table. For this microbData object, the following are available:\n",
        paste(paste0('  "', names(mD@Features), '"'), collapse = "\n")
      )
    )
  }

  assignments <- unique(mD@Features[[level]])
  abund.dt <- get.abundances(mD = mD, as.DT = T) %>% setkeyv(mD@Sample.col)
  res.dt <- copy(abund.dt[, .SD, .SDcols = mD@Sample.col])
  new.features <- NULL
  for (assignment in assignments) {
    a.features <- mD@Features[get(level) == assignment][[mD@Feature.col]]
    res.dt <- abund.dt[, .(get(mD@Sample.col), rowSums(.SD)), .SDcols = a.features] %>%
      set_names(c(mD@Sample.col, assignment)) %>%
      setkeyv(mD@Sample.col) %>%
      merge(x = res.dt, y = ., all.x = T)
    if (is.null(decreasing.assignments)) {
      keep.col <- names(mD@Features)
    } else {
      if (decreasing.assignments) {
        keep.cols <- names(mD@Features)[1:str_which(names(mD@Features), level)]
      } else {
        keep.cols <- names(mD@Features)[str_which(names(mD@Features), level):ncol(mD@Features)]
      }
    }
    new.features <- rbind(
      new.features,
      mD@Features[get(level) == assignment][1, ..keep.cols]
    )
  }

  res.mat <- as.matrix(res.dt[, -1]) %>%
    set_rownames(res.dt[[mD@Sample.col]]) %>%
    .[, names(sort(colSums(.), decreasing = T))]
  setkeyv(new.features, level)

  mD.return <- microbData(
    metadata = mD@Metadata,
    abundances = res.mat,
    features = new.features
  )
  if (!is.null(mD@Other.data)) {
    mD.return@Other.data <- mD@Other.data
  }
  mD.return <- add.other.data(mD = mD.return, name = "Features glommed", x = level)
  return(mD.return)
}
