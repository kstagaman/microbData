#' @name glom.assignments
#' @title Agglomerate Feature Abundances by Level
#' @description Sum feature counts by the specified assignment level per sample. E.g., sum ASVs by their Genus assignment.
#' @param mD required; the \code{microbData} object to be glommed. It must contain a Assignments table.
#' @param level required; a character matching a column name in the Assignments table by which to glom (sum) abundances.
#' @param decreasing.assignments logical or NULL; in the Assignments table, are assignments decreasing in hierarchy (e.g. Phylum down to Genus), or increasing (e.g. KO up to Pathway). If you do not want the returning Feature table trimmed of lower level assignments, set to NULL. Default is TRUE.
#' @details This function will sum feature abundances for any assignment level that is present in the Assignments table. Note the \code{decreasing.assignments} argument. It is assumed that the Assignments table will be organized in some sort of hierarchy, and will only keep the columns of that feature table that are the same level as an higher than the column passed to \code{level}.
#' @details Also note that this function will drop any phylogenetic tree that was present as it cannot accurately reflect the phylogeny of the glommed assignments. Also note that things like count transformation and alpha-, beta-diversity, and beta-dispersion results already added to Metadata, or notes about them in Other.data will remain there, even if no longer accurate. Therefore it is important to run this function on a \code{microbData} object with raw counts and no other analyses conducted on it, beyond, perhaps, filtering.
#' @returns a \code{microbData} object with a new Abundances and Feature table that reflect the glomming process.
#' @export

glom.assignments <- function(mD, level, decreasing.assignments = TRUE) {
  if (class(mD) != "microbData") {
    rlang::abort("Argument `mD' must be an object of class `microbData'")
  }
  if (is.null(mD@Assignments)) {
    rlang::abort("There is no Assignments table in this microbData object. Please add a Assignments table and retry.")
  }
  if (!{level %in% names(mD@Assignments)}) {
    rlang::abort(
      paste0(
        "The argument supplied to `level' must match a column name in the Assignments table. For this microbData object, the following are available:\n",
        paste(paste0('  "', names(mD@Assignments), '"'), collapse = "\n")
      )
    )
  }

  assignments <- unique(mD@Assignments[[level]])
  abund.dt <- get.abundances(mD = mD, as.DT = T) %>% setkeyv(mD@Sample.col)
  res.dt <- copy(abund.dt[, .SD, .SDcols = mD@Sample.col])
  new.assignments <- NULL
  for (assignment in assignments) {
    a.assignments <- mD@Assignments[get(level) == assignment][[mD@Feature.col]]
    res.dt <- abund.dt[, .(get(mD@Sample.col), rowSums(.SD)), .SDcols = a.assignments] %>%
      set_names(c(mD@Sample.col, assignment)) %>%
      setkeyv(mD@Sample.col) %>%
      merge(x = res.dt, y = ., all.x = T)
    if (is.null(decreasing.assignments)) {
      keep.col <- names(mD@Assignments)
    } else {
      if (decreasing.assignments) {
        keep.cols <- names(mD@Assignments)[1:str_which(names(mD@Assignments), level)]
      } else {
        keep.cols <- names(mD@Assignments)[str_which(names(mD@Assignments), level):ncol(mD@Assignments)]
      }
    }
    new.assignments <- rbind(
      new.assignments,
      mD@Assignments[get(level) == assignment][1, ..keep.cols]
    )
  }

  res.mat <- as.matrix(res.dt[, -1]) %>%
    set_rownames(res.dt[[mD@Sample.col]]) %>%
    .[, names(sort(colSums(.), decreasing = T))]
  setkeyv(new.assignments, level)

  mD.return <- microbData(
    metadata = mD@Metadata,
    abundances = res.mat,
    assignments = new.assignments
  )
  if (!is.null(mD@Other.data)) {
    mD.return@Other.data <- mD@Other.data
  }
  mD.return <- add.other.data(mD = mD.return, name = "Assignments glommed", x = level)
  return(mD.return)
}
