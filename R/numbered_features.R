#' @name numbered.features
#' @title Replace feature IDs with an arbitrary prefix and numbers
#' @description This function replaces the feature names in a \code{microbData} object with IDs including a prefix and then numerals, e.g. "ASV0001", "ASV0002", ..., "ASV1234".
#' @param mD required; \code{microbData} object from which data will be read.
#' @param prefix required; character string to append to numeric IDs, e.g., "ASV", "KO".
#' @param old.IDs.file character; name of file to write old IDs (in FASTA, CSV, or RDS format inferred from extention: must be ".fa", ".fasta", ".csv", or ".rds") to. NULL will not save IDs to a file. Default is NULL.
#' @seealso \code{\link[microbData]{update_feature_names}}
#' @export

numbered.features <- function(mD, prefix, old.IDs.file = NULL) {
  if (prefix == "" | length(prefix) == 0 | is.logical(prefix)) {
    rlang::abort(
      "Argument `prefix' must be a non-numeric-only string of length greater than 0 and not logical (TRUE/FALSE)"
    )
  }
  old.ids <- copy(mD@Feature.names)
  n.digits <- nchar(length(old.ids))
  new.ids <- sapply(1:length(old.ids), function(d) {
    zeroes <- paste(rep("0", n.digits - nchar(d)), collapse = "")
    return(paste0(zeroes, d))
  }) %>%
    unlist() %>%
    paste0(prefix, .)
  ids.dt <- data.table(Old = old.ids, New = new.ids)
  setkey(ids.dt, Old)
  if (!is.null(old.IDs.file)) {
    if (str_detect(old.IDs.file, "\\.rds$")) {
      names(old.ids) <- new.ids
      saveRDS(old.ids, file = old.IDs.file)
    } else if (str_detect(old.IDs.file, "\\.csv")) {
      write.table(ids.dt, file = old.IDs.file, sep = ",", row.names = FALSE)
    } else if (str_detect(old.IDs.file, "\\.fas*t*a*$")) {
      seqinr::write.fasta(sequences = old.ids, names = new.ids)
    } else {
      rlang::abort(
        "Argument `old.IDs.file' must have an extension that is '.fa', '.fasta', '.csv', or '.rds'"
      )
    }
  }
  update_feature_names(mD = mD, new.names = new.ids) %>%
    return()
}
