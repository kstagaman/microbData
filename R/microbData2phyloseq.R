#' @name microbData2phyloseq
#' @title MicrobData to Phyloseq
#' @description Convert a \code{microbData} object to a phyloseq object
#' @param mD required; the \code{microbData} object to be converted to a phyloseq object.
#' @seealso \code{\link[phyloseq]{phyloseq}}
#' @export

microbData2phyloseq <- function(mD) {
  require(phyloseq)
  abundances <- mD@Abundances
  metadata <- as.data.frame(mD@Metadata)
  row.names(metadata) <- mD@Metadata[[mD@Sample.col]]
  metadata <- metadata[rownames(abundances), ]

  ps <- phyloseq(
    sample_data(metadata),
    otu_table(abundances)
  )
  if (!is.null(mD@Phylogeny)) {
    phy_tree(ps) <- mD@Phylogeny
  }
  if (!is.null(mD@Assignments)) {
    assignments <- as.data.frame(mD@Assignments)
    row.names(assignments) <- mD@Assignments[[mD@Feature.col]]
    tax_table(ps) <- assignments
  }
  return(ps)
}
