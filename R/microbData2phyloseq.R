#' @name microbData2phyloseq
#' @title MicrobData to Phyloseq
#' @description Convert a \code{microbData} object to a phyloseq object
#' @param mD required; the \code{microbData} object to be converted to a phyloseq object.
#' @seealso \code{\link[phyloseq]{phyloseq}}
#' @export

microbData2phyloseq <- function (mD) {
  require(phyloseq)
  abundances <- mD@Abundances
  metadata <- as.data.frame(mD@Metadata)
  row.names(metadata) <- mD@Metadata[[mD@Sample.col]]
  metadata <- metadata[rownames(abundances), ]
  ps <- phyloseq(
    sample_data(metadata), 
    otu_table(abundances, taxa_are_rows = FALSE)
  )
  if (!is.null(mD@Phylogeny)) {
    phy_tree(ps) <- mD@Phylogeny
  }
  if (!is.null(mD@Assignments)) {
    assignments <- as.matrix(mD@Assignments)
    rownames(assignments) <- mD@Assignments[[mD@Feature.col]]
    identical(sort(rownames(assignments)), sort(taxa_names(ps)))
    tax_table(ps) <- assignments[taxa_names(ps), ]
  }
  return(ps)
}

