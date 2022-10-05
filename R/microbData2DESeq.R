#' @name microbData2DESeq
#' @title MicrobData to DESeq2
#' @description Convert a \code{microbData} object to a DESeqDataSet object.
#' @param mD required; the \code{microbData} obejct to be converted.
#' @param expt.design required; from \code{\link[DESeq2]{DESeqDataSetFromMatrix}}: a formula or matrix. the formula expresses how the counts for each gene depend on the variables in colData. Many R formula are valid, including designs with multiple variables, e.g., ~ group + condition, and designs with interactions, e.g., ~ genotype + treatment + genotype:treatment. See results for a variety of designs and how to extract results tables. By default, the functions in this package will use the last variable in the formula for building results tables and plotting. ~ 1 can be used for no design, although users need to remember to switch to another design for differential testing.
#' @seealso \code{\link[DESeq2]{DESeqDataSetFromMatrix}}
#' @export

microbData2DEseq <- function(mD, expt.design) {
  require(DESeq2)
  DESeq2::DESeqDataSetFromMatrix(
    countData = t(mD@Abundances),
    colData = mD@Metadata,
    design = expt.design
  ) %>%
    return()
}
