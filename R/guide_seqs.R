#' @name guide.seqs
#' @title Dealing with Guide Seqs
#' @description Functions for add and remove guide seqs from FASTAs, alignments, and trees.
#' @param asv.seqs.file required; FASTA file containing ASV sequences.
#' @param guide.seqs.file required; FASTA file containing guide sequences.
#' @param joined.file required; FASTA file to write concatenated sequences to.
#' @param keep.guide.ids logical; whether to return a vector containing guide sequence IDs for later removal from output. Default is TRUE
#' @param guide.seq.ids required; a charater vector of the guide sequence IDs, usually from running `add.guide.seqs`.
#' @param full.tree.file required; the file containing the phylogenetic tree containing both ASV and guide sequence IDs in Newick format.
#' @param asvs.only.tree.file required; the file to write the pruned tree of just ASV IDs to. Will be written in Newick format.
#' @seealso \code{\link{generate.tool.command}}, \code{\link{system}}

#' @rdname guide.seqs
#' @name add.guide.seqs
#' @export

add.guide.seqs <- function(asv.seqs.file, guide.seqs.file, joined.file, keep.guide.ids = TRUE) {
  cmd <- paste("cat", asv.seqs.file, guide.seqs.file, ">", joined.file)
  cat(paste("Executing command:", cmd), sep = "\n")
  system(cmd)
  if (keep.guide.ids) {
    readLines(con = guide.seqs.file) %>%
      str_subset("^>") %>%
      str_remove("^> *") %>%
      return()
  }
}

#' @rdname guide.seqs
#' @name remove.guide.seqs
#' @export

remove.guide.seqs <- function(guide.seq.ids, full.tree.file, asvs.only.tree.file) {
  full.tree <- read_tree(full.tree.file)
  full.tree.rooted <- phangorn::midpoint(full.tree)

  asvs.tree.rooted <- ape::drop.tip(phy = full.tree.rooted, tip = guide.seq.ids)
  write.tree(asvs.tree.rooted, file = asvs.only.tree.file)
}
