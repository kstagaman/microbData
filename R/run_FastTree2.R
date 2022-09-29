#' @name run.FastTree2
#' @title Build a Phylogenetic Tree with FastTree2
#' @description Build an approximately-maximum-likelihood with FastTree2 (http://www.microbesonline.org/fasttree/)
#' @param aligned.seqs.file required; FASTA file containing aligned sequences to build tree with.
#' @param output.file required; file name to write output, should have the '.nwk' extension. If not, it will be added.
#' @param fasttree.path haracter; path to `mothur` executable if it is not in your system PATH variable. Default is NULL.
#' @seealso \code{\link{generate.tool.command}}, \code{\link{system}}
#' @export
#' @examples
#'

run.FastTree2 <- function(aligned.seqs.file, output.file, fasttree.path = NULL, log.file = NULL, ncores = 1) {
  if (ncores < 1) {
    rlang::abort("Argument `ncores' must be an integer ≥1")
  } else {
    ncores <- as.integer(ncores)
  }
  if (is.null(fasttree.path)) {
    sys.fasttree.path <- ifelse(ncores == 1, Sys.which("FastTree"), Sys.which("FastTreeMP"))
    if (!str_detect(sys.fasttree.path, "FastTree")) {
      rlang::abort(
        "Argument `fasttree.path' is NULL, but the FastTree/FastTreeMP executable is not in your PATH variable. Please add FastTree/FastTreeMP to PATH or supply the path to the executable with `fasttree.path'"
      )
    } else {
      fasttree.path <- sys.fasttree.path[1]
    }
  }
  if (!str_detect(output.file, "\\.nwk$")) { output.file <- paste0(output.file, ".nwk") }

  base.cmd <- generate.tool.command(
    tool.path = fasttree.path,
    tool.syntax = "bash",
    nt = "flag",
    nosupport = "flag",
    quote = "flag",
    gtr = "flag",
    gamma = "flag",
    log = log.file,
    no.double.dashes = TRUE
  ) %>%
    paste(aligned.seqs.file, ">", output.file)
  cmd <- ifelse(ncores > 1, paste0("export OMP_NUM_THREADS=", ncores, "; ", base.cmd), base.cmd)
  cat(paste("Executing command:", cmd), sep = "\n")
  system(cmd)
}
