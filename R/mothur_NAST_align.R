#' @name mothur.NAST.align
#' @title NAST Alignment with \code{mothur}
#' @description Uses the \code{align.seqs} function from \code{mothur} (requires install: https://mothur.org/wiki/installation/)
#' @param fasta required; the path and name of the FASTA file containing the sequences to be aligned
#' @param template.file required; the path and name of the alignment file to be used. See https://mothur.org/wiki/alignment_database/ to download a database of your choice.
#' @param mothur.path character; path to \code{mothur} executable if it is not in your system PATH variable. Default is NULL.
#' @param output.dir character; path to output directory. Default is \code{getwd()}.
#' @param ncores integer; the number of cores for \code{mothur}.
#' @seealso \code{\link{generate.tool.command}}, \code{\link{system}}
#' @export

mothur.NAST.align <- function(fasta, template.file, mothur.path = NULL, output.dir = getwd(), ncores = 1) {
  if (is.null(mothur.path)) {
    sys.mothur.path <- Sys.which("mothur")
    if (!str_detect(sys.mothur.path, "mothur")) {
      rlang::abort(
        "Argument `mothur.path' is NULL, but the mothur executable is not in your PATH variable. Please add mothur to PATH or supply the path to the executable with `mothur.path'"
      )
    } else {
      mothur.path <- sys.mothur.path[1]
    }
  }
  cmd <- generate.tool.command(
    tool.path = mothur.path,
    function.name = "align.seqs",
    tool.syntax = "mothur",
    fasta = fasta,
    reference = template.file,
    flip = "t",
    keepdots = "t",
    processors = ncores,
    outputdir = output.dir
  )
  cat(paste("Executing command:", cmd), sep = "\n")
  system(cmd)
}
