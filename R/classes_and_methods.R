################################################################################
# Taken from phyloseq:
#' S3 class placeholder definition (list) for phylogenetic trees.
#'
#' The ape package does not export a version of its \code{\link[ape]{phylo}}-class,
#' partly because it is not really defined formally anywhere.
#' Instead, it is an S3 class extended from the base class, \code{\link{list}} --
#' this is a very common and easy approach --
#' and proper behavior of any method taking an instance of this class
#' requires exact naming conventions for element names of the components.
#' The phyloseq package does not provide any validity checks that a given phylo
#' instance is valid (conforms to the conventions in the ape package). Yet.
#' If problems arise, this might be considered, and they could be defined
#' judiciously and within phyloseq.
#' Similarly, if a formal definition for the the phylo-class is ever exported
#' by ape, the current philosophy of phyloseq would be to remove this
#' internal definition and import the former. Note that there is still some
#' work going on for the phylobase package, which is addressing these same
#' exact issues for S4 phylogenetic tree interaction.
#' A very large number of packages (around 60 at my last count), depend on ape,
#' making it easily the de facto standard for representing phylogenetic trees in R;
#' and the phyloseq team would prefer to use any exported definitions from
#' the ape package if possible and available.
#'
#' @seealso
#' \code{\link[ape]{phylo}}
#'
#' @keywords internal
phylo <- structure(list(), class = "phylo")
################################################################################
# If this ever works
# @importClassesFrom ape phylo
################################################################################
#' An S4 placeholder of the main phylogenetic tree class from the ape package.
#'
#' See the \code{\link[ape]{ape}} package for details about this type of
#' representation of a phylogenetic tree.
#' It is used throughout the ape package.
#'
#' @seealso \code{\link[ape]{phylo}}, \code{\link{setOldClass}}
#'
#' @name phylo-class
#' @rdname phylo-class
#' @exportClass phylo
setOldClass("phylo")
################################################################################
# Taken from phyloseq:
# Use setClassUnion to define the unholy NULL-data union as a virtual class.
# This is a way of dealing with the expected scenarios in which one or more of
# the component data classes is not available, in which case NULL will be used
# instead.
################################################################################
setOldClass("dist")
#' @keywords internal
setClassUnion("data.tableOrNULL", c("data.table", "NULL"))
#' @keywords internal
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
#' @keywords internal
setClassUnion("phyloOrNULL", c("phylo", "NULL"))
#' @keywords internal
setClassUnion("listOrNULL", c("list", "NULL"))
library
#' @keywords internal
setClassUnion("listOrDistOrNULL", c("list", "dist", "NULL"))
#' @keywords internal
setClassUnion("characterOrNULL", c("character", "NULL"))
################################################################################
#' @title The class for microbData
#' @name microbData-class
#' @rdname microbData-class
#' @exportClass microbData

setClass(
  Class = "microbData",
  slots = c(
    Metadata = "data.tableOrNULL",
    Abundances = "matrixOrNULL",
    Assignments = "data.tableOrNULL",
    Phylogeny = "phyloOrNULL",
    Sample.names = "characterOrNULL",
    Feature.names = "characterOrNULL",
    Distance.matrices = "listOrDistOrNULL",
    Sample.col = "characterOrNULL",
    Feature.col = "characterOrNULL",
    Other.data = "listOrNULL"
  )
)
################################################################################

#' @name show-microbData
#' @title Method extensions to `show` for microbData object
#' @seealso \code{\link[methods]{show}}
#' @inheritParams methods::show
#' @rdname show-methods
#' @export

setMethod(
  "show",
  "microbData",
  function(object) {
    display.vector <- function(vec) {
      if (length(vec) > 4) {
        to.show <- paste0(paste(vec[1:3], collapse = ", "), ",..., ", tail(vec, 1))
      } else {
        to.show <- paste(vec, collapse = ", ")
      }
      return(to.show)
    }
    display.table <- function(tbl, n = 6) {
      n.col <- min(c(n, ncol(tbl)))
      n.row <- min(c(n, nrow(tbl)))
      print(tbl[1:n.row, 1:n.col])
    }
    display.list <- function(lst) {
      to.show <- paste("List of", length(lst))
      if (is.null(names(lst))) {
        to.show <- paste0(to.show, "; unnamed")
      } else {
        to.show <- paste0(to.show, ";\n\tnames: ", display.vector(names(lst)))
      }
      return(to.show)
    }

    cat(
      paste(
        "Sample Metadata:",
        nrow(object@Metadata),
        "samples with",
        ncol(object@Metadata) - 1,
        "covariates"
      ),
      sep = "\n"
    )
    cat("Preview:", sep = "\n")
    display.table(object@Metadata, 4)
    cat("", sep = "\n")

    cat(paste("Feature Abundances:", ncol(object@Abundances), "features"), sep = "\n")
    cat("Preview:", sep = "\n")
    display.table(object@Abundances, 4)
    cat("", sep = "\n")

    if (!is.null(object@Assignments)) {
      cat(
        paste(
          "Feature Assignments:",
          ncol(object@Assignments) - 1,
          "levels assigned"
        ),
        sep = "\n"
      )
      cat("Preview:", sep = "\n")
      display.table(object@Assignments, 4)
      cat("", sep = "\n")
    }

    if (!is.null(object@Phylogeny)) {
      cat(
        paste(
          "Phylogentic Tree:",
          ncol(object@Abundances), "tips",
          object@Phylogeny$Nnode, "internal nodes"
        ),
        sep = "\n"
      )
      cat("", sep = "\n")
    }

    if (!is.null(object@Distance.matrices)) {
      cat(
        paste("Distance Matrices:", display.list(object@Distance.matrices)),
        sep = "\n"
      )
      cat("", sep = "\n")
    }

    cat(paste("Sample Names:", display.vector(object@Sample.names)), sep = "\n")
    cat("", sep = "\n")
    cat(paste("Feature Names:", display.vector(object@Feature.names)), sep = "\n")
    cat("", sep = "\n")

    if (!is.null(object@Other.data)) {
      cat("Other Data:", sep = "\n")
      to.print <- lapply(seq_along(object@Other.data), function(i) {
        item <- names(object@Other.data)[i]
        vals <- object@Other.data[[i]]
        item.line <- paste0("  ", item, ":")
        if (is.vector(vals) & !is.list(vals)) {
          cat(paste(item.line, display.vector(vals)), sep = "\n")
        } else if (is.list(vals)) {
          cat(paste(item.line, display.list(vals)), sep = "\n")
        } else if (class(vals) == "function") {
          cat(paste(item.line, paste(deparse(vals), collapse = "")), sep = "\n")
        } else {
          cat(paste(item.line, "1", paste(class(vals), collapse = "/")), sep = "\n")
        }
      })
    }
  }
)
