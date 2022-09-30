#' @name generate.tool.command
#' @title Generate tool command
#' @description This function generates a character string of a command given a tool and its arguments.
#' @param tool.path required; path to executable.
#' @param function.name character; name of function if building a mothur command. Default is NULL.
#' @param tool.syntax character; what syntax to use for building the tool command. Options include "python", "bash", or "mothur". Default is "python".
#' @param no.double.dashes logical; some command line tools don't use double dashes for flag names longer than 1 character. Setting this to true gives all flags a single dash prefix. Default is FALSE.
#' @param ... other commands to pass to appropriate tool. Names must match short or long version found in that tool's help page. If the flag takes no argument in the tool, pass "flag" to the argument. E.g. for bash command `ls -l` you could run `generate.tool.command(tool = "ls", l = "flag")`.
#' @details This function is primarily run internally, but is available to the user as a wrapper for running other command line tools outside of R.
#' @seealso \code{\link{system}}, \code{\link{list2}}
#' @export

generate.tool.command <- function(
    tool.path,
    function.name = NULL,
    tool.syntax = "python",
    no.double.dashes = FALSE,
    ...
    ) {
  require(magrittr)
  require(stringr)

  if (is.null(function.name) & tool.syntax == "mothur") {
    rlang::abort("If building a mothur command, please provide the name of the function")
  }

  vargs <- rlang::list2(...)
  vargs <- vargs[!sapply(vargs, is.null)]

  if (tool.syntax %in% c("bash", "python")) {
    cmd <- build.py.bash.command(cmd.base = tool.path, flags = vargs, no.double.dashes = no.double.dashes)
  } else if (tool.syntax == "mothur") {
    cmd <- build.mothur.command(mothur.path = tool.path, mothur.function = function.name, args = vargs)
  }
  return(cmd)
}

#' @rdname generate.tool.command
#' @export

build.mothur.command <- function(mothur.path, mothur.function, args) {
  cmd.base <- paste0(mothur.path, " \"#", mothur.function, "(")
  cmd.args <- sapply(seq_along(args), function(arg) {
    paste0(names(args)[[arg]], "=", args[[arg]]) %>%
      return()
  }) %>% paste(collapse = ", ")
  paste(cmd.base, cmd.args, ")\"") %>% return()
}

#' @rdname generate.tool.command
#' @export

build.py.bash.command <- function(cmd.base, flags, no.double.dashes = FALSE) {
  cmd.args <- sapply(seq_along(flags), function(arg) {
    arg.name <- names(flags)[[arg]]
    arg.val <- ifelse(flags[[arg]] == "flag", "", flags[[arg]])
    if (no.double.dashes) {
      dashes <- 1
    } else {
      dashes <- ifelse(nchar(arg.name) > 1, 2, 1)
    }
    paste0(
      paste(rep("-", dashes), collapse = ""),
      str_replace_all(arg.name, "\\.", "-"),
      ifelse(flags[[arg]] == "flag", "", " "),
      arg.val
    ) %>% return()
  }) %>% paste(collapse = " ")
  paste(cmd.base, cmd.args) %>% return()
}
