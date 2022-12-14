% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_tool_command.R
\name{generate.tool.command}
\alias{generate.tool.command}
\alias{build.mothur.command}
\alias{build.py.bash.command}
\title{Generate tool command}
\usage{
generate.tool.command(
  tool.path,
  function.name = NULL,
  tool.syntax = "python",
  no.double.dashes = FALSE,
  ...
)

build.mothur.command(mothur.path, mothur.function, args)

build.py.bash.command(cmd.base, flags, no.double.dashes = FALSE)
}
\arguments{
\item{tool.path}{required; path to executable.}

\item{function.name}{character; name of function if building a mothur command. Default is NULL.}

\item{tool.syntax}{character; what syntax to use for building the tool command. Options include "python", "bash", or "mothur". Default is "python".}

\item{no.double.dashes}{logical; some command line tools don't use double dashes for flag names longer than 1 character. Setting this to true gives all flags a single dash prefix. Default is FALSE.}

\item{...}{other commands to pass to appropriate tool. Names must match short or long version found in that tool's help page. If the flag takes no argument in the tool, pass "flag" to the argument. E.g. for bash command \code{ls -l} you could run \code{generate.tool.command(tool = "ls", l = "flag")}.}
}
\description{
This function generates a character string of a command given a tool and its arguments.
}
\details{
This function is primarily run internally, but is available to the user as a wrapper for running other command line tools outside of R.
}
\seealso{
\code{\link{system}}, \code{\link{list2}}
}
