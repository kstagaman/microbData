% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes_and_methods.R
\docType{class}
\name{microbData-class}
\alias{microbData-class}
\title{The class for microbData}
\description{
The class for microbData
}
\section{Slots}{

\describe{
\item{\code{Metadata}}{a \code{\link{data.table}} of sample metadata, must be keyed by sample name}

\item{\code{Abundances}}{a \code{\link{matrix}} of abundance counts by sample}

\item{\code{Assignments}}{a \code{\link{data.table}} of assignments/annotations for the features}

\item{\code{Phylogeny}}{a \code{\link{phylo}} class phylogeny of features}

\item{\code{Sample.names}}{a \code{character} vector of sample names/IDs}

\item{\code{Feature.names}}{a \code{character} vector of feature names/IDs}

\item{\code{Distance.matrices}}{= a \code{\link{dist}} object (or \code{list} of them) of distance metrics pairwise by sample}

\item{\code{Sample.col}}{= the name of the sample column in the \code{Metadata} slot (\code{character})}

\item{\code{Feature.col}}{= the name of the sample column in the \code{Assignments} slot (\code{character})}

\item{\code{Other.data}}{= a \code{list} of any other information to be associated with the microbData object}
}}

