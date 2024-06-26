% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copy.attr.R
\name{copy_attributes}
\alias{copy_attributes}
\alias{copy_attributes.default}
\alias{copy_attributes.generic_spct}
\alias{copy_attributes.generic_mspct}
\alias{copy_attributes.waveband}
\title{Copy attributes}
\usage{
copy_attributes(x, y, which, ...)

\method{copy_attributes}{default}(x, y, which = NULL, ...)

\method{copy_attributes}{generic_spct}(x, y, which = NULL, which.not = NULL, copy.class = FALSE, ...)

\method{copy_attributes}{generic_mspct}(x, y, which = NULL, which.not = NULL, copy.class = FALSE, ...)

\method{copy_attributes}{waveband}(x, y, which = NULL, ...)
}
\arguments{
\item{x, y}{R objects}

\item{which}{character Names of attributes to copy, if NULL all those
relevant according to the class of \code{x} is used as defaul,}

\item{...}{not used}

\item{which.not}{character Names of attributes not to be copied. The
names passed here are removed from the list for \code{which}, which
is most useful when we want to modify the default.}

\item{copy.class}{logical If TRUE class attributes are also copied.}
}
\value{
A copy of \code{y} with additional attributes set.
}
\description{
Copy attributes from \code{x} to \code{y}. Methods defined for spectral
and waveband objects of classes from package 'photobiology'.
}
\section{Methods (by class)}{
\itemize{
\item \code{copy_attributes(default)}: Default for generic function

\item \code{copy_attributes(generic_spct)}: 

\item \code{copy_attributes(generic_mspct)}: 

\item \code{copy_attributes(waveband)}: 

}}
