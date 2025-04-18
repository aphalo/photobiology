% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.utils.r
\name{drop_user_cols}
\alias{drop_user_cols}
\alias{drop_user_cols.default}
\alias{drop_user_cols.generic_spct}
\alias{drop_user_cols.source_spct}
\alias{drop_user_cols.response_spct}
\alias{drop_user_cols.object_spct}
\alias{drop_user_cols.filter_spct}
\alias{drop_user_cols.reflector_spct}
\alias{drop_user_cols.solute_spct}
\alias{drop_user_cols.chroma_spct}
\alias{drop_user_cols.calibration_spct}
\alias{drop_user_cols.cps_spct}
\alias{drop_user_cols.raw_spct}
\alias{drop_user_cols.generic_mspct}
\title{Drop user columns}
\usage{
drop_user_cols(x, keep.also, ...)

\method{drop_user_cols}{default}(x, keep.also = NULL, ...)

\method{drop_user_cols}{generic_spct}(x, keep.also, ...)

\method{drop_user_cols}{source_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{response_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{object_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{filter_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{reflector_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{solute_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{chroma_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{calibration_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{cps_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{raw_spct}(x, keep.also = NULL, ...)

\method{drop_user_cols}{generic_mspct}(x, keep.also = NULL, ...)
}
\arguments{
\item{x}{An R object}

\item{keep.also}{character Additionlal columns to preserve.}

\item{...}{needed to allow derivation.}
}
\value{
A copy of \code{x} possibly with some columns removed.
}
\description{
Remove from spectral object additional columns that are user defined.
}
\section{Methods (by class)}{
\itemize{
\item \code{drop_user_cols(default)}: 

\item \code{drop_user_cols(generic_spct)}: 

\item \code{drop_user_cols(source_spct)}: 

\item \code{drop_user_cols(response_spct)}: 

\item \code{drop_user_cols(object_spct)}: 

\item \code{drop_user_cols(filter_spct)}: 

\item \code{drop_user_cols(reflector_spct)}: 

\item \code{drop_user_cols(solute_spct)}: 

\item \code{drop_user_cols(chroma_spct)}: 

\item \code{drop_user_cols(calibration_spct)}: 

\item \code{drop_user_cols(cps_spct)}: 

\item \code{drop_user_cols(raw_spct)}: 

\item \code{drop_user_cols(generic_mspct)}: 

}}
\seealso{
Other experimental utility functions: 
\code{\link{collect2mspct}()},
\code{\link{thin_wl}()},
\code{\link{uncollect2spct}()}
}
\concept{experimental utility functions}
