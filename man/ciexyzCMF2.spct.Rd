% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ciexyz.data.r
\docType{data}
\name{ciexyzCMF2.spct}
\alias{ciexyzCMF2.spct}
\title{Linear energy CIE xyz colour matching function (CMF) 2 deg data}
\format{
A chroma_spct object with 441 rows and 4 variables
}
\usage{
ciexyzCMF2.spct
}
\description{
A dataset containing wavelengths at a 1 nm interval (390 nm to
  830 nm) and the corresponding x, y, and z 2 degrees CMF values. Derived
  from proposed CIE 2006 standard. Original data from
  \url{http://www.cvrl.org/} downloaded on 2014-04-29 The variables are as
  follows:

\itemize{
  \item w.length (nm)
  \item x
  \item y
  \item z }
}
\note{
These data are not from the official CIE on-line distribution but
  are retained for backwards compatibility. It is recommended to download
  the latest version from \url{https://cie.co.at/data-tables}.
}
\examples{
ciexyzCMF2.spct

}
\seealso{
Other Visual response data examples: 
\code{\link{beesxyzCMF.spct}},
\code{\link{ciev10.spct}},
\code{\link{ciev2.spct}},
\code{\link{ciexyzCC10.spct}},
\code{\link{ciexyzCC2.spct}},
\code{\link{ciexyzCMF10.spct}},
\code{\link{cone_fundamentals10.spct}}
}
\author{
CIE
}
\concept{Visual response data examples}
\keyword{datasets}
