% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.metadata.r
\name{trimInstrDesc}
\alias{trimInstrDesc}
\title{Trim the "instr.desc" attribute}
\usage{
trimInstrDesc(
  x,
  fields = c("time", "spectrometer.name", "spectrometer.sn", "bench.grating",
    "bench.slit", "entrance.optics")
)
}
\arguments{
\item{x}{a generic_spct object}

\item{fields}{a character vector with the names of the fields to keep,
or if first member is `"-"`, the names of fields to delete; "*" as
first member of the vector makes the function a no-op, leaving the spectrum
object unaltered.}
}
\value{
x
}
\description{
Function to trim the "instr.desc" attribute of an existing generic_spct
object, discarding all fields except for `spectrometer.name`,
`spectrometer.sn`, `bench.grating`, `bench.slit`, and calibration name.
}
\note{
This function alters x itself by reference and in addition
  returns x invisibly. If x is not a generic_spct object, x is not
  modified.
}
\seealso{
Other measurement metadata functions: 
\code{\link{add_attr2tb}()},
\code{\link{getFilterProperties}()},
\code{\link{getHowMeasured}()},
\code{\link{getInstrDesc}()},
\code{\link{getInstrSettings}()},
\code{\link{getSoluteProperties}()},
\code{\link{getWhatMeasured}()},
\code{\link{getWhenMeasured}()},
\code{\link{getWhereMeasured}()},
\code{\link{get_attributes}()},
\code{\link{isValidInstrDesc}()},
\code{\link{isValidInstrSettings}()},
\code{\link{select_spct_attributes}()},
\code{\link{setFilterProperties}()},
\code{\link{setHowMeasured}()},
\code{\link{setInstrDesc}()},
\code{\link{setInstrSettings}()},
\code{\link{setSoluteProperties}()},
\code{\link{setWhatMeasured}()},
\code{\link{setWhenMeasured}()},
\code{\link{setWhereMeasured}()},
\code{\link{spct_attr2tb}()},
\code{\link{spct_metadata}()},
\code{\link{subset_attributes}()},
\code{\link{trimInstrSettings}()}
}
\concept{measurement metadata functions}
