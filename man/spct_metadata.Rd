% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spct.metadata.r
\name{spct_metadata}
\alias{spct_metadata}
\title{Access metadata}
\usage{
spct_metadata(
  x,
  col.names = NULL,
  idx = "spct.idx",
  na.rm = is.null(col.names),
  unnest = TRUE
)
}
\arguments{
\item{x}{generic_mspct or generic_spct Any collection of spectra or spectrum.}

\item{col.names}{named character vector Name(s) of column(s) to create.}

\item{idx}{character Name of the column with the names of the members of the
collection of spectra.}

\item{na.rm}{logical Flag controlling deletion of columns containing only NA
values.}

\item{unnest}{logical Flag controlling if metadata attributes that are lists
of values should be returned in a list column or in separate columns.}
}
\value{
A tibble With the metadata attributes and an index column.
}
\description{
Return metadata attributes from a single spectrum or a collection of spectra
as a tibble.
}
\details{
Attributes are returned as columns in a tibble. If the argument to
  \code{col.names} is a named vector, with the names of members matching the
  names of attributes, then the values are used as names for the columns
  created. This permits setting any valid name for the new columns. If the
  vector passed to \code{col.names} has no names, then the values are
  interpreted as the names of the attributes to add, and also used as names
  for the new columns.

  Some metadata values are stored in lists or data frames, these can be
  returned as a list columns or the individual fields unnested into separate
  columns.
}
\examples{

my.mspct <- source_mspct(list(sun1 = sun.spct, sun2 = sun.spct * 2))

spct_metadata(my.mspct)

spct_metadata(sun.spct)

spct_metadata(my.mspct, na.rm = TRUE)

spct_metadata(sun.spct, na.rm = TRUE)

spct_metadata(my.mspct, col.names = c(geocode = "geo", "instr.desc"))

spct_metadata(sun.spct, col.names = c(geocode = "geo", "instr.desc"))

spct_metadata(sun.spct, col.names = "where.measured")$where.measured

}
\seealso{
\code{\link{add_attr2tb}} for more details.

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
\code{\link{subset_attributes}()},
\code{\link{trimInstrDesc}()},
\code{\link{trimInstrSettings}()}
}
\concept{measurement metadata functions}
