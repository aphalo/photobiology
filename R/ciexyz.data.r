#' @title CIE xyz chromaticity coordinates 2 deg data
#'
#' @description A dataset containing wavelengths at a 1 nm interval (390 nm to
#'   830 nm) and the corresponding x, y, and z chromaticity coordinates.
#'   According to proposed CIE 2006 standard. Original data from
#'   \url{http://www.cvrl.org/} downloaded on 2014-04-28 The variables are as
#'   follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item x
#'   \item y
#'   \item z }
#'
#' @author CIE
#'
#' @docType data
#' @keywords datasets
#' @format A chroma_spct object with 441 rows and 4 variables
#' @family Visual response data examples
#'
#' @examples
#' ciexyzCC2.spct
#'
"ciexyzCC2.spct"

#' @title CIE xyz chromaticity coordinates (CC) 10 deg data
#'
#' @description A dataset containing wavelengths at a 1 nm interval (390 nm to
#'   830 nm) and the corresponding x, y, and z chromaticity coordinates. Derived
#'   from proposed CIE 2006 standard. Original data from
#'   \url{http://www.cvrl.org/} downloaded on 2014-04-29 The variables are as
#'   follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item x
#'   \item y
#'   \item z }
#'
#' @author CIE
#'
#' @docType data
#' @keywords datasets
#' @format A chroma_spct object with 441 rows and 4 variables
#' @family Visual response data examples
#'
#' @examples
#' ciexyzCC10.spct
#'
"ciexyzCC10.spct"

#' @title Linear energy CIE xyz colour matching function (CMF) 2 deg data
#'
#' @description A dataset containing wavelengths at a 1 nm interval (390 nm to
#'   830 nm) and the corresponding x, y, and z 2 degrees CMF values. Derived
#'   from proposed CIE 2006 standard. Original data from
#'   \url{http://www.cvrl.org/} downloaded on 2014-04-29 The variables are as
#'   follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item x
#'   \item y
#'   \item z }
#'
#' @author CIE
#'
#' @docType data
#' @keywords datasets
#' @format A chroma_spct object with 441 rows and 4 variables
#' @family Visual response data examples
#'
#' @examples
#' ciexyzCMF2.spct
#'
"ciexyzCMF2.spct"

#' @title Linear energy CIE xyz colour matching function (CMF) 10 deg data
#'
#' @description A dataset containing wavelengths at a 1 nm interval (390 nm to
#'   830 nm) and the corresponding x, y, and z 10 degrees CMF values. Derived
#'   from proposed CIE 2006 standard. Original data from
#'   \url{http://www.cvrl.org/} downloaded on 2014-04-29 The variables are as
#'   follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item x
#'   \item y
#'   \item z }
#'
#' @author CIE
#'
#' @docType data
#' @keywords datasets
#' @format A chroma_spct object with 441 rows and 4 variables
#' @family Visual response data examples
#'
#' @examples
#' ciexyzCMF10.spct
#'
"ciexyzCMF10.spct"

#' @title Linear energy CIE 2008 luminous efficiency function 10 deg data
#'
#' @description A dataset containing wavelengths at a 1 nm interval (390 nm to
#'   830 nm) and the corresponding response values for a 10 degrees target.
#'   Original data from \url{http://www.cvrl.org/} downloaded on 2014-04-29 The
#'   variables are as follows:
#'
#' \itemize{
#'   \item w.length (nm)
#'   \item x
#'   \item y
#'   \item z }
#'
#' @author CIE
#'
#' @docType data
#' @keywords datasets
#' @format A chroma_spct object with 441 rows and 4 variables
#' @family Visual response data examples
#'
#' @examples
#' ciev10.spct
#'
"ciev10.spct"

#' @title Linear energy CIE 2008 luminous efficiency function 2 deg data
#'
#' @description A dataset containing wavelengths at a 1 nm interval (390 nm to
#'   830 nm) and the corresponding response values for a 2 degrees target.
#'   Original data from \url{http://www.cvrl.org/} downloaded on 2014-04-29 The
#'   variables are as follows:
#'
#' @details
#' \itemize{
#'   \item w.length (nm)
#'   \item x
#'   \item y
#'   \item z }
#'
#' @author CIE
#'
#' @docType data
#' @keywords datasets
#' @format A chroma_spct object with 441 rows and 4 variables
#' @family Visual response data examples
#'
#' @examples
#' ciev2.spct
#'
"ciev2.spct"

#' @title Ten-degree cone fundaamentals
#'
#' @description A dataset containing wavelengths at a 1 nm interval (390 nm to
#'   830 nm) and the corresponding response values for a 2 degrees target.
#'   Original data from \url{http://www.cvrl.org/} downloaded on 2014-04-29 The
#'   variables are as follows:
#'
#' @details
#' \itemize{
#'   \item w.length (nm)
#'   \item x
#'   \item y
#'   \item z }
#'
#' @author CIE
#'
#' @return A \code{chroma_spct} object.
#' @docType data
#' @keywords datasets
#' @format A chroma_spct object with 440 rows and 4 variables
#' @family Visual response data examples
#'
#' @examples
#' cone_fundamentals10.spct
#'
"cone_fundamentals10.spct"


#' @rdname cone_fundamentals10.spct
#'
#' @return A \code{response_mspct} object containing the same data in three
#' \code{response_spct} objects.
#'
"cone_fundamentals10.mspct"
