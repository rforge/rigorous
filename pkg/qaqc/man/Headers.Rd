\name{compareHeader}
\alias{compareHeader}
\alias{multipleHeader}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test DICOM Headers}
\description{
  ...
}
\usage{
compareHeader(dcm, dicomHeader, strings, htmlfile, numeric = FALSE)
multipleHeader(dcm, dicomHeader, htmlfile, numeric = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dcm}{DICOM object.}
  \item{dicomHeader}{Text string denoting the DICOM header field.}
  \item{strings}{Vector of text strings for comparison.}
  \item{htmlfile}{'hwriter' object.}
  \item{numeric}{Logical parameter (default = \code{FALSE}) for \code{matchHeader}.}
}
%\details{}
\value{
  A expression to be evaluated is produced, either a
  \code{expression(next)} or an assignment; e.g., \code{expression(a <- b)}. 
}
%\references{}
\author{Brandon Whitcher \email{bjw34032@gsk.com}}
%\note{}
\seealso{\code{\link{matchHeader}}}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
