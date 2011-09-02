\name{verifyHeader}
\alias{verifyHeader}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Verify DICOM Header Information
}
\description{
  Check specific DICOM header values against a user-specified list of
  fields and reference values.
}
\usage{
verifyHeader(hdr, reference, nr = NULL, nc = 4, out.names = c("Field", "Value", "Reference Value", "Result"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hdr}{is the DICOM header information in a list structure.}
  \item{reference}{is the path to a user-specified file containing the
    specific DICOM header fields and reference values.  See below for
    instructions on the format of this file.}
  \item{nr}{is the number of rows of the data frame (default =
    \code{NULL}).}
  \item{nc}{is the number of columns of the data frame (default =
    \code{4}).}
  \item{out.names}{are the column names of the output data frame
    (default provided).}
}
\details{
  TBD.
}
\value{
  A data frame is produced where each row is a DICOM header field.
}
%\references{}
\author{Brandon Whitcher \email{bjw34032@gsk.com}}
%\note{}
\seealso{\code{\link{qaqc}}}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }

