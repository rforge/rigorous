\name{qaqc}
\alias{qaqc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QAQC of DICOM Header Information}
\description{
  Quality assurance and quality control checks are applied to a
  collection of DICOM data sets.
}
\usage{
qaqc(subject, visit, exam, setup, pixelData = TRUE, domain = "/usr/cic/apps/bin/domain", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subject}{is a list of subject IDs (obtained from Domain).}
  \item{visit}{is a list of visit IDs (obtained from Domain).}
  \item{exam}{is a list of exam IDs (obtained from Domain).}
  \item{setup}{is a list of parameter values required for input/output
    of the QAQC process.  The following parameters are provided:
    \describe{
      \item{today}{is the date in YYYY-MM-DD format.}
      \item{domDate}{is the date in DD.MM.YYY format.}
      \item{study}{is a character string for the study ID.}
      \item{logFile}{is a character string for the log file.}
      \item{domFile}{is the path to the.}
      \item{dicom}{is the path to the Domain download directory.}
      \item{nifti}{is the path to the NIfTI directory.  All
	dicom-to-nifti conversions will be placed in this directory.}
      \item{shared.area}{is an alternative location for the HTML output
	(intended to be a Windows-based shared area).}
      \item{compare}{is the user-specified file for
	\code{\link{verifyHeader}}.}
      \item{emailList}{is an email address for the notification to be
	sent.}
      \item{mailText}{is a character string containg the body of the
	email message.}
      \item{css}{is the name of the CSS file.}
    }
  }
  \item{pixelData}{is a flag to determine if pixel data will be uploaded
    from the DICOM files (default = \code{TRUE}).}
  \item{domain}{is the path to the Domain executable.}
  \item{verbose}{is a flag to provide verbose output.}
  \item{\dots}{is for additional parameters to be sent to dicomSeparate().}
}
\details{
  TBD.
}
\value{
  No object is returned from this function.  An HTML file will be
  produced for each (subject,visit,exam) triplet and the corresponding
  paths to these files are emailed to the user.
}
%\references{}
\author{Brandon Whitcher \email{bjw34032@gsk.com}}
%\note{}
\seealso{\code{\link[oro.dicom]{dicom2nifti}},
  \code{\link[oro.dicom]{dicomSeparate}},
  \code{\link{verifyHeader}}}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
