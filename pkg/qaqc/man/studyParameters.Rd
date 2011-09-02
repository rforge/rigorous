\name{studyParameters}
\alias{studyParameters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial Parameters for QAQC}
\description{
  Define a number of important input parameters for the QAQC process
  using text-based input from the command line.
}
\usage{
studyParameters(argList, date2process = NULL, emailList = NULL, study = NULL, inPath = file.path("", "cic_data", "Oncology"), outPath = file.path("", "scratch"), shared.area = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{argList}{Text string from the \dQuote{--args} call to \code{Rscript}.}
  \item{date2process}{Date (in DD-MM-YYYY format) to be processes.  If
    no date is supplied, then today's date is used.}
  \item{emailList}{Text string of valid email addresses separated by commas.}
  \item{study}{Study ID.}
  \item{inPath}{Path name for input parameters.}
  \item{outPath}{Path name for output parameters.}
  \item{shared.area}{Path name for the Windows \sQuote{shared area}
    (default = \code{NULL} for local directory).}
}
%\details{}
\value{
  List with the following elements:
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
%\references{}
\author{Brandon Whitcher \email{bjw34032@gmail.com}}
%\note{}
\seealso{
  \code{\link{commandArgs}}, \code{\link{qaqc}}
}
\examples{
args <- "--args study=\"MET111645\" date2process=\"22-03-2010\" inPath=\"/cic_data/Oncology\" outPath=\"/scratch\" sharedArea=\"//hlwsfshr10/sa200400214\" emailList=\"bjw34032@gsk.com,cf984373@gsk.com\""
unlist(studyParameters(args))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
