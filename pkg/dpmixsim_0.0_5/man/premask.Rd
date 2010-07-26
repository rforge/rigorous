\name{premask}
\alias{premask}
\concept{Mask mri data}
\title{ Data preparation }
\description{
  \code{premask} 
  applies a pre-defined mask to a MRI slice in order to select regions of interest
  (ROIs) for processing 
}
\usage{
  premask(slicedata, subsamp=TRUE)
}
\arguments{
  \item{slicedata}{ list as output by \code{read.sliceimg} }
  \item{subsamp}{ logical variable; if \samp{TRUE} a downsampled image by a factor
  of 2 is used in the MCMC simulation, otherwise the full-sized image is taken.
  After parameter estimation, the full-sized image should be used for clustering
  and image segmentation.  The use of downsampled images can substantially
  reduce runtime, with little quality degradation. }  
}
\value{
  \item{xv}{ processed data vector }
}
\author{A. Ferreira da Silva, Universidade Nova de Lisboa,
Faculdade de Ciencias e Tecnologia,\cr
\email{afs@fct.unl.pt}.
}

\seealso{
  \code{\link{dpmixsim}},
  \code{\link{readsliceimg}}
}

\examples{
\dontrun{
  slicedata <- readsliceimg(fbase="t1_pn3_rf0", swap=FALSE)
  x0 <- premask(slicedata, subsamp=TRUE)
  x  <- prescale(x0) 
  print(str(x))
}}

\keyword{utilities}
