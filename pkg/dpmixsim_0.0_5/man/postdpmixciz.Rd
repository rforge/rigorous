\name{postdpmixciz}
\alias{postdpmixciz}
\concept{Summary statistics and cluster estimation}

\title{Summary statistics and cluster estimation}

\description{
\code{postdpmixciz} computes post-simulation summary statistics,
  and estimates cluster partition.
}
\usage{postdpmixciz(x, res, kmax=30, rec=300, ngrid=200, plot=TRUE)}

\arguments{
  \item{x}{ data used in the simulation }
  \item{kmax}{ maximum number of clusters }
  \item{res}{ output of the MCMC simulation }
  \item{rec}{ number of recorded iteration steps }
  \item{ngrid}{ dimension of the grid used in density estimation }
  \item{plot}{ logical variable to omit plots (default = \samp{TRUE} }
}
\value{
  \item{z}{cluster partition estimation}
}
\references{ 
Adelino Ferreira da Silva, A Dirichlet process mixture model for
 brain MRI tissue classification,
 \emph{Medical Image Analysis} 11 (2007) 169-182.

Adelino Ferreira da Silva, Bayesian mixture models of variable
 dimension for image segmentation,
 \emph{Comput. Methods Programs Biomed.} 94 (2009) 1-14.
}

\author{A. Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia,\cr
  \email{afs@fct.unl.pt}.}

\seealso{
  \code{\link{dpmixsim}}
}

\examples{
\dontrun{
## Example: MRI brain image segmentation
  slicedata <- readsliceimg(fbase="t1_pn3_rf0", swap=FALSE)
  image(slicedata$niislice, col=gray((0:255)/256), main="original image")
  x0 <- premask(slicedata, subsamp=TRUE)
  x  <- prescale(x0) 
  rec <- 3000
  res <- dpmixsim(x, M=1, a=1, b=2, upalpha=1, maxiter=4000,
      rec=rec, nclinit=8)
  ## post-simulation
  ngrid <- 200
  z <- postdpmixciz(x, res=res, rec=rec, ngrid=ngrid, plot=TRUE)
}}

\keyword{utilities}
