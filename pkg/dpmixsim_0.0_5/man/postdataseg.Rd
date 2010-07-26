\name{postdataseg}
\alias{postdataseg}
\concept{data segmentation}

\title{Data segmentation}

\description{
  \code{postdataseg} 
  performs data segmentation based on labelled cluster estimates.
}
\usage{postdataseg(x, z, ngrid, dbg=FALSE)}
\arguments{
  \item{x}{ full-sized scaled image data prepared by \code{premask} }
  \item{z}{ cluster labels produced by \code{postdpmixciz} }
  \item{ngrid}{ dimension of the grid used in estimation }
  \item{dbg}{ logical variable to show debugging output (default = \samp{FALSE}) }
}

\details{
  Once the distributions of the indicator variables \eqn{z_i} are
  calculated we can separate the components of the mixture.
  Individual components are selected according to the most probable
  \eqn{z_i} value in a given region of the distributional space,
  leading to a partition of this space into regions.
  Intensity threshold values are associated with the partition of the
  distributional space to drive the image segmentation.
  In brief, the partition of the distributional space induced by the
  \eqn{z} values is used to segment the data space.
  From a computational point of view, the use of these two separate 
  spaces enables us to optimize the MCMC implementation.
}

\value{
  \item{cx}{vector of image cluster values}
}

\author{A. Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia,\cr
  \email{afs@fct.unl.pt}.}

\seealso{
  \code{\link{dpmixsim}},
  \code{\link{readsliceimg}},
  \code{\link{premask}},
  \code{\link{postdpmixciz}}
}

\examples{
\dontrun{
## see Example 2 in dpmixsim.
}}

\keyword{utilities}
