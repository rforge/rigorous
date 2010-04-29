\name{postimgcomps}
\alias{postimgcomps}
\concept{Display cluster components}
\title{Display cluster components}
\description{
  \code{postimgcomps} displays the components of the segmented image
  with the estimated number of components
}

\usage{postimgcomps(mask, cx)}
\arguments{
  \item{mask}{ scaled masked full-sized image data prepared by
  \code{premask} } \item{cx}{ data segmentation prepared by
  \code{postdataseg} }
}

\details{
  Display components based on the estimated number of clusters.
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
  \code{\link{dpmixsim}},
  \code{\link{readsliceimg}},
  \code{\link{premask}},
  \code{\link{postdpmixciz}},
  \code{\link{postdataseg}},
  \code{\link{postimgclgrp}}
}

\examples{
\dontrun{
## see Example 2 in `dpmixsim'.
}}

\keyword{dplot}
