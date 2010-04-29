\name{postimgclgrp}
\alias{postimgclgrp}
\concept{Display segmented image}

\title{Segment image with the estimated number of  components}

\description{
  \code{postimgclgrp} displays the segmented image with the estimated number
  of components
}

\usage{postimgclgrp(mask, cx, palcolor=TRUE)}

\arguments{
  \item{mask}{ full-sized scaled image data prepared by \code{premask}}
  \item{cx}{ data segmentation prepared by \code{postdataseg} }
  \item{palcolor}{ logical variable for selecting colored/grey image
  visualization (default = \samp{TRUE}) }
}

\details{
  Display image segmentation with the estimated number of components.
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
  \code{\link{postdataseg}}
}
\examples{
\dontrun{
## see Example 2 in `dpmixsim'.
}}
\keyword{dplot}
