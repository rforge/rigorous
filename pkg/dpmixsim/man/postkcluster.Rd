\name{postkcluster}
\alias{postkcluster}
\concept{Re-cluster and display image with a reduced number components}

\title{Segmentation with a fixed number of clusters} 
\description{
  \code{postkcluster} re-clusters the data with a user-specified number
  of components, and displays the segmented image.
}

\usage{postkcluster(mask, cx, clk=4, plot=TRUE)}
\arguments{
  \item{mask}{ masked full-sized image data prepared by \code{premask}}
  \item{cx}{ data segmentation prepared by \code{postdataseg} }
  \item{clk}{ desired fixed number of components, including the background
  component, to use in the data segmentation; default \samp{clk=4}:
  gray matter (GM), white matter (WM), CSF, and background} 
  \item{plot}{ logical variable; enables suspension of output images
  (default = \samp{TRUE}) }  
}
\details{
  Partitioning clustering around medoids (PAM) is applied to the classes
  simulated from \code{dpmixsim} as a post-processing step.
  This procedure may be applied to merge clusters, and reduce the number
  of clusters to the specified value \samp{clk}.
  \code{postkcluster} computes a \code{clara} object using \pkg{cluster}
  (see Struyf et.al.), a list representing a clustering of the data into
  \samp{clk} clusters.
}

\references{ 
Adelino Ferreira da Silva, A Dirichlet process mixture model for
 brain MRI tissue classification,
 \emph{Medical Image Analysis} 11 (2007) 169-182.

Adelino Ferreira da Silva, Bayesian mixture models of variable
 dimension for image segmentation,
 \emph{Comput. Methods Programs Biomed.} 94 (2009) 1-14.
}

\references{ 
Anja Struyf, Mia Hubert & Peter J. Rousseeuw (1996): Clustering in an
  Object-Oriented Environment. Journal of Statistical Software, 1.
  http://www.stat.ucla.edu/journals/jss/
}

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs@fct.unl.pt}.
}

\seealso{
  \code{\link{dpmixsim}},
  \code{\link{readsliceimg}},
  \code{\link{premask}},
  \code{\link{postdpmixciz}},
  \code{\link{postdataseg}},
  \code{\link{postimgcomps}}
}

\examples{
\dontrun{
## see Example 2 in `dpmixsim'.
}}

\keyword{cluster}
