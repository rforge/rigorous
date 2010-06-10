\name{dpmixsim}
\alias{dpmixsim}
\concept{bayes}
\concept{MCMC}
\concept{Gibbs Sampling}
\concept{Dirichlet Process Mixture models}

\title{ Dirichlet Process mixture model for clustering and image segmentation }
\description{
  \code{dpmixsim} implements a Dirichlet Process mixture (DPM) model.
  The DPM model is a Bayesian nonparametric methodology that relies on MCMC simulations
  for exploring mixture models with an unknown number of components.
  The function implements conjugate models with normal structure (conjugate normal-normal
  DP mixture model).

}

\usage{
dpmixsim(x, M=1, a=1, b=2, upalpha=1, a0=2, b0=2, maxiter=4000, rec=1000,
  fsave=NA, kmax=30)
}
\arguments{
  \item{x}{ scaled input data as vector in range \{0,1\} }
  \item{M}{ DP precision hyperparameter }
  \item{a}{ Gamma prior hyperparameter }
  \item{b}{ Gamma prior hyperparameter }
  \item{upalpha}{ is a logical variable (default = \samp{TRUE}) for
  simulations with \{automatic,fixed\} calibration of the precision
  hyperparameter M }
  \item{a0}{ Gamma prior hyperparameter for M }
  \item{b0}{ Gamma prior hyperparameter for M } 
  \item{maxiter}{ maximum number of MCMC iteration steps }
  \item{rec}{ record the last \samp{rec} iteration steps }
  \item{fsave}{ filename for saving the MCMC simulation (def: \samp{NULL} do not save) }
  \item{kmax}{ maximum number of clusters in the simulation, default 30 }
}

\details{
  Consider n observations \eqn{x_1},...,\eqn{x_n} which we regard as
  exchangeable.  We model the distribution from which the \eqn{x_i} are
  drawn as a mixture of distributions. Dirichlet process mixture models
  are based on Dirichlet process priors for the primary parameters \eqn{\theta_i}.
  DP mixture models assume that the prior distribution function \eqn{G}
  itself is uncertain, drawn from a Dirichlet process
  \eqn{G} \eqn{\sim}{~} \eqn{DP(M G_0)},
  with base prior \eqn{G_0} and precision parameter \eqn{M}.
  This specification may be expressed by the hierarchical model:

  \eqn{x_i} \eqn{\sim}{~} \eqn{N(.|\theta_i, \sigma^2)} \cr
  \eqn{\theta_i} \eqn{\sim}{~} G  \cr
  \eqn{G} \eqn{\sim}{~} \eqn{DP(M  N(0,1))}  \cr
  \eqn{\sigma^{-2}} \eqn{\sim}{~} \eqn{Gamma(a,b)} 
}

\value{
  simulation output as a list of draws containing:
  \item{krec}{ cluster indicator variables }
  \item{wrec}{ cluster weights }
  \item{phirec}{ theta cluster parameters }
  \item{varrec}{ sigma cluster parameters }
}


\references{ 
Adelino Ferreira da Silva, A Dirichlet process mixture model for
 brain MRI tissue classification,
 \emph{Medical Image Analysis} 11 (2007) 169-182.

Adelino Ferreira da Silva, Bayesian mixture models of variable
 dimension for image segmentation,
 \emph{Comput. Methods Programs Biomed.} 94 (2009) 1-14.
}

\author{ Adelino Ferreira da Silva, Universidade Nova de Lisboa,
  Faculdade de Ciencias e Tecnologia, Portugal, \email{afs@fct.unl.pt}.
}

\seealso{
  \code{\link{readsliceimg}},
  \code{\link{postdataseg}},
  \code{\link{postdpmixciz}},
  \code{\link{postimgclgrp}},
  \code{\link{postimgcomps}},
  \code{\link{postkcluster}},
  \code{\link{premask}},
  \code{\link{readsliceimg}}
}

\examples{
## Example 1: simple test using 'galaxy' data
   data("galaxy")
   x0 <- galaxy$speed
   x <- prescale(x0) 
   res <- dpmixsim(x, M=1, a=1, b=0.001, upalpha=1, maxiter=4000, rec=3000)
## post-simulation
   z <- postdpmixciz(x, res=res, rec=3000, ngrid=100, plot=TRUE)
\dontrun{
## Example 2: MRI segmentation
## Testing note: this example should reproduce the equivalent segmented images
## used in the author's references 
  slicedata <- readsliceimg(fbase="t1_pn3_rf0", swap=FALSE)
  image(slicedata$niislice, col=gray((0:255)/256), main="original image")
  x0 <- premask(slicedata, subsamp=TRUE)
  x  <- prescale(x0) 
  fsave <- "/tmp/segtest2.sav"
  answ <- readline("Run MCMC simulation first ? ")
  run <- FALSE
  if (substr(answ, 1, 1) == "y") { run <- TRUE }
  rec <- 2000
  ngrid <- 100
  if(run) {
    res <- dpmixsim(x, M=1, a=1, b=2, upalpha=1, maxiter=4000, rec=rec,
      fsave=fsave)
  } else {
    load(file=fsave)
    cat("loaded",fsave,"\n")
  }
  ## post-simulation
  z <- postdpmixciz(x, res=res, rec=rec, ngrid=ngrid, plot=TRUE)
  x0 <- premask(slicedata, subsamp=FALSE) # use full-sized image after estimation 
  x  <- prescale(x0) 
  cx   <- postdataseg(x, z, ngrid=ngrid)
  cat("*** view grouped segmentations:\n")
  postimgclgrp(slicedata$mask, cx, palcolor=FALSE)
  cat("*** display all clusters:\n")
  postimgcomps(slicedata$mask, cx)
  cat("*** re-cluster with 4 clusters:\n")
  postkcluster(slicedata$mask, cx, clk=4)
}
}

\keyword{models}
