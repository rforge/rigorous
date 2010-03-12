% Part of the oro.dti package for R
% Distributed under the BSD license: see oro.dti/COPYING
%
% $Id: $

\name{Diffusion Tensor Imaging}
\alias{tensor.slice}
\alias{tensor.volume}
\alias{tensor.slice.multiple}
\title{Parameter Estimation for the Diffusion Tensor Model}
\description{
  ...
}
\usage{
tensor.slice(slice, G, b, weight=FALSE, mask=NULL, boot=FALSE, Nboot=99,
             HC=2, verbose=TRUE, delta=0.01)
tensor.slice.multiple(slice, G, b, weight=TRUE, scan.seq=NULL,
                      mask=NULL, boot=FALSE, Nboot=99, verbose=TRUE)
tensor.volume(volume, G, b, mask, weight=FALSE, boot=FALSE, Nboot=99,
              HC=2, verbose=TRUE)
}
\arguments{
  \item{slice}{}
  \item{volume}{}
  \item{G}{}
  \item{b}{}
  \item{weight}{}
  \item{mask}{}
  \item{boot}{}
  \item{Nboot}{}
  \item{HC}{}
  \item{verbose}{}
  \item{delta}{}
  \item{scan.seq}{}
}
\details{
}
\value{
}
%\references{}
%\seealso{}
%\examples{}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\keyword{}
