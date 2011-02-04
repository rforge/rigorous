\name{aslModels}
\alias{asl2p}
\alias{asl3p}
\title{Kinetic Models for Arterial Spin Labeling}
\description{
  Theoretical models for the acquisition of pulsed ASL data.
}
\usage{
asl2p(beta, TI, T1b = 1.3, T1 = 1.0, tau = 0.7, alpha = 0.9,
      lambda = 0.9, Mob = 4095)
asl3p(beta, TI, T1b = 1.3, T1 = 1.0, alpha = 0.9, lambda = 0.9,
      Mob = 100)
}
\arguments{
  \item{beta}{is a vector (lenth = 2 or 3) of kinetic parameters.}
  \item{TI}{is a vector of inversion times (TIs) in seconds.}
  \item{T1b}{is the longitudinal relaxation time of blood.}
  \item{T1}{is the longitudinal relaxation time of tissue.}
  \item{tau}{is the first TI in the second paper (?).}
  \item{alpha}{is the fraction of maximum possible change in the
    longitudinal magnetization that was achieved.}
  \item{lambda}{is the equilibrium tissue/blood partition coefficient of
    water.} 
  \item{Mob}{is the equilibrium magnetization of arterial blood.}
}
\details{
  See Buxton etal. (1998) for more details.
}
\value{
  A theoretical curve of pulsed ASL signal as a function of time.
}
\references{
  Buxton, R.B., Frank, L.R., Wong, E.C., Siewert, B., Warach, S. and
  Edelman, R.R. (1998) A General Kinetic Model for Quantitative
  Perfusion Imaging with Arterial Spin Labeling, \emph{Magnetic
    Resonance in Medicine}, \bold{40}, 383-396.
}
%\seealso{}
\examples{
TI <- seq(0, 3, length=256)
beta <- c(f=0.8, deltaT=0.5)
plot(TI, asl2p(beta, TI), type="l", lwd=2,
     xlab="Time (sec)", ylab="Signal", main="Pulsed ASL: Standard Model")
abline(v=1.5, lwd=2, col=2)
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{models}
