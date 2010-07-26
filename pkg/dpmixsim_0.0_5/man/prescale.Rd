\name{prescale}
\alias{prescale}
\concept{data preparation}
\title{ Data preparation }
\description{
  \code{prescale} scales data to be in the range \{0,1\}, as a preparation 
  for simulation.
}
\usage{
  prescale(xv)
}
\arguments{
  \item{xv}{ unscaled data vector }
}
\value{
  \item{x}{ scaled data vector }
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
  data("galaxy")
  x0 <- galaxy$speed
  x  <- prescale(x0) 
  print(range(x))
}}

\keyword{utilities}
