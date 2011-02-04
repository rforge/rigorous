asl2p <- function(beta, TI, T1b=1.3, T1=1.0, tau=1.0, alpha=0.9,
                  lambda=0.9, Mob=4095) {
  
  ## explain: this function is actually PASL function.

  if (length(beta) != 2) {
    stop("There must be two parameters in the beta vector.")
  }
  f <- beta[1]
  deltaT <- beta[2]

  epoch1 <- TI <= deltaT
  epoch2 <- TI > deltaT & TI < tau + deltaT
  epoch3 <- TI >= tau + deltaT
  
  k <- 1/T1b - 1/T1 - f/lambda
  
  deltaM <- rep(0, length(TI))
  ## fun[epoch1] <- 0
  deltaM[epoch2] <- 2 * Mob * f * (TI[epoch2] - deltaT) * alpha *
    exp(-TI[epoch2] / T1b) * exp(k * TI[epoch2]) *
      (exp(-k * deltaT) - exp(-k * TI[epoch2])) / (k * (TI[epoch2] - deltaT))
  deltaM[epoch3] <- 2 * Mob * f * tau * alpha * exp(-TI[epoch3] / T1b) *
    exp(k * TI[epoch3]) * (exp(-k * deltaT) - exp(-k * (tau + deltaT))) /
      (k * tau) 

  return(deltaM)
  
}

asl3p <- function(beta, TI, T1b=1.3, T1=1.0, alpha=0.9, lambda=0.9,
                  Mob=100) {

  ## explain: this function is actually PASL function.

  if (length(beta) != 3) {
    stop("There must be three parameters in the beta vector.")
  }
  f <- beta[1]
  deltaT <- beta[2]
  tau <- beta[3]

  ## TI=linspace(0,3.5,701)
  epoch1 <- TI <= deltaT
  epoch2 <- TI > deltaT & TI < (tau + deltaT)
  epoch3 <- TI >= (tau + deltaT)

  k <- 1/T1b - 1/T1 - f/lambda

  fun <- rep(0, length(TI))
  ## fun(epoch1) = 0
  fun[epoch2] <- 2 * Mob * f * (TI[epoch2] - deltaT) * alpha *
    exp(-TI[epoch2] / T1b) * exp(k * TI[epoch2]) *
      (exp(-k * deltaT) - exp(-k * TI[epoch2])) / (k * (TI[epoch2] - deltaT)) 
  fun[epoch3] <- 2 * Mob * f * tau * alpha * exp(-TI[epoch3] / T1b) *
    exp(k * TI[epoch3]) * (exp(-k * deltaT) - exp(-k * (tau + deltaT))) /
      (k * tau)

  return(fun)

}
