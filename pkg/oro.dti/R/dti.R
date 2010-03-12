##
## Copyright (c) 2010, Brandon Whitcher
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * Neither the name of Rigorous Analytics Ltd. nor the names of
##       its contributors may be used to endorse or promote products 
##       derived from this software without specific prior written 
##       permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
## $Id: $
##

tensor.slice <- function(slice, G, b, weight=FALSE, mask=NULL,
                         boot=FALSE, Nboot=99, HC=2, verbose=TRUE,
                         delta=0.01) {
  ##
  ## Estimates diffusion tensor (D) from data provided by D. Tuch
  ##
  M <- nrow(slice)
  N <- ncol(slice)
  p <- 7
  ## Construct design matrix
  Z <- cbind(-b*cbind(G$x^2, G$y^2, G$z^2, 2*G$x*G$y, 2*G$x*G$z, 2*G$y*G$z),
             rep(1, nrow(G))) 
  ## Format response into a matrix
  ## Y <- matrix(unlist(log(slice + delta)), nrow(G), byrow=TRUE)
  Y <- matrix(unlist(log(slice+delta)), nrow(G), byrow=TRUE)
  ## Estimation via multivariate multiple linear regression
  beta <- qr.coef(qr(Z), Y)
  ## Hat matrix
  H <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  if(weight) {
    ## Compute weight matrix
    V <- exp(H %*% Y)
    if(!is.null(mask))
      V <- V[, as.logical(unlist(mask))]
    else
      stop("Cannot estimate weighting matrix without mask")
    V <- rowMeans(V, na.rm=TRUE)
    ## Weighted regression
    wbeta <- qr.coef(qr(V * Z), V * Y)
    beta <- wbeta
  }
  ## Fitted values
  Yhat <- Z %*% beta
  ## Residuals
  u <- Y - Yhat
  if(boot) {
    ## Distribution from Mammen (1993)
    F1 <- function(n) {
      p <- (sqrt(5)+1)/2/sqrt(5)
      sample(c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2), n, replace=TRUE, prob=c(p, 1-p))
    }
    ## Distribution from Davidson and Flachaire (2001)
    F2 <- function(n) {
      sample(c(-1,1), n, replace=TRUE)
    }
    ## Prepare parameters for the residuals
    m <- nrow(Y)
    n <- ncol(Y)
    k <- 7
    a <- switch(HC, rep(sqrt(m/(m-k)), m), 1/sqrt(1-diag(H)), 1/(1-diag(H)))
    ## Perform the wild bootstrap
    beta.boot <- array(NA, c(k,n,Nboot+1))
    beta.boot[,,1] <- beta
    for(i in 1:Nboot+1) {
      if(!(i %% 10)) cat("  i =", i, fill=TRUE)
      ################################
      epsilon <- matrix(F2(m*n), m, n)
      ################################
      estar <- a * u * epsilon
      Ystar <- Yhat + estar
      beta.boot[,,i] <- qr.coef(qr(Z), Ystar)
    }
  } else beta.boot <- NULL
  ## Residual sum of squares (RSS)
  ## RSS <- t(Y) %*% Y - t(Yhat) %*% Yhat
  diag.RSS <- colSums(u^2) / (nrow(G) - 7)
  beta.var <- array(NA, c(M,N,7))
  for(i in 1:7)
    beta.var[,,1] <- matrix(diag.RSS, M, N) * solve(t(Z) %*% Z)[i,i]
  ## Cook's Distance ## 
  hi <- matrix(diag(H), nrow(G), M*N)
  sigma.hat <- sqrt(matrix(diag.RSS, nrow(G), M*N, byrow=TRUE))
  ## Studentized residuals...
  r <- u / (sigma.hat * sqrt(1 - hi))
  CD <- r^2 / p * hi / (1 - hi)
  ## Diffusion tensor (D) is stored as an array.  The third dimension
  ## contains the estimated tensor values.
  DT <- array(c(t(beta)), c(M,N,p))
  CD <- array(c(t(CD)), c(M,N,nrow(G)))
  list(tensor = DT, se = sqrt(beta.var), boot = beta.boot, CooksD = CD)
}

tensor.volume <- function(volume, G, b, mask, weight=FALSE, boot=FALSE,
                          Nboot=99, HC=2, verbose=TRUE) {
  M <- nrow(volume)
  N <- ncol(volume)
  P <- nsli(volume)
  p <- 7
  delta <- 0.01
  ## Construct design matrix
  Z <- cbind(-b*cbind(G$x^2, G$y^2, G$z^2, 2*G$x*G$y, 2*G$x*G$z, 2*G$y*G$z),
             rep(1, nrow(G))) 
  ## Hat matrix
  H <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  ## Format response into a matrix from all unmasked voxels
  mask4D <- array(as.logical(mask), dim=dim(volume))
  vec <- volume[mask4D]
  Y <- matrix(log(vec + delta), nrow(G), byrow=TRUE)
  ## Estimation via multivariate multiple linear regression
  beta <- qr.coef(qr(Z), Y)
  ## Diffusion tensor (D) is stored as an array.  The third dimension
  ## contains the estimated tensor values.
  DT <- array(0, dim=c(M,N,P,p))
  DT[mask4D[,,,1:p]] <- unlist(t(beta))

  ## Fitted values
  Yhat <- Z %*% beta
  ## Residuals
  u <- Y - Yhat
  if(boot) {
    ## Distribution from Mammen (1993)
    F1 <- function(n) {
      p <- (sqrt(5) + 1)/2/sqrt(5)
      sample(c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2), n, replace=TRUE, prob=c(p, 1-p))
    }
    ## Distribution from Davidson and Flachaire (2001)
    F2 <- function(n) {
      sample(c(-1,1), n, replace=TRUE)
    }
    ## Prepare parameters for the residuals
    m <- nrow(Y)
    n <- ncol(Y)
    k <- 7
    a <- switch(HC, rep(sqrt(m/(m-k)), m), 1/sqrt(1-diag(H)), 1/(1-diag(H)))
    ## Perform the wild bootstrap
    beta.boot <- array(NA, c(k,n,Nboot+1))
    beta.boot[,,1] <- beta
    for(i in 1:Nboot+1) {
      if(!(i %% 10)) cat("  i =", i, fill=TRUE)
      epsilon <- matrix(F2(m*n), m, n)
      estar <- a * u * epsilon
      Ystar <- Yhat + estar
      beta.boot[,,i] <- qr.coef(qr(Z), Ystar)
    }
  } else beta.boot <- NULL
  list(D = DT, b = beta.boot)
}

tensor.slice.multiple <- function(slice, G, b, weight=TRUE, scan.seq=NULL,
                                  mask=NULL, boot=FALSE, Nboot=99,
                                  verbose=TRUE) {
  ##
  ## Estimates diffusion tensor (D) from data provided by D. Tuch
  ##
  M <- nrow(slice)
  N <- ncol(slice)
  delta <- .001
  ## Construct design matrix
  Z <- cbind(-b*cbind(G$x^2, G$y^2, G$z^2, 2*G$x*G$y, 2*G$x*G$z, 2*G$y*G$z),
             rep(1, nrow(G)))
  p <- ncol(Z)
  ## Format response into a matrix
  Y <- matrix(unlist(log(slice + delta)), nrow(G), byrow=TRUE)
  ## Estimation via multivariate multiple linear regression
  beta <- qr.coef(qr(Z), Y)
  ## Hat matrix
  ## H <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  if(weight) {
    if(is.null(scan.seq))
      stop("Cannot estimate weighting matrix without scanning sequence")
    ## Compute weight matrix
    V7 <- Y[,as.logical(unlist(mask))]
    V <- numeric(nrow(Y))
    for(i in unique(scan.seq)) {
      index <- which(i == scan.seq)
      V[index] <- mad(V7[index,])
    }
    ## Weighted regression
    qrVZ <- qr((1/V) * Z)
    wbeta <- qr.coef(qrVZ, (1/V) * Y)
    beta <- wbeta
  }
  ## Fitted values
  Yhat <- Z %*% beta
  ## Residuals
  u <- Y - Yhat
  if(boot) {
    ## Perform the non-parametric bootstrap
    beta.boot <- array(NA, c(nrow(beta),ncol(beta),Nboot+1))
    beta.boot[,,1] <- beta
    for(j in 2:(Nboot+1)) {
      Yb <- matrix(NA, nrow(Y), ncol(Y))
      ## Resample the b=0 images
      s0 <- sample(which(scan.seq == 0), replace=TRUE)
      Yb[scan.seq == 0,] <- Y[s0,]
      ## Resample the chunks of diffusion-weighted images
      s1 <- sample(1:10, replace=TRUE)
      for(k in 1:10) {
        ind <- (7 * (k-1)):(7 * k - 2) + 2
        inds <- (7 * (s1[k] - 1)):(7 * s1[k] - 2) + 2
        Yb[ind,] <- Y[inds,] 
      }
      ## Weighted regression
      beta.boot[,,j] <- qr.coef(qrVZ, (1/V) * Yb)
    }
  }
  else
    beta.boot <- NULL
  ## Residual sum of squares (RSS)
  ## RSS <- t(Y) %*% Y - t(Yhat) %*% Yhat
  ## diag.RSS <- colSums(u^2) / (nrow(G) - p)
  beta.var <- array(NA, c(M,N,p))
  DT <- array(c(t(beta)), c(M,N,p))
  list(tensor = DT, se = sqrt(beta.var), boot=beta.boot)
}

eigen.slice <- function(DT, analytic=TRUE, debug=FALSE) {
  ##
  ## Computes the eigenvalue decomposition for the diffusion tensor
  ##

  M <- nrow(DT)
  N <- ncol(DT)

  if(analytic) {
    ## Analytical Computation of the Eigenvalues and Eigenvectors in DT-MRI
    ## Hasan, Basser, Parker and Alexander (2001)
    ## Journal of Magnetic Resonance 152, 41-47.
    
    ## Invariants...
    I1 <- DT[,,1] + DT[,,2] + DT[,,3]
    I2 <- DT[,,1]*DT[,,2] + DT[,,1]*DT[,,3] + DT[,,2]*DT[,,3] -
      (DT[,,4]^2 + DT[,,5]^2 + DT[,,6]^2)
    I3 <- DT[,,1]*DT[,,2]*DT[,,3] + 2*DT[,,4]*DT[,,5]*DT[,,6] -
      (DT[,,3]*DT[,,4]^2 + DT[,,2]*DT[,,5]^2 + DT[,,1]*DT[,,6]^2)
    ## Determination of the eigenvalues
    v <- (I1/3)^2 - I2/3
    s <- (I1/3)^3 - I1*I2/6 + I3/2
    phi <- acos(s/v/sqrt(v))/3
    lambda <- array(NA, c(M,N,3))
    lambda[,,1] <- I1/3 + 2*sqrt(v) * cos(phi)
    lambda[,,2] <- I1/3 - 2*sqrt(v) * cos(pi/3 + phi)
    lambda[,,3] <- I1/3 - 2*sqrt(v) * cos(pi/3 - phi)
    ## Determination of the eigenvectors
    A <- array(DT[,,1], c(M,N,3)) - lambda
    B <- array(DT[,,2], c(M,N,3)) - lambda
    C <- array(DT[,,3], c(M,N,3)) - lambda
    ex <- ey <- ez <- array(NA, c(M,N,3))
    for(i in 1:3) {
      ex[,,i] <- (DT[,,4]*DT[,,6] - B[,,i]*DT[,,5]) *
        (DT[,,5]*DT[,,6] - C[,,i]*DT[,,4])
      ey[,,i] <- (DT[,,5]*DT[,,6] - C[,,i]*DT[,,4]) *
        (DT[,,5]*DT[,,4] - A[,,i]*DT[,,6])
      ez[,,i] <- (DT[,,4]*DT[,,6] - B[,,i]*DT[,,5]) *
        (DT[,,5]*DT[,,4] - A[,,i]*DT[,,6])
    }
    e <- vector("list", 3)
    for(i in 1:3) {
      e[[i]] <- array(NA, c(M,N,3))
      denom <- sqrt(ex[,,i]^2 + ey[,,i]^2 + ez[,,i]^2)
      e[[i]][,,1] <- ex[,,i] / denom
      e[[i]][,,2] <- ey[,,i] / denom
      e[[i]][,,3] <- ez[,,i] / denom
    }
    return(list(value=lambda, vector=e))
  } else {
    ## Function that creates a 3x3 matrix for the diffusion tensor
    sym.mat <- function(x)
      matrix(c(x[1], x[4], x[5], x[4], x[2], x[6], x[5], x[6], x[3]), 3, 3)
    
    DT.value <- array(NA, c(M,N,3))
    DT.vector <- vector("list", 3)
    for(i in 1:3)
      DT.vector[[i]] <- array(NA, c(M,N,3))
    for(i in 1:M) {
      if(debug) cat(" i =", i, fill=TRUE)
      for(j in 1:N) {
        if(any(is.na(DT[i,j,]))) {
          DT.vector[[1]][i,j,] <- DT.vector[[2]][i,j,] <-
            DT.vector[[3]][i,j,] <- rep(NA,3)
        } else {
          temp <- eigen(sym.mat(DT[i,j,]), symmetric=TRUE)
          DT.value[i,j,] <- temp$values
          if(all(is.na(temp$values))) {
            DT.vector[[1]][i,j,] <- DT.vector[[2]][i,j,] <-
              DT.vector[[3]][i,j,] <- rep(NA,3)
          } else {
            DT.vector[[1]][i,j,] <- as.vector(temp$vectors[,1])
            DT.vector[[2]][i,j,] <- temp$vectors[,2]
            DT.vector[[3]][i,j,] <- temp$vectors[,3]
          }
        }
      }
    }
    return(list(value=DT.value, vector=DT.vector))
  }
}

eigen.volume <- function(DT, debug=FALSE) {
  ##
  ## Computes the eigenvalue decomposition for the diffusion tensor
  ##
  M <- nrow(DT)
  N <- ncol(DT)
  P <- dim(DT)[3]
  
  ## Analytical Computation of the Eigenvalues and Eigenvectors in DT-MRI
  ## Hasan, Basser, Parker and Alexander (2001)
  ## Journal of Magnetic Resonance 152, 41-47.
    
  ## Invariants...
  I1 <- DT[,,,1] + DT[,,,2] + DT[,,,3]
  I2 <- DT[,,,1]*DT[,,,2] + DT[,,,1]*DT[,,,3] + DT[,,,2]*DT[,,,3] -
    (DT[,,,4]^2 + DT[,,,5]^2 + DT[,,,6]^2)
  I3 <- DT[,,,1]*DT[,,,2]*DT[,,,3] + 2*DT[,,,4]*DT[,,,5]*DT[,,,6] -
    (DT[,,,3]*DT[,,,4]^2 + DT[,,,2]*DT[,,,5]^2 + DT[,,,1]*DT[,,,6]^2)
  ## Determination of the eigenvalues
  v <- (I1/3)^2 - I2/3
  s <- (I1/3)^3 - I1*I2/6 + I3/2
  phi <- acos(s/v/sqrt(v))/3
  lambda <- array(NA, c(M,N,P,3))
  lambda[,,,1] <- I1/3 + 2*sqrt(v) * cos(phi)
  lambda[,,,2] <- I1/3 - 2*sqrt(v) * cos(pi/3 + phi)
  lambda[,,,3] <- I1/3 - 2*sqrt(v) * cos(pi/3 - phi)
  ## Determination of the eigenvectors
  A <- array(DT[,,,1], c(M,N,P,3)) - lambda
  B <- array(DT[,,,2], c(M,N,P,3)) - lambda
  C <- array(DT[,,,3], c(M,N,P,3)) - lambda
  ex <- ey <- ez <- array(NA, c(M,N,P,3))
  for(i in 1:3) {
    ex[,,,i] <- (DT[,,,4]*DT[,,,6] - B[,,,i]*DT[,,,5]) *
      (DT[,,,5]*DT[,,,6] - C[,,,i]*DT[,,,4])
    ey[,,,i] <- (DT[,,,5]*DT[,,,6] - C[,,,i]*DT[,,,4]) *
      (DT[,,,5]*DT[,,,4] - A[,,,i]*DT[,,,6])
    ez[,,,i] <- (DT[,,,4]*DT[,,,6] - B[,,,i]*DT[,,,5]) *
      (DT[,,,5]*DT[,,,4] - A[,,,i]*DT[,,,6])
  }
  e <- vector("list", 3)
  for(i in 1:3) {
    e[[i]] <- array(NA, c(M,N,P,3))
    denom <- sqrt(ex[,,,i]^2 + ey[,,,i]^2 + ez[,,,i]^2)
    e[[i]][,,,1] <- ex[,,,i] / denom
    e[[i]][,,,2] <- ey[,,,i] / denom
    e[[i]][,,,3] <- ez[,,,i] / denom
  }
  return(list(values=lambda, vectors=e))
}

eigen.slice.boot <- function(mat, debug=FALSE) {
  ## Analytical Computation of the Eigenvalues and Eigenvectors in DT-MRI
  ## Hasan, Basser, Parker and Alexander (2001)
  ## Journal of Magnetic Resonance 152, 41-47.

  ## Invariants...
  I1 <- mat[1,,] + mat[2,,] + mat[3,,]
  I2 <- mat[1,,]*mat[2,,] + mat[1,,]*mat[3,,] + mat[2,,]*mat[3,,] -
    (mat[4,,]^2 + mat[5,,]^2 + mat[6,,]^2)
  I3 <- mat[1,,]*mat[2,,]*mat[3,,] + 2*mat[4,,]*mat[5,,]*mat[6,,] -
    (mat[3,,]*mat[4,,]^2 + mat[2,,]*mat[5,,]^2 + mat[1,,]*mat[6,,]^2)
  ## Determination of the eigenvalues
  v <- (I1/3)^2 - I2/3
  s <- (I1/3)^3 - I1*I2/6 + I3/2
  phi <- acos(s/v/sqrt(v))/3
  dim1 <- c(3,dim(mat)[2],dim(mat)[3])
  lambda <- array(NA, dim1)
  lambda[1,,] <- I1/3 + 2*sqrt(v) * cos(phi)
  lambda[2,,] <- I1/3 - 2*sqrt(v) * cos(pi/3 + phi)
  lambda[3,,] <- I1/3 - 2*sqrt(v) * cos(pi/3 - phi)
  ## Determination of the eigenvectors
  dim2 <- c(dim(mat)[2],dim(mat)[3],3)
  A <- aperm(array(mat[1,,], dim2), c(3,1,2)) - lambda
  B <- aperm(array(mat[2,,], dim2), c(3,1,2)) - lambda
  C <- aperm(array(mat[3,,], dim2), c(3,1,2)) - lambda
  ex <- ey <- ez <- array(NA, dim1)
  for(i in 1:3) {
    ex[i,,] <- (mat[4,,]*mat[6,,] - B[i,,]*mat[5,,]) *
      (mat[5,,]*mat[6,,] - C[i,,]*mat[4,,]) # D_{xz} or D_{xy}
    ey[i,,] <- (mat[5,,]*mat[6,,] - C[i,,]*mat[4,,]) *
      (mat[5,,]*mat[4,,] - A[i,,]*mat[6,,])
    ez[i,,] <- (mat[4,,]*mat[6,,] - B[i,,]*mat[5,,]) *
      (mat[5,,]*mat[4,,] - A[i,,]*mat[6,,])
  }
  e <- vector("list", 3)
  for(i in 1:3) {
    e[[i]] <- array(NA, dim1)
    denom <- sqrt(ex[i,,]^2 + ey[i,,]^2 + ez[i,,]^2)
    e[[i]][1,,] <- -ex[i,,] / denom
    e[[i]][2,,] <- -ey[i,,] / denom
    e[[i]][3,,] <- -ez[i,,] / denom
    }
  list(value=lambda, vector=e)
}

fractional.anisotropy <- function(x) {
  ## Mean diffusivity
  MD <- rowMeans(x, dims=2)
  ## numerator
  n <- sqrt((x[,,1] - MD)^2 + (x[,,2] - MD)^2 + (x[,,3] - MD)^2)
  ## denominator
  d <- sqrt(x[,,1]^2 + x[,,2]^2 + x[,,3]^2)
  return(sqrt(3/2) * n/d)     # According to D. Tuch
}

FA.volume <- function(x) {
  ## Mean diffusivity
  MD <- rowMeans(x, dims=3)
  ## numerator
  n <- sqrt((x[,,,1] - MD)^2 + (x[,,,2] - MD)^2 + (x[,,,3] - MD)^2)
  ## denominator
  d <- sqrt(x[,,,1]^2 + x[,,,2]^2 + x[,,,3]^2)
  return(sqrt(3/2) * n/d)     # According to D. Tuch
}

fa.boot <- function(x) {
  ## Mean diffusivity
  MD <- colMeans(x, dims=1)
  ## numerator
  n <- sqrt((x[1,,] - MD)^2 + (x[2,,] - MD)^2 + (x[3,,] - MD)^2)
  ## denominator
  d <- sqrt(x[1,,]^2 + x[2,,]^2 + x[3,,]^2)
  return(sqrt(3/2) * n/d)     # According to D. Tuch
}

FA <- function(x) {
  ## Mean diffusivity
  MD <- mean(x)
  ## numerator
  n <- sqrt((x[1] - MD)^2 + (x[2] - MD)^2 + (x[3] - MD)^2)
  ## denominator
  d <- sqrt(x[1]^2 + x[2]^2 + x[3]^2)
  return(sqrt(3/2) * n/d)     # According to D. Tuch
}

boot.wild <- function(G, y, b, HC=2, Nboot=99, weight=FALSE) {
  ## Distribution from Mammen (1993)
  F1 <- function(n) {
    p <- (sqrt(5)+1)/2/sqrt(5)
    sample(c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2), n, replace=TRUE, prob=c(p, 1-p))
  }
  ## Distribution from Davidson and Flachaire (2001)
  F2 <- function(n)
    sample(c(-1,1), n, replace=TRUE)

  ## Perform initial multiple linear regression
  X <- cbind(-b*cbind(G$x^2, G$y^2, G$z^2, 2*G$x*G$y, 2*G$x*G$z, 2*G$y*G$z),
             rep(1, nrow(G)))
  beta <- qr.coef(qr(X), log(y))
  ## Hat matrix
  h <- X %*% solve(t(X) %*% X) %*% t(X)
  if(weight) {
    ## Compute weight matrix
    V <- as.vector(exp(h %*% log(y)))
    ## Weighted regression
    wbeta <- qr.coef(qr(V * X), V * log(y))
    beta <- wbeta
  }
  yhat <- as.vector(X %*% beta)
  u <- log(y) - yhat
  k <- 7
  s2 <- t(u) %*% u / (nrow(G) - k)
  
  ## Prepare parameters for the residuals
  n <- length(y)
  a <- switch(HC, rep(sqrt(n/(n-k)), n), 1/sqrt(1-diag(h)), 1/(1-diag(h)))

  ## Perform the wild bootstrap
  beta.boot <- matrix(NA, Nboot+1, k)
  beta.boot[1,] <- beta
  for(i in 1:Nboot + 1) {
    epsilon <- F2(n)
    estar <- a * u * epsilon
    ystar <- yhat + estar
    if(weight)
      beta.boot[i,] <- qr.coef(qr(V * X), V * ystar)
    else
      beta.boot[i,] <- qr.coef(qr(X), ystar)
  }
  list(boot = beta.boot, stan = sqrt(s2 * diag(solve(t(X)%*%X))))
}

logm <- function(A) {
  S <- svd(A)
  S$v %*% log(diag(S$d)) %*% solve(S$v)
}

expm <- function(A) {
  S <- svd(A)
  S$v %*% exp(diag(S$d)) %*% solve(S$v)
}
