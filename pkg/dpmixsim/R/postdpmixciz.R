##
## post-analysis after dpmixsim
##
postdpmixciz <- function(x, res, kmax=30, rec=300, ngrid=200, plot=TRUE) {
  ## relabel to keep values in sequence
  relabel <- function(z) {
    u <- sort(unique(z))
    v <- 1:length(u)
    for (k in 1:length(u)) {
      w <- which(z == u[k])
      z[w] <- v[k]
    }
    invisible(z)
  }
  krec <- res$krec
  wrec <- matrix(res$w, nc=kmax, byrow=TRUE)
  phirec <- matrix(res$phirec, nc=kmax, byrow=TRUE)
  varrec <- res$varrec
  ## Histogram for k
  khist <- numeric(kmax)
  for (i in 1:length(krec)) {
    khist[krec[i]] <- khist[krec[i]] + 1  
  }
  kfreq <- which(khist == max(khist))[1] # take the first if more than one
  cat("most frequent k (kfreq):", kfreq, "\n") 
  khist <- khist/length(krec)
  ##------------------------------------------
  ## par(ask=TRUE)
  if (plot) {
    old.par <- par(no.readonly=TRUE)
    ## x11()
    ## par(ask=TRUE)
    ## plot(khist, ty="h", xlim=c(0,kmax), ylim=c(0,1))
    ## barplot(khist, xlim=c(0,kmax), ylim=c(0,1))
    narg1 <- 3
    narg2 <- min(length(khist), 12)
    par(mfrow=c(2,2), mar=c(5,1,2,1) + 0.1)
    barplot(khist[narg1:narg2], names.arg=narg1:narg2, main="DPM",
            axes=FALSE, xlab="k")
  }
  ##------------------------------------------
  ## Density estimate using kfreq simulated components at each iteration
  ## cidens     <- numeric(ngrid)
  t <- seq(0, 1, length=ngrid)
  cidens <- matrix(0, nr=kfreq, nc=ngrid)
  ydens <- numeric(ngrid)
  ## muk <- NULL
  ## vark <- NULL
  km <- 0
  for (i in 1:length(krec)) {
    nk <- krec[i]
    w  <- wrec[i,1:nk] 
    mu <- phirec[i,1:nk]
    var <- varrec[i]
    ## densities for simulations with kfreq components
    ## !!! beware label switching
    ## !!! guarantee of equal number of components (with sort) is not enough 
    ## !!! many mu simulated values are centered on different ranges
    if (nk == kfreq) {
      km <- km + 1
      imu <- sort(mu, ind=TRUE)
      loc <- imu$ix
      ## reordering of mu simulated values
      ## muk <- cbind(muk,mu)
      ## vark <- c(vark,var)
      mu <- mu[loc]
      w <- w[loc]
      ## muk <- cbind(muk,mu)
      dens <- matrix(0, nr=nk, nc=ngrid) 
      for (j in 1:nk) {
        dens[j,] <- dnorm(t, mu[j], sqrt(var)) # unique var
      }
      ydens <- ydens + as.vector(w%*%dens)    # y density for k=kfreq
      cidens <- cidens + dens*w               # component densities
    }
  }
  cidens <- cidens/km
  ydens <- ydens/km
  cat("n.iterations with kfreq:",km,"\n") 
  ##-------------------------
  if (plot) {
    ## Density plots
    par(mar=c(5,3,2,1) + 0.1)
    hist(x, probability=TRUE, breaks="Freedman-Diaconis",
         xlab="scaled intensity", ylab="", main="DPM")
    ## lines(t, ydens)
    lines(t, ydens, col="red")
    ##--------------
    par(mar=c(5,3,2,1) + 0.1)
    ##
    for (k in 1:kfreq) {
      if (k == 1) {
        ## plot(t, cidens[k,], ty="l", lty=k, xlim=c(0,1),ylim=c(-0.2,4),
        plot(t, cidens[k,], type="l", col=k, xlim=c(0,1),
             ylim=c(0, max(cidens)), axes=FALSE,
             xlab="scaled intensity", ylab="", main="DPM")
      } else {
        lines(t, cidens[k,], col=k)
      }
    }
    axis(1, seq(0,1,by=0.2), seq(0,1,by=0.2))
    axis(2, seq(0,4,by=0.5), seq(0,4,by=0.5))
  }
  ## dev <- dev.cur()
  ##--------------------------------------------
  ## z values : evaluate from w-,mu-,var-results 
  zall <- matrix(0, nr=kfreq, nc=ngrid) 
  for (k in 1:kfreq) {
    zi <- cidens[k,]/ydens 
    zall[k,] <- zi
    if (plot && k == 1) {
      plot(t, zi, lty="dashed", col=k, type="l", xlim=c(0,1), ylim=c(0,1),
           xlab="", ylab="", main="Cluster partition")
      ## !!! applying round of zi may leave holes in z
      ## lines(t,round(zi),lty="dotted",col=k)
      ## par(new=TRUE)
    } else {
      lines(t, zi, lty="dashed", col=k)
    }
  }
  ##--------------------------------------------
  ## better: choosing z by the most probable zi
  choose <- function(x) { return( which(x == max(x))) }
  z <- apply(zall, 2, choose)
  z <- relabel(z)
  if (plot) { 
    points(t, rep(0,ngrid), col=z, pch=20)
    ## apply also to previous plot
    ## cdev <- dev.cur()
    ## dev.set(dev)
    ## for (i in 1:max(z)) {
    ##   jj <- range(which(z == i))
    ##   lines(t[jj], c(-0.18, -0.18), lty=i+1, lwd=2.2)
    ##   points(t[jj], c(-0.18, -0.18), pch=22)
    ## }
    ## dev.set(dev) # restore dev
    par(old.par)
  }
  ##--------------------------------------------
  ## par(ask=FALSE)
  invisible(z)
}

