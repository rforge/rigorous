##
## Run DPM model
##

dpmixsim <-
function(x, M=1, a=1, b=2, upalpha=1, a0=2, b0=2, maxiter=4000, rec=1000,  fsave=NA, kmax=30)
{
    ## run dpmodel
    runif(1) # just for generating the seed if it doesn't exit
    n    <- length(x)
    cat("simulation length:",n,"\n")
    ptm  <- proc.time()
    res  <- .C("gibbsdpm",
        as.double(x),
        as.integer(n),
        as.double(M),
        as.double(a),
        as.double(b),
        as.integer(upalpha),
        as.double(a0),
        as.double(b0),
        as.integer(maxiter),
        as.integer(rec),
        as.integer(kmax),
        krec   = as.integer(rep(0,rec)), # number of simulated components
        wrec   = as.double (rep(0,rec*kmax)), # weights
        phirec = as.double (rep(0,rec*kmax)),
        varrec = as.double (rep(0,rec)))
    cat("\ntime of gibbdpm: ", (proc.time() - ptm)[1]/60,"\n")
    res <- list(krec=res$krec,wrec=res$wrec,phirec=res$phirec,varrec=res$varrec)
    if(!is.na(fsave)) {
      cat("saving simulation ", fsave, "...")
      save(res, file = fsave)
      cat("\n")
    }
    invisible(res)
}

