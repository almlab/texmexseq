dpoilog <- function(n, mu, sig, trunc=FALSE) {
  if (length(mu)>1 | length(sig)>1) stop('vectorization of mu and sig is not supported') 
  if (any((n[n!=0]/trunc(n[n!=0]))!=1)) stop('all n must be integers')
  if (!all(is.finite(c(mu,sig)))) stop('all parameters should be finite')
  if (sig<=0) stop('sig is not larger than 0')
  if (trunc & all(n == 0)) warning('all n will be truncated')
  
  pdfs <- .dpoilog(n, mu, sig)
  if (trunc) {
    pdfs <- pdfs / (1.0 - .dpoilog(0, mu, sig))
    pdfs[n == 0] <- 0.0
  }
  pdfs
}

.dpoilog <- function(n, mu, sig) {
  .C('poilog', as.integer(n), as.double(mu), as.double(sig),
     n_obs=as.integer(length(n)), val=double(length(n)))$val
}

rpoilog <- function(S, mu, sig, condS=FALSE, keep0=FALSE){
   sim <- function(nr){
     lamx <- rnorm(nr)
     x <- rpois(nr,exp(sig*lamx+mu))
     if (!keep0) x <- x[x>0]
     return(x)
   }
   
   if (S<1) stop('S is not positive')
   if (!is.finite(S)) stop('S is not finite')
   if ((S/trunc(S))!=1) stop('S is not an integer')
   if (sig<0) stop('sig is not positive')
   
   if (condS) {
     simVec <- vector('numeric',0)
     fac <- 2
     nr  <- S
     while (length(simVec)<S){
       simvals <- sim(nr*fac)
       simVec <- c(simVec,simvals)
       fac <- (1/(length(simvals)/(nr*fac)))*2
       fac <- ifelse(is.finite(fac),fac,1000)
       nr <- S-length(simvals)
     }
     simVec <- simVec[1:S]
   }
   
   else simVec <- sim(S)
   return(simVec)
}

poilogMLE <- function(n, start.mu, start.sig, trunc=TRUE, method='L-BFGS-B',
  control=list(fnscale=length(n)), ...) {

  if (is.matrix(n) | (is.data.frame(n))) {
    stop(paste('n has',ncol(n),'colums, supply a vector',sep=' ')) 
  }
  
  # truncate the input
  if (trunc) n <- n[n > 0]
  
  # guess start values
  startVals=c(mu=start.mu, sig=start.sig)

  # dereplicate
  un <- unique(n) # unique N values (counts)
  nr <- rep(NA,length(un)) # number of replicates each unique N has
  for (i in 1:length(un)){ nr[i] <- sum(n%in%un[i]) }
  
  # negative log likelihood is the objective function
  lnL <- function(z) {
    -sum((log(dpoilog(un, z[1], exp(z[2]), trunc=trunc)))*nr)
  }
  fit <- optim(startVals, lnL, control=control, method=method, lower=-20, upper=20, ...)
  
  if (fit$convergence!=0){
    if (fit$convergence==1) stop('the iteration limit has been reached!   try different startVals or increase maxit') 
    if (fit$convergence==10) stop('degeneracy of the Nelder Mead simplex ....')
    else stop(paste('unknown error in optimization', fit$message))
  } 
  
  fit$par <- c(as.numeric(fit$par), 1 - dpoilog(0, fit$par[1], exp(fit$par[2])))
  res <- list('par'=c('mu'=fit$par[1],'sig'=exp(fit$par[2])),'p'=fit$par[3],'logLval'=-fit$value,'gof'=NULL,boot=NULL)
  
  return(res)
}

texmex.fit <- function(n, start.mus=c(-2.0, -1.0, 0.0, 1.0, 2.0), start.sigs=rep(1.0, times=5)) {
    if (length(start.mus) != length(start.sigs)) stop('must provide same number of starting mu and sigma values')

    for (i in seq(from=1, to=(length(start.mus)))) {
         start.mu <- start.mus[i]
         start.sig <- start.sigs[i]
         tryCatch({
             res <- poilogMLE(n, start.mu, start.sig)
             return(res)
         }, error = function(e) warning(paste("fit", i, "failed", sep=" ")))
    }

    stop("all fit attempts failed")
}
