# Modify the original bacon function source code to control RNG of bpapply.
# There are only two modifications to the source code: 
#     1. The addition of a variable for setting a seed to control the random 
#        number generation of the BiocParallel bpapply() function (it ignores 
#        global 'set.seed()'). This variable is called 'rng' and the default
#        value is set to NULL.
#     2. The inclusion of BPOPTIONS in the bpapply() call that sets the 
#        RNGseed = rng.

.bacon <- function(i, object, niter, nbins, trim, level, verbose, priors, globalSeed){
  if(!is.null(globalSeed))
    set.seed(globalSeed)
  ##TODO add some kind of trimming?
  tstats <- tstat(object)[,i]
  
  if (object@na.exclude) 
    tstats <- tstats[!is.na(tstats)]
  
  medy <- median(tstats)
  mady <- mad(tstats)
  binned <- !is.null(nbins)
  
  if(verbose & binned)
    message("Use multinomial weighted sampling...")
  else if(verbose & !binned)
    message("Use fast random weighted sampling...")
  
  if(!is.null(nbins)){
    ##q <- range(tstats) ##sensitive to outlying values
    q <- quantile(tstats, prob=c(1 - trim, trim))
    tstats <- tstats[tstats > q[1] & tstats < q[2]]
    breaks <- seq(q[1], q[2], length = nbins + 1)
    x <- breaks[-c(nbins+1)] + 0.5*(q[2] - q[1])/(nbins) #identical to h$mids
    h <- hist(tstats, breaks = breaks, plot=FALSE)
    w <- h$counts
    tstats <- h$mids
  }
  else
    w <- 1.0+0*tstats
  
  output <- .C("bacon",
               y = as.vector(tstats, mode="double"),
               w = as.vector(w, mode="double"),
               medy = as.double(medy),
               mady = as.double(mady),
               n = as.integer(length(tstats)),
               niter = as.integer(niter),
               level = as.double(level),
               binned = as.integer(binned),
               verbose = as.integer(verbose),
               gibbsmu = vector(3*niter, mode="double"),
               gibbssig = vector(3*niter, mode="double"),
               gibbsp = vector(3*niter, mode="double"),
               alpha = as.double(priors$sigma$alpha),
               beta = as.double(priors$sigma$beta),
               lambda = as.vector(priors$mu$lambda, mode="double"),
               tau = as.vector(priors$mu$tau, mode="double"),
               gamma = as.vector(priors$epsilon$gamma, mode="double"),
               PACKAGE="bacon")
  
  gibbsp <- matrix(output$gibbsp, ncol=3, byrow=TRUE,
                   dimnames=list(NULL, paste0("p.", 0:2)))
  gibbsmu <- matrix(output$gibbsmu, ncol=3, byrow=TRUE,
                    dimnames=list(NULL, paste0("mu.", 0:2)))
  gibbssig <- matrix(output$gibbssig, ncol=3, byrow=TRUE,
                     dimnames=list(NULL, paste0("sigma.", 0:2)))
  
  return(cbind(gibbsp, gibbsmu, gibbssig))
}
bacon <- function(teststatistics=NULL, effectsizes=NULL, standarderrors=NULL,
                  niter=5000L, nburnin = 2000L, nbins=1000, trim =0.999, level=0.05, 
                  na.exclude=FALSE, verbose=FALSE,
                  priors = list(sigma = list(alpha = 1.28,
                                             beta = 0.36), ##original uses 0.36*mad(teststatistics)
                                mu = list(lambda = c(0.0, 3.0, -3.0),
                                          tau = c(1000.0, 100.0, 100.0)),
                                epsilon = list(gamma = c(90.0, 5.0, 5.0))),
                  globalSeed = 42, # if set to NULL, randomization will occur for sequential and parallel bacon calls
                  parallelSeed = 42){ # if input statistics are a matrix and globalSeed=NULL, setting parallelSeed=NULL will allow randomization across parallel bacon calls and across separate calls to bacon. 
  ##create new Bacon-object
  object <- new("Bacon", teststatistics = teststatistics,
                effectsizes = effectsizes,
                standarderrors = standarderrors,
                na.exclude = na.exclude,
                niter = niter,
                nburnin = nburnin,
                priors = priors)
  
  nset <- ncol(tstat(object))    
  ##run the Gibbs Sampler
  if(ncol(tstat(object)) > 1){
    nworkers <- bpworkers(bpparam())
    if(nworkers <= 1) {
      if(nset > 1)
        message("Did you registered a biocparallel back-end?\n Continuing serial!")
      for(i in 1:ncol(tstat(object)))
        object@traces[,,i] <- .bacon(i, object, niter, nbins, trim, level, verbose, priors, globalSeed)
    }
    else{
      nworkers <- min(c(nset, nworkers))
      message("Detected ", nworkers, " workers!\n Running in parallel!")
      if(!is.null(parallelSeed)){
        ret <- bplapply(1:ncol(tstat(object)), .bacon, object=object,
                      niter=niter, nbins=nbins, trim=trim, level=level, verbose=verbose, priors=priors, 
                      globalSeed=globalSeed, BPOPTIONS = bpoptions(RNGseed = parallelSeed))
      } else {
        ret <- bplapply(1:ncol(tstat(object)), .bacon, object=object,
                        niter=niter, nbins=nbins, trim=trim, level=level, verbose=verbose, priors=priors, 
                        globalSeed=globalSeed)
      }
      object@traces <- simplify2array(ret)
    }
  } else
    object@traces[,,1] <- .bacon(1, object, niter, nbins, trim=trim, level, verbose, priors, globalSeed)
  
  ##summarize traces
  for(i in 1:ncol(tstat(object)))
    object@estimates[i,] <- apply(object@traces[-c(1:object@nburnin),,i], 2, mean)
  
  return(object)
}