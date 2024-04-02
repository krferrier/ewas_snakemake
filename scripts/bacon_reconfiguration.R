# Script for running BACON

# Import libraries
library(argparse)
library(BiocParallel)
library(dplyr)
library(data.table)
library(bacon)


set.seed(12345)
biases <- runif(6, -0.2, 0.2)
inflations <- runif(6, 1, 1.3)
es <- matrix(nrow=5000, ncol=6)
for(i in 1:6)
  es[,i] <- rnormmix(5000, c(0.9, biases[i], inflations[i], 0, 4, 1), shuffle=FALSE)
se <- replicate(6, 0.8*sqrt(4/rchisq(5000,df=4)))
colnames(es) <- colnames(se) <- LETTERS[1:ncol(se)]
rownames(es) <- rownames(se) <- 1:5000
head(rownames(es))

# run bacon 
bc1 <- bacon(teststatistics = NULL,
            effectsizes = es,
            standarderrors = se,
            verbose = T,
            trim = 0.999)
# run bacon again on the same exact data
bc2 <- bacon(teststatistics = NULL,
             effectsizes = es,
             standarderrors = se,
             verbose = T,
             trim = 0.999)
# check if outputs are the same
test.1 <- list(inflation(bc1), bias(bc1), pval(bc1))
test.2 <- list(inflation(bc1), bias(bc2), pval(bc2))
identical(test.1[[1]], test.2[[1]])
identical(test.1[[2]], test.2[[2]])
identical(test.1[[3]], test.2[[3]])

RNGseed = 3
.bacon <- function(i, object, niter, nbins, trim, level, verbose, priors){
  
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
                                epsilon = list(gamma = c(90.0, 5.0, 5.0)))){
  
  
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
        object@traces[,,i] <- .bacon(i, object, niter, nbins, trim, level, verbose, priors)
    }
    else{
      nworkers <- min(c(nset, nworkers))
      message("Detected ", nworkers, " workers!\n Running in parallel!")
      ret <- bplapply(1:ncol(tstat(object)), .bacon, object=object,
                      niter=niter, nbins=nbins, trim=trim, level=level, verbose=verbose, priors=priors, BPOPTIONS = bpoptions(RNGseed = RNGseed))
      object@traces <- simplify2array(ret)
    }
  } else
    object@traces[,,1] <- .bacon(1, object, niter, nbins, trim=trim, level, verbose, priors)
  
  ##summarize traces
  for(i in 1:ncol(tstat(object)))
    object@estimates[i,] <- apply(object@traces[-c(1:object@nburnin),,i], 2, mean)
  
  return(object)
}

bacon_bmi_ewas <- list()
set.seed(2974)
for (i in 1:length(bmi_ewas)){
  name <- names(bmi_ewas)[i]
  subset <- bmi_ewas[[name]]
  # Run bacon
  bacon_res <- bacon(teststatistics = subset$statistic,
                     verbose = T,
                     trim = 0.999)
  bacon_effects <- bacon(teststatistics = NULL,
                         effectsizes = subset$estimate,
                         standarderrors = subset$std.error,
                         verbose = T,
                         trim = 0.999)
  # Add bacon results to subset dataframe
  subset$bacon.pval <- bacon::pval(bacon_res)
  subset$bacon.statistic <- bacon::tstat(bacon_res)
  subset$bacon.es <- bacon::es(bacon_effects)
  subset$bacon.se <- bacon::se(bacon_effects)
  b.lambda <- QCEWAS::P_lambda(subset$bacon.pval)
  lambda <- QCEWAS::P_lambda(subset$p.value)
  subset$lambda <-lambda
  subset$b.lambda <- b.lambda
  bacon_bmi_ewas[[name]] <- subset
  cat("lambda before bacon: ", lambda, "\n", "lambda after bacon: ", b.lambda, "\n")
  # Run performance tests
  traces(bacon_res)
  posteriors(bacon_res)
  fit(bacon_effects)
  qqman::qq(subset$p.value, main = "Before Bacon")
  qqman::qq(subset$bacon.pval, main = "After Bacon")
}
