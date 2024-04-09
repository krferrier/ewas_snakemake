# Modify how bacon initializes:
setClass("Bacon",
         representation(teststatistics = "matrix",
                        effectsizes    = "matrix",
                        standarderrors = "matrix",
                        traces         = "array",
                        estimates      = "matrix",
                        priors         = "list",
                        niter          = "integer",
                        nburnin        = "integer",
                        na.exclude     = "logical"),
         prototype(teststatistics = matrix(1),
                   effectsizes    = matrix(1),
                   standarderrors = matrix(1),
                   traces         = array(1),
                   estimates      = matrix(1),
                   priors         = list(),
                   niter          = integer(1),
                   nburnin        = integer(1),
                   na.exclude     = logical(1)),
         validity = function(object){
           if(niter <= nburnin)
             return("niter should be > nburnin!")
           else
             return(TRUE)}
)

setMethod("initialize", "Bacon",
          function(.Object, teststatistics, effectsizes, standarderrors, 
                   niter, nburnin, priors, na.exclude) {
            # Change the first conditional statement to also check if effectsize and standarderrors are NULL before running with only tstats
            if(!is.null(teststatistics) & is.null(effectsizes) & is.null(standarderrors)){
              .Object@teststatistics <- as.matrix(teststatistics)
            }
            else if(is.null(teststatistics) & is.null(effectsizes) & is.null(standarderrors))
              stop("Need to provide test-statistics or both effect-sizes and standard errors!")
            # Add conditional statement that checks if effectsize and standarderrors were given, but no tstats, and in this case set tstats = effectsize/stderror
            else if(is.null(teststatistics) & !is.null(effectsizes) & !is.null(standarderrors)) {
              .Object@effectsizes <- as.matrix(effectsizes)
              .Object@standarderrors <- as.matrix(standarderrors)
              .Object@teststatistics <- as.matrix(effectsizes/standarderrors)
            }
            # If all three are given, use all three provided statistics
            else {
              .Object@effectsizes <- as.matrix(effectsizes)
              .Object@standarderrors <- as.matrix(standarderrors)
              .Object@teststatistics <- as.matrix(teststatistics)
            }
            if(!all(is.finite(.Object@teststatistics)) & !na.exclude)
              stop("Non finite value(s) in test statistics and na.exclude = FALSE!")
            
            .Object@estimates <- matrix(nrow=ncol(.Object@teststatistics), ncol=9,
                                        dimnames=list(colnames(.Object@teststatistics),
                                                      paste0(rep(c("p.", "mu.", "sigma."), each=3), 0:2)))
            .Object@traces <- array(dim=c(niter, 9, ncol(.Object@teststatistics)),
                                    dimnames=list(NULL,
                                                  paste0(rep(c("p.", "mu.", "sigma."), each=3), 0:2),
                                                  colnames(.Object@teststatistics)))
            
            .Object@niter <- niter
            .Object@nburnin <- nburnin
            .Object@priors <- priors
            .Object@na.exclude <- na.exclude
            .Object
          })