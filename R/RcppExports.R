# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#'@export
glPricer_cpp <- function(strikeMat, mkt, glWts, glNodes, cfVals, Nfactors, alpha = 0, sigmaRef = 1.0) {
    .Call('_transformOptionPricer_glPricer_cpp', PACKAGE = 'transformOptionPricer', strikeMat, mkt, glWts, glNodes, cfVals, Nfactors, alpha, sigmaRef)
}

