#' @rdname cosPricing
#' @title Cosine Transform pricing
#' @description Implementation of Fang and Oosterlee's (2008) Fourier-cosine expansion method of pricing European derivatives.
#' @param strikeMat array of size \code{TxKxS} of relative log-strikes
#' @param N number of points of integration.
#' @param mkt data.frame with \code{T} rows and fields: r -- risk-free rate, q -- dividend yield, t -- option maturity.
#' @param intLim \code{numeric(2)} integration limits in terms of log-strikes, see Fang and Oosterlee for some optimal choices.
#' @param payCoeffFoo function which accepts arguments \code{intLim, Nterms, strikeMat} and calculates the cosine expansion coefficents of a given derivative contract. This function determines what you are pricing. Use \code{callCosCoeffs} to price call (and put, via parity) options.
#' @param N.factors integer, number of stochastic volatility factors, argument for \code{charFun}. If your \code{charFun} doesn't accept such an argument, for example you're pricing in the Black-Scholes model, use \code{N.factors = 0}.
#' @param preCalc optional, array of size \code{NxTxKxS} with pre-calculated characteristic function values for integration. If you provide it, you can ommit \code{charFun} as it will not be called.
#' @param ... arguments to charFun required for pricing (state variables, parameters, etc.)
#' @return list of arrays of size \code{T x K x S} with relative prices of call options, put options, out of the money options.
#' @details The integration method is described in detail in F. Fang and C.W. Oosterlee, ``A Novel Pricing Method for European Options Based on Fourier-Cosine Series Expansions'', accessible at http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf.
#' @export

cosTransformPricer <- function(strikeMat, mkt, N=120, intLim, payCoeffFoo = callCosCoeffs, N.factors = 3, charFun = affineCF, preCalc = NULL, ...) {
  
  dimStrikeMat <- dim(strikeMat)
  T <- dimStrikeMat[1]
  K <- dimStrikeMat[2]
  if(length(dimStrikeMat)!=3){
    S <- 1
    strikeMat <- array(as.numeric(strikeMat), dim = c(T,K,S))
  } else {
    S <- dimStrikeMat[3]
  }
  
  cf.u <- 0:(N-1) * pi / diff(intLim) * 1i
  if(is.null(preCalc)){
    cf.values <- charFun(u = cbind(cf.u,matrix(0,nrow=N,ncol=N.factors)), t.vec = NULL, N.factors = N.factors, mkt = mkt, ...)
    cf.values <- array(cf.values, dim = c(N,T,S,K))
    cf.values <- aperm(a = cf.values, perm = c(1,2,4,3))
  } else {
    cf.values <- preCalc
  }
  
  disc.factor <- exp(-(mkt$r)* mkt$t)
  disc.factor <- array(disc.factor, dim = c(T,K,S))
  
  pay.coeffs <- payCoeffFoo(intLim, N, strikeMat)
  eij.mat <- eijMat(intLim = intLim, Nterms = N, strikeMat = strikeMat)
  
  call <- cf.values * eij.mat * pay.coeffs
  loc.array <- colSums(call[-1,,,,drop=FALSE])
  call <- Re(0.5*call[1,,,,drop=FALSE] + array(colSums(call[-1,,,]), dim = c(1,dim(loc.array))))
  call <- drop(call)
  call <- array(call, dim = c(T,K,S))
  call <- disc.factor * call
  
  df1 <- array(rep(exp(-mkt$r*mkt$t),K*S),c(T,K,S))
  df2 <- array(rep(exp(-mkt$q*mkt$t),K*S),c(T,K,S))
  
    # now use put-call parity for put prices
  put <- call + exp(strikeMat)*df1 - 1*df2
  
  # now get otm
  otm <- call
  otm[strikeMat < 0] <- put[strikeMat < 0]
  
  return(list(call = call, put = put, otm = otm))
}

#' @describeIn cosPricing
#' @export

callCosCoeffs <- function(intLim, Nterms = 120, strikeMat){
  
  # Internal helper functions
  psiFoo <- function(N,c,d){
    kk <- 1:(N-1)
    psi <- matrix(0,nrow = N, ncol = length(c))
    psi[1,] <- (d - c) * rep(1,length(c))
    psi[-1,] <- kronecker(matrix((b-a)/(kk*pi),ncol = 1), matrix(1, ncol = length(c))) * (kronecker(matrix(sin(kk*pi*(d-a)/(b-a)), ncol = 1),matrix(1,ncol=length(c))) - sin(kronecker(matrix(kk*pi,ncol=1),matrix((c-a)/(b-a),ncol=length(c)))))
    return(psi)
  }
  
  chiFoo <- function(N,c,d){
    kk <- 0:(N-1)
    chi0 <- kronecker(matrix(1/(1 + (kk*pi/(b-a))^2),ncol=1), matrix(1,ncol = length(c)))
    chi1 <- kronecker(matrix(cos(kk * pi * (d-a)/(b-a))*exp(d),ncol=1),matrix(1,ncol=length(c)))
    chi2 <- cos(kronecker(matrix(kk * pi,ncol=1), matrix((c-a)/(b-a),ncol=length(c)))) * kronecker(matrix(1,nrow=length(kk)),matrix(exp(c),nrow=1))
    chi3 <- kronecker(matrix(sin(kk * pi * (d-a)/(b-a))*exp(d)*kk*pi/(b-a),ncol=1),matrix(1,ncol=length(c)))
    chi4 <- sin(kronecker(matrix(kk*pi,ncol=1), matrix((c-a)/(b-a),nrow=1))) * kronecker(matrix(kk*pi/(b-a),ncol=1), matrix(exp(c),nrow=1));
    chi <- chi0 * (chi1-chi2+chi3-chi4)
    return(chi)
  }
    
  a <- intLim[1]
  b <- intLim[2]
  
  dimStrikeMat <- dim(strikeMat)
  T <- dimStrikeMat[1]
  K <- dimStrikeMat[2]
  if(length(dimStrikeMat)!=3){
    S <- 1
  } else {
    S <- dimStrikeMat[3]
  }
  
  strikeMat <- as.numeric(strikeMat)
  
  VkVec <- 2/(b-a) * (chiFoo(Nterms,strikeMat,b)-psiFoo(Nterms,strikeMat,b) * kronecker(matrix(1,nrow=Nterms),exp(t(strikeMat))))
  
  VkVec <- array(data = as.numeric(VkVec), dim = c(Nterms,T,K,S))
  
  return(VkVec)
}

#' @describeIn cosPricing
#' @export

putCosCoeffs <- function(intLim, Nterms = 120, strikeMat){
  
  # Internal helper functions
  psiFoo <- function(N,c,d){
    kk <- 1:(N-1)
    psi <- matrix(0,nrow = N, ncol = length(c))
    psi[1,] <- (d - c) * rep(1,length(c))
    psi[-1,] <- kronecker(matrix((b-a)/(kk*pi),ncol = 1), matrix(1, ncol = length(c))) * (kronecker(matrix(sin(kk*pi*(d-a)/(b-a)), ncol = 1),matrix(1,ncol=length(c))) - sin(kronecker(matrix(kk*pi,ncol=1),matrix((c-a)/(b-a),ncol=length(c)))))
    return(psi)
  }
  
  chiFoo <- function(N,c,d){
    kk <- 0:(N-1)
    chi0 <- kronecker(matrix(1/(1 + (kk*pi/(b-a))^2),ncol=1), matrix(1,ncol = length(c)))
    chi1 <- kronecker(matrix(cos(kk * pi * (d-a)/(b-a))*exp(d),ncol=1),matrix(1,ncol=length(c)))
    chi2 <- cos(kronecker(matrix(kk * pi,ncol=1), matrix((c-a)/(b-a),ncol=length(c)))) * kronecker(matrix(1,nrow=length(kk)),matrix(exp(c),nrow=1))
    chi3 <- kronecker(matrix(sin(kk * pi * (d-a)/(b-a))*exp(d)*kk*pi/(b-a),ncol=1),matrix(1,ncol=length(c)))
    chi4 <- sin(kronecker(matrix(kk*pi,ncol=1), matrix((c-a)/(b-a),nrow=1))) * kronecker(matrix(kk*pi/(b-a),ncol=1), matrix(exp(c),nrow=1));
    chi <- chi0 * (chi1-chi2+chi3-chi4)
    return(chi)
  }
  
  a <- intLim[1]
  b <- intLim[2]
  
  dimStrikeMat <- dim(strikeMat)
  T <- dimStrikeMat[1]
  K <- dimStrikeMat[2]
  if(length(dimStrikeMat)!=3){
    S <- 1
  } else {
    S <- dimStrikeMat[3]
  }
  
  strikeMat <- as.numeric(strikeMat)
  
  VkVec <- 2/(b-a) * (-chiFoo(Nterms,strikeMat,b)+psiFoo(Nterms,strikeMat,b) * kronecker(matrix(1,nrow=Nterms),exp(t(strikeMat))))
  
  VkVec <- array(data = as.numeric(VkVec), dim = c(Nterms,T,K,S))
  
  return(VkVec)
}

#' @describeIn cosPricing
#' @export

digitalCosCoeffs <- function(intLim, Nterms, strikeMat){

  # Internal helper functions
  psiFoo <- function(N,c,d){
    kk <- 1:(N-1)
    psi <- matrix(0,nrow = N, ncol = length(c))
    psi[1,] <- (d - c) * rep(1,length(c))
    psi[-1,] <- kronecker(matrix((b-a)/(kk*pi),ncol = 1), matrix(1, ncol = length(c))) * (kronecker(matrix(sin(kk*pi*(d-a)/(b-a)), ncol = 1),matrix(1,ncol=length(c))) - sin(kronecker(matrix(kk*pi,ncol=1),matrix((c-a)/(b-a),ncol=length(c)))))
    return(psi)
  }
  
  a <- intLim[1]
  b <- intLim[2]
  
  dimStrikeMat <- dim(strikeMat)
  T <- dimStrikeMat[1]
  K <- dimStrikeMat[2]
  if(length(dimStrikeMat)!=3){
    S <- 1
  } else {
    S <- dimStrikeMat[3]
  }
  
  strikeMat <- as.numeric(strikeMat)
  
  VkVec <- 2/(b-a) * psiFoo(Nterms,strikeMat,b)
  
  VkVec <- array(data = VkVec, dim = c(Nterms,T,K,S))
  
  return(VkVec)
}



eijMat <- function(intLim, Nterms, strikeMat){
  
  a <- intLim[1]
  b <- intLim[2]
  
  dimStrikeMat <- dim(strikeMat)
  T <- dimStrikeMat[1]
  K <- dimStrikeMat[2]
  if(length(dimStrikeMat)!=3){
    S <- 1
  } else {
    S <- dimStrikeMat[3]
  }
  
  vInt <- 0:(Nterms-1)
  
  eij <- exp(-1i * vInt * pi * a / (b-a))
  eij <- array(eij, dim = c(Nterms,T,K,S))
  
  return(eij)
}