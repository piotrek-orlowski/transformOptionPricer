#' @title Laguerre quadrature pricing.
#' @description Calculate option prices using a Laguerre quadrature of the difference between a reference characteristic function (Black-Scholes for high volatility) and user-defined characteristic function. Only European call and put options.
#' @param strikeMat array of size \code{TxKxS} of relative log-strikes
#' @param N number of points of integration.
#' @param sigma.ref variance (volatility squared) value for the reference characteristic function, length \code{S}. If not provided, \code{...} will be checked for existence of a state matrix and \col{rowSums} of variance states will be taken.
#' @param mkt data.frame with \code{T} rows and fields: r -- risk-free rate, q -- dividend yield, t -- option maturity.
#' @param alpha parameter of the laguerre quadrature.
#' @param N.factors integer, number of stochastic volatility factors, argument for \code{charFun}. If your \code{charFun} doesn't accept such an argument, for example you're pricing in the Black-Scholes model, use \code{N.factors = 0}.
#' @param preCalc optional a list containing preCalculated values of the characteristic function, and the quadrature parameters (useful if derivatives are reevaluated many times for different states, but the same parameter)
#' @param ... arguments to charFun required for pricing (state variables, parameters, etc.)
#' @return list of arrays of size \code{T x K x S} with relative prices of call options, put options, out of the money options.
#' @details In extensive tests this pricer yields results comparable to the Fourier-cosine series pricer, with fewer characteristic function evaluations. Care should be taken at very high values of variance factors.
#' @export
#' 
gaussLaguerrePricer <- function(strikeMat, mkt, alpha=0, N=64, sigma.ref = NULL, N.factors = 3, charFun, preCalc = NULL, ...) {
  
  dimStrikeMat <- dim(strikeMat)
  T <- dimStrikeMat[1]
  K <- dimStrikeMat[2]
  if(length(dimStrikeMat)!=3){
    S <- 1
    strikeMat <- array(as.numeric(strikeMat), dim = c(T,K,S))
  } else {
    S <- dimStrikeMat[3]
  }
	
	# if CF values are not pre-calculated then calculate them
	if (is.null(preCalc)) {
	  
	  require('statmod')
	  weightAndnodes <- gauss.quad(N,'laguerre',alpha=alpha)
	  uVec <- weightAndnodes$nodes    
	  # cfVal <- twoFactorJumpODEsSolve(u=cbind(matrix(1i*uVec,ncol=1),matrix(0,N,N.factors)),params.Q,mkt,rtol=rtol,N.factors=N.factors)
	  cfVal <- charFun(u = cbind(matrix(1i*uVec,ncol=1),matrix(0,N,N.factors)), t.vec = mkt$t, N.factors = N.factors, ...)
	  
	} else {
	  
	  weightAndnodes <- preCalc$weightAndnodes
	  uVec <- weightAndnodes$nodes
	  cfVal <- preCalc$cfVal
	}
	
	F <- exp(mkt$t*(mkt$r-mkt$q))
	
	if(!is.null(sigma.ref)){
	  V <- sigma.ref
	} else {
	  dots.arg <- as.list(substitute(list(...)))[-1L]
	  stateVec <- eval(dots.arg[["v.0"]])
	  V <- if(!is.null(dim(stateVec))) rowSums(stateVec) else stateVec
	  V[which(V > 0.5)] <- 0.5
	}
	
	# create B-S c.f. array TxKxSxN dimensions
	M1 <- matrix(log(F),T,1) %*% matrix(1,1,S)
	M2 <- matrix(mkt$t,T,1) %*% matrix(V,1,S)
	
	M3 <- array(rep(as.vector(2*M1-M2),N)*rep(uVec,each = T*S),c(T,S,N))
	M4 <- array(rep(as.vector(-M2),N)*rep(uVec^2,each = T*S),c(T,S,N))
	
	cfBS <- aperm(array(exp(0.5 * rep(as.vector(M3*1i+M4),K)),c(T,S,N,K)),c(1,4,2,3))

  
	# create B-S price matrix
	s <- aperm(array(rep(as.vector(M2),K),c(T,S,K)),c(1,3,2))
	m <- array(rep(F,S*K),c(T,K,S))

	d1 <-1/sqrt(s)*(log(m)-strikeMat+0.5*s)
	d2 <- d1 - sqrt(s)

	pBS <- exp(array(-rep(mkt$t*mkt$r,K*S),c(T,K,S)))*(pnorm(d1)*array(rep(F,K*S),c(T,K,S)) - exp(strikeMat)*pnorm(d2))

# 	# create affine CF array TxKxSxN dimensions
# 	cfODE <- aperm(array(rep(as.vector(cfVal[,,"a"]),K*S),c(N,T,K,S)),c(2,3,4,1))
	cfVal <- aperm(array(cfVal, dim = c(N,T,S,K)), c(2,4,3,1))
#   for (nn in 1:N.factors) {
# 	  cfODE <- cfODE + aperm(array(rep(rep(as.vector(cfVal[,,paste0("b",nn)]),S)*rep(stateVec[,nn],each=N*T),K),c(N,T,S,K)),c(2,4,3,1))
#   }
# 	cfODE <- exp(cfODE)
	
	dCF <-  (cfVal - cfBS) / (1i * aperm(array(rep(uVec^(alpha+1),T*K*S),c(N,T,K,S)),c(2,3,4,1))* (1-1i* aperm(array(rep(uVec,T*K*S),c(N,T,K,S)),c(2,3,4,1))))

	cExp <- -1i*array(rep(as.vector(strikeMat),N),c(T,K,S,N))
	cExp <- cExp*(aperm(array(rep(uVec,T*K*S),c(N,T,K,S)),c(2,3,4,1)))
	cExp <- exp(cExp)

	# create u Mat
	dC <- dCF*cExp
	
	dC <- dC * aperm(array(rep(weightAndnodes$weights*exp(uVec),T*K*S),c(N,T,K,S)),c(2,3,4,1))
	call <- rowSums(array(dC[,,,],c(T,K,S,N)),dims=3)
	
	df1 <- array(rep(exp(-mkt$r*mkt$t),K*S),c(T,K,S))
	df2 <- array(rep(exp(-mkt$q*mkt$t),K*S),c(T,K,S))
		
	call <- -exp(strikeMat)*df1/pi * Re(call) + pBS
	
	
	# now use put-call parity for put prices
	put <- call + exp(strikeMat)*df1 - 1*df2
	
	# now get otm
	otm <- call
	otm[strikeMat < 0] <- put[strikeMat < 0]
	
	if("day" %in% names(mkt)){
		return(list(day = mkt$day, call = call, put = put, otm = otm))
	}
	else{
		return(list(call = call, put = put, otm = otm))
	}
}