#' Calculate derivatives of option prices with respect to state values, using a laguerre quadrature
#' @param params.Q The Q-parameters of the stock-price process
#' @param strikeMat TxKxS array, containing the strikes at which the options are valued
#' @param mkt list describing the market structure (times to maturity, interest rate, dividend yield, current stock price)
#' @param stateVec SxN.factors matrix containing the states at which the pricing/derivatives are calculated
#' @param alpha parameter of the laguerre quadrature
#' @param N number of quadrature points to take
#' @param preCalc a list containing preCalculated values of the ODE solutions (useful if derivatives are reevaluated many times for different states, but the same parameter)
#' @param N.factors The number of stochastic factors we have
#' @export
#' @return A list containing elements call, put and otm, which in turn contain two lists of arrays of size TxKxS which show the sensitivity to state 1 or 2 at the given maturity,strike,state

gaussLaguerre.stateDeriv <- function(params.Q, strikeMat, mkt, stateVec, rtol=1e-13, alpha=0, N=64, preCalc = NULL, N.factors=3) {
  # sanity check
  stopifnot(ncol(stateVec) == N.factors)
  
	T <- dim(strikeMat)[1]
	K <- dim(strikeMat)[2]
	S <- dim(strikeMat)[3]
	
	# if CF values are not pre-calculated then calculate them
	
	if (is.null(preCalc)) {
		
		require('statmod')
		weightAndnodes <- gauss.quad(N,'laguerre',alpha=alpha)
		uVec <- weightAndnodes$nodes    
		cfVal <- twoFactorJumpODEsSolve(u=cbind(matrix(1i*uVec,ncol=1),matrix(0,N,N.factors)),params.Q,mkt,rtol=rtol,N.factors=N.factors)
		
	} else {
		
		weightAndnodes <- preCalc$weightAndnodes
		uVec <- weightAndnodes$nodes
		cfVal <- preCalc$cfVal
	}
	
	F <- exp(mkt$t*(mkt$r-mkt$q))
	
	V <- rowSums(stateVec)
	
	
	# create affine CF array TxKxSxN dimensions
	cfODE <- aperm(array(rep(as.vector(cfVal[,,"a"]),K*S),c(N,T,K,S)),c(2,3,4,1))
  for (nn in 1:N.factors) {
    cfODE <- cfODE + aperm(array(rep(rep(as.vector(cfVal[,,paste0("b",nn)]),S)*rep(stateVec[,nn],each=N*T),K),c(N,T,S,K)),c(2,4,3,1))
  }
	cfODE <- exp(cfODE)
	
	# now create derivatives for both factors
	callList <- putList <- otmList <- list()
	
	for (f in 1:N.factors) {
		cfODE.scaled <- cfODE * aperm(array(rep(rep(as.vector(cfVal[,,paste0("b",f)]),S),K),c(N,T,S,K)),c(2,4,3,1))
		
		dCF <-  (cfODE.scaled) / (1i * aperm(array(rep(uVec^(alpha+1),T*K*S),c(N,T,K,S)),c(2,3,4,1))* (1-1i* aperm(array(rep(uVec,T*K*S),c(N,T,K,S)),c(2,3,4,1))))
		
		cExp <- -1i*array(rep(as.vector(strikeMat),N),c(T,K,S,N))
		cExp <- cExp*(aperm(array(rep(uVec,T*K*S),c(N,T,K,S)),c(2,3,4,1)))
		cExp <- exp(cExp)
		
		# create u Mat
		dC <- dCF*cExp
		
		dC <- dC * aperm(array(rep(weightAndnodes$weights*exp(uVec),T*K*S),c(N,T,K,S)),c(2,3,4,1))
		
		call <- rowSums(array(dC[,,,],c(T,K,S,N)),dims=3)
		
		df1 <- array(rep(exp(-mkt$r*mkt$t),K*S),c(T,K,S))
		df2 <- array(rep(exp(-mkt$q*mkt$t),K*S),c(T,K,S))
		
		call <- -exp(strikeMat)*df1/pi * Re(call)
		
		# now use put-call parity for put prices
		put <- call 
		
		# now get otm
		otm <- call
		otm[strikeMat < 0] <- put[strikeMat < 0]
		
		otmList[[f]] <- otm
		callList[[f]] <- call
		putList[[f]] <- put
	}
	
	if("day" %in% names(mkt)){
		return(list(day = mkt$day, call = callList, put = putList, otm = otmList))
	}
	else{
		return(list(call = callList, put = putList, otm = otmList))
	}
}

