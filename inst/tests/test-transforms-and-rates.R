### TRANSFORM PRICERS AND HANDLING INTEREST RATES / DIVIDEND YIELDS ###
library(transformOptionPricer)

# define Black-Scholes characteristic function
bscf <- function(u,sig,t.vec,N.factors,mkt){
  tau <- mkt$t
  return(exp((mkt$r - mkt$q - 0.5*sig^2)*u*tau + 0.5*sig^2*u^2*tau))
}

library(RQuantLib)

# define maturity, interest rate, dividend yield
t <- 1/12
r <- 0.03
q <- 0.1

# define number of integration points and limits for cosine transform pricing
N <- 384
intLim <- c(-5,5)

# number of pricing days
S <- 1

# BS IV
vol <- 1

# Strike, handling interest rates and yields
strikeMat <- t(matrix(seq(0,1,by=0.01),length(seq(0,1,by=0.01)),1))

mkt <- data.frame(p = 1, q = q, r= r , t= t)

mkt.zr <- data.frame(p=1,q=0.0,r=0.0,t= t)

# helper number
K <- dim(strikeMat)[2]

# calculate CF values for the cosine pricer with mkt structure containing interest rates
cf.u <- 0:(N-1) * pi / diff(intLim) * 1i
cf.values <- bscf(u = cbind(cf.u,matrix(0,nrow=N,ncol=1)), t.vec = NULL, N.factors = 0, mkt = mkt, sig = vol)
cf.values <- array(cf.values, dim = c(N,T,S,K))
cf.values <- aperm(a = cf.values, perm = c(1,2,4,3))

# calculate CF values for the cosine pricer while correcting for interest rates afterwards
cf.values.zr <- bscf(u = cbind(cf.u,matrix(0,nrow=N,ncol=1)), t.vec = NULL, N.factors = 0, mkt = mkt.zr, sig = vol) * exp(mkt$t*(mkt$r - mkt$q)*cf.u)
cf.values.zr <- array(cf.values.zr, dim = c(N,T,S,K))
cf.values.zr <- aperm(a = cf.values.zr, perm = c(1,2,4,3))

# price option with implicitly calculated interest rates effects
pr <- cosTransformPricer(strikeMat = strikeMat, mkt = mkt, N = dim(cf.values)[1], intLim = intLim, payCoeffFoo = callCosCoeffs, N.factors = 0, charFun = bscf, sig = vol, preCalc = cf.values)$otm[1]

# price option with explicitly added interest rates effects
pr.zr <- exp(-r*t) * cosTransformPricer(strikeMat = strikeMat, mkt = mkt.zr, N = dim(cf.values)[1], intLim = intLim, payCoeffFoo = callCosCoeffs, N.factors = 0, charFun = bscf, sig = vol, preCalc = cf.values.zr)$otm[1]

# true option price
pr.true <- EuropeanOption(type = 'call', underlying = 1, strike = exp(strikeMat[1]), dividendYield = q, riskFreeRate = r, maturity = t, volatility = vol)$value

# true option price
pr.true.zr <- exp(-r*t) * EuropeanOption(type = 'call', underlying = exp(t*(r-q)), strike = exp(strikeMat[1]), dividendYield = 0, riskFreeRate = 0, maturity = t, volatility = vol)$value

library(testthat)
sapply(list(pr, pr.zr, pr.true.zr), function(x){
  expect_equal(object = x, expected = pr.true)
})

### BATES 2006 model ###
load("d:/git/diveRgence/data/affine-parameters-bates2006.RData")

# Define CF with internally calculated interest rate effects
batesCF <- function(u,t.vec,N.factors,mkt,...){
  return(affineCF(u = u, params.Q = par.bates.svj1$Q, params.P = NULL, t.vec = NULL, v.0 = matrix(par.bates.svj1$P$`1`$eta), jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, CGF = F, mod.type = 'standard', mkt = mkt))
}

# Define CF with explicitly corrected interest rate effects
batesCF.zr <- function(u,t.vec,N.factors,mkt,...){
  res <- affineCF(u = u, params.Q = par.bates.svj1$Q, params.P = NULL, t.vec = NULL, v.0 = matrix(par.bates.svj1$P$`1`$eta), jumpTransform = getPointerToJumpTransform('expNormJumpTransform')$TF, N.factors = 1, CGF = F, mod.type = 'standard', mkt = mkt)
  res <- res * exp(u[,1] * t * (r - q))
  return(res)
}

N <- 1200

# Price by cosine transform with internally calculated interest rate effects
pr <- cosTransformPricer(strikeMat = strikeMat, mkt = mkt, N = dim(cf.values)[1], intLim = intLim, payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = batesCF)$otm

# Price by cosine transform with explicitly corrected interest rate effects
pr.zr <- exp(-r*t) * cosTransformPricer(strikeMat = strikeMat, mkt = mkt.zr, N = dim(cf.values)[1], intLim = intLim, payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = batesCF.zr)$otm

expect_equal(pr.zr, pr)

# Move to Gauss-Laguerre integration

# Number of SV factors in the Bates model
N.factors <- 1

# Number of integration points
N <- 128

# Gauss-Laugerre pricer and internally calculated interest rate effects 
pr.gl <- gaussLaguerrePricer(strikeMat = strikeMat, mkt = mkt, alpha = 0, N = N, sigma.ref = 0.2, N.factors = 1, charFun = batesCF, preCalc = NULL)$otm

# Gauss-Laugerre pricer with pre-calculated CF values
# GL quadrature points
weightAndnodes <- gauss.quad(N,'laguerre',alpha=0)
uVec <- weightAndnodes$nodes    

# evaluate bates CF with internal calculation of interest rate effects
cfVal <- batesCF(u = cbind(matrix(1i*uVec,ncol=1),matrix(0,N,N.factors)), t.vec = NULL, N.factors = N.factors, mkt = mkt)

# evaluate bates CF with explicit correction for interest rates
cfVal.zr <- batesCF.zr(u = cbind(matrix(1i*uVec,ncol=1),matrix(0,N,N.factors)), t.vec = NULL, N.factors = N.factors, mkt = mkt.zr)

precalc.gl.zr <- list(weightAndnodes = weightAndnodes, cfVal = cfVal.zr)

pr.gl.zr <- gaussLaguerrePricer(strikeMat = strikeMat, mkt = mkt, alpha = 0, N = N, sigma.ref = 0.2, N.factors = 1, charFun = NULL, preCalc = precalc.gl.zr)$otm

pr.gl.zr.cpp.top <- transformOptionPricer::glPricer_cpp(strikeMat = array(strikeMat,c(1,length(seq(0,1,by=0.01)),1)), mkt = as.matrix(mkt[,c(4,3,2,1)]), glWts = weightAndnodes$weights, glNodes = weightAndnodes$nodes, cfVals = precalc.gl.zr$cfVal, Nfactors = 1, sigmaRef = 0.2)

library(divergenceModelR)
pr.gl.zr.cpp.dm <- divergenceModelR::glPricer_cpp(strikeMat = array(strikeMat,c(1,length(seq(0,1,by=0.01)),1)), mkt = as.matrix(mkt[,c(4,3,2,1)]), glWts = weightAndnodes$weights, glNodes = weightAndnodes$nodes, cfVals = precalc.gl.zr$cfVal, Nfactors = 1, sigmaRef = 0.2)
