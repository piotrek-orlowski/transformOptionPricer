### test cos transform pricer with BS ###
library(transformOptionPricer)

bscf <- function(u,sig,t.vec,N.factors){
  tau <- t.vec
  return(exp(-1/2*sig^2*u*tau + 1/2*sig^2*u^2*tau))
}

library(RQuantLib)

strikeMat <- t(matrix(log(c(0.9,1.0,1.1)),3,1))

mkt <- data.frame(p=0,q=0.02,r=0.03,t=c(0.25,1/2))

pr <- cosTransformPricer(strikeMat = strikeMat, mkt = mkt, N = 200, intLim = c(-2,2), payCoeffFoo = callCosCoeffs, N.factors = 0, charFun = bscf, sig = 1)

pr.true <- list()
pr.true$call <- sapply(exp(as.numeric(strikeMat)), function(x) EuropeanOption(type = 'call', underlying = 1, strike = x, dividendYield = 0, riskFreeRate = 0, maturity = mkt$t, volatility = 1)$value)
pr.true$put <- sapply(exp(as.numeric(strikeMat)), function(x) EuropeanOption(type = 'put', underlying = 1, strike = x, dividendYield = 0, riskFreeRate = 0, maturity = mkt$t, volatility = 1)$value)


dotFoo <- function(...){
  dots.arg <- as.list(substitute(list(...)))[-1L]
  return(1)
}

pr.gl <- gaussLaguerrePricer(strikeMat = strikeMat, mkt = mkt, N = 64, alpha = 0, sigma.ref = 0.5, N.factors = 0, charFun = bscf, preCalc = NULL, sig = 1)

### test cos transform pricer with some Heston ###
library(affineModelR)
data("heston.params")

parListHestonLev2 <- parListHestonLev
parListHestonLev2$Q[["2"]] <- list(kpp = 3, eta = 1, phi = 0.1, rho = -0.9, lmb = 1)
parListHestonLev2$Q$jmp$lprop <- c(0,0)

parListHestonLev3 <- parListHestonLev
parListHestonLev3$Q[["2"]] <- list(kpp = 3, eta = 1, phi = 0.1, rho = -0.9, lmb = 1)
parListHestonLev3$Q$jmp<- list(lvec = 1, lprop = c(3,3), rhoc = -0.1, muYc = -0.02, muSc = 0.03, sigmaYc = 0.03)


stateVec <- matrix(c(3.7,1,1.3), nrow = 3)
library(abind)
strikeMatHeston <- abind(0.1*strikeMat, 0.15 * strikeMat, 0.2 * strikeMat)
strikeMatHeston <- array(strikeMatHeston, c(1,2,3))

pr.he.gl <- gaussLaguerrePricer(strikeMat = strikeMatHeston[,,1,drop=F], mkt = mkt, N = 64, alpha = 0, sigma.ref = NULL, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec[1,1,drop=F])

pr.he.cos <- cosTransformPricer(strikeMat = strikeMatHeston[,,1,drop=F], mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec[1,1,drop=F])

pr.he.gl.3 <- gaussLaguerrePricer(strikeMat = strikeMatHeston, mkt = mkt, N = 64, alpha = 0, sigma.ref = NULL, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec)

pr.he.cos.3 <- cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec)

stateVec2 <- cbind(c(3.7,1,1.3), c(1,1.5,0.2))

pr.he2.gl.3 <- gaussLaguerrePricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 64, alpha = 0, sigma.ref = NULL, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev2$Q, v.0 = stateVec2)

pr.he2.cos.3 <- cosTransformPricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev2$Q, v.0 = stateVec2)


pr.he3.gl.3 <- gaussLaguerrePricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 64, alpha = 0, sigma.ref = NULL, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev3$Q, v.0 = stateVec2)

pr.he3.cos.3 <- cosTransformPricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev3$Q, v.0 = stateVec2)
