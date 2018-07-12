### test cos transform pricer with BS ###
library(transformOptionPricer)

bscf <- function(u,sig,t.vec,N.factors,mkt){
  tau <- mkt$t
  return(exp(-1/2*sig^2*u*tau + 1/2*sig^2*u^2*tau))
}

library(RQuantLib)


strikeMat <- t(matrix(log(c(0.9,1.0,1.1)),3,1))

mkt <- data.frame(p=0,q=0.0,r=0.0,t=c(1/2))

sig <- 0.3

intLim <- mkt$r * mkt$t + c(-10,10) * sqrt(sig^2*mkt$t)

pr <- cosTransformPricer(strikeMat = strikeMat, mkt = mkt, N = 2e5, intLim = c(-5,5), payCoeffFoo = callCosCoeffs, N.factors = 0, charFun = bscf, sig = sig)
pr <- cosTransformPricer(strikeMat = strikeMat, mkt = mkt, N = 2e5, intLim = intLim, payCoeffFoo = callCosCoeffs, N.factors = 0, charFun = bscf, sig = sig)

pr.true <- list()
pr.true$call <- sapply(exp(as.numeric(strikeMat)), function(x) EuropeanOption(type = 'call', underlying = 1, strike = x, dividendYield = 0, riskFreeRate = 0, maturity = mkt$t, volatility = sig)$value)
pr.true$put <- sapply(exp(as.numeric(strikeMat)), function(x) EuropeanOption(type = 'put', underlying = 1, strike = x, dividendYield = 0, riskFreeRate = 0, maturity = mkt$t, volatility = sig)$value)


dotFoo <- function(...){
  dots.arg <- as.list(substitute(list(...)))[-1L]
  return(1)
}

pr.gl <- gaussLaguerrePricer(strikeMat = strikeMat, mkt = mkt, N = 64, alpha = 0, sigma.ref = 0.5, N.factors = 0, charFun = bscf, preCalc = NULL, sig = sig)

### test cos transform pricer with some Heston ###
library(affineModelR)
data("heston-parameters")

stateVec <- matrix(c(3.7,1,0.15), nrow = 3)

cumulants <- affineCF(u = cbind(c(1,2,4),0), params.Q = parListHeston$Q, t.vec = c(0.04,0.25,0.5,1), v.0 = stateVec, N.factors = 1)
cumulants <- log(Re(cumulants))

intLimsLeft <- cumulants[1,,] - 10*sqrt(cumulants[2,,] + sqrt(cumulants[3,,]))
intLimsRight <- cumulants[1,,] + 10*sqrt(cumulants[2,,] + sqrt(cumulants[3,,]))

library(abind)
strikeMatHeston <- seq(-0.3,0.3,by=0.025)
strikeMatHeston <- array(strikeMatHeston, c(1,length(strikeMatHeston),1))

mkt$q <- 0.01
mkt$r <- 0.03

N_sequence <- seq(32,1024,by=64)

pr_cos_1 <- sapply(N_sequence, function(nn) cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = nn, intLim = c(intLimsLeft[3,1], intLimsRight[3,1]), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec[1,1,drop=F])$otm)
pr_cos_2 <- sapply(N_sequence, function(nn) cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = nn, intLim = c(intLimsLeft[3,2], intLimsRight[3,2]), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec[2,1,drop=F])$otm)
pr_cos_3 <- sapply(N_sequence, function(nn) cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = nn, intLim = c(intLimsLeft[3,3], intLimsRight[3,3]), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec[3,1,drop=F])$otm)

plot(log(pr_cos_1[,1]), ylim = range(log(pr_cos_1[pr_cos_1 > 0])))
points(log(pr_cos_1[,4]), col = 'blue', pch = 4)
points(log(pr_cos_1[,8]), col = 'brown', pch = 4)
points(log(pr_cos_1[,12]), col = 'darkgreen', pch = 4)
points(log(pr_cos_1[,16]), col = 'orange', pch = 4)

plot(log(pr_cos_2[,1]), ylim = range(log(pr_cos_2[pr_cos_2 > 0])))
points(log(pr_cos_2[,4]), col = 'blue', pch = 4)
points(log(pr_cos_2[,8]), col = 'brown', pch = 4)
points(log(pr_cos_2[,12]), col = 'darkgreen', pch = 4)
points(log(pr_cos_2[,16]), col = 'orange', pch = 4)

plot(log(pr_cos_3[,1]), ylim = range(log(pr_cos_3[pr_cos_3 > 0])))
points(log(pr_cos_3[,4]), col = 'blue', pch = 4)
points(log(pr_cos_3[,8]), col = 'brown', pch = 4)
points(log(pr_cos_3[,12]), col = 'darkgreen', pch = 4)
points(log(pr_cos_3[,16]), col = 'orange', pch = 4)

#### HEston with Leverage ####

parListHestonLev <- parListHeston
parListHestonLev$Q$`1`$rho <- -0.95

stateVec <- matrix(c(3.7,1,0.15), nrow = 3)

cumulants <- affineCF(u = cbind(c(1,2,4,6),0), params.Q = parListHestonLev$Q, t.vec = c(0.04,0.25,0.5,1), v.0 = stateVec, N.factors = 1)
cumulants <- log(Re(cumulants))

intLimsLeft <- cumulants[1,,] - 10*sqrt(cumulants[2,,] + sqrt(cumulants[3,,] + sqrt(cumulants[4,,])))
intLimsRight <- cumulants[1,,] + 10*sqrt(cumulants[2,,] + sqrt(cumulants[3,,] + sqrt(cumulants[4,,])))

library(abind)
strikeMatHeston <- seq(-0.3,0.3,by=0.025)
strikeMatHeston <- array(strikeMatHeston, c(1,length(strikeMatHeston),1))

mkt$q <- 0.01
mkt$r <- 0.03

mkt$t <- 0.04

N_sequence <- seq(64,1092,by=64)

pr_cos_1 <- sapply(N_sequence, function(nn) cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = nn, intLim = c(intLimsLeft[1,1], intLimsRight[1,1]), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHestonLev$Q, v.0 = stateVec[1,1,drop=F])$otm)
pr_cos_2 <- sapply(N_sequence, function(nn) cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = nn, intLim = c(intLimsLeft[1,2], intLimsRight[1,2]), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHestonLev$Q, v.0 = stateVec[2,1,drop=F])$otm)
pr_cos_3 <- sapply(N_sequence, function(nn) cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = nn, intLim = c(intLimsLeft[1,3], intLimsRight[1,3]), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHestonLev$Q, v.0 = stateVec[3,1,drop=F])$otm)

plot(log(pr_cos_1[,1]), ylim = range(log(pr_cos_1[pr_cos_1 > 0])))
points(log(pr_cos_1[,4]), col = 'blue', pch = 4)
points(log(pr_cos_1[,8]), col = 'brown', pch = 4)
points(log(pr_cos_1[,12]), col = 'darkgreen', pch = 4)
points(log(pr_cos_1[,16]), col = 'orange', pch = 4)

plot(log(pr_cos_2[,1]), ylim = range(log(pr_cos_2[pr_cos_2 > 0])))
points(log(pr_cos_2[,4]), col = 'blue', pch = 4)
points(log(pr_cos_2[,8]), col = 'brown', pch = 4)
points(log(pr_cos_2[,12]), col = 'darkgreen', pch = 4)
points(log(pr_cos_2[,16]), col = 'orange', pch = 4)

plot(log(pr_cos_3[,1]), ylim = range(log(pr_cos_3[pr_cos_3 > 0])))
points(log(pr_cos_3[,4]), col = 'blue', pch = 4)
points(log(pr_cos_3[,8]), col = 'brown', pch = 4)
points(log(pr_cos_3[,12]), col = 'darkgreen', pch = 4)
points(log(pr_cos_3[,16]), col = 'orange', pch = 4)

#### Legacy ####

pr.he.gl <- gaussLaguerrePricer(strikeMat = strikeMatHeston[,,1,drop=F], mkt = mkt, N = 80, alpha = 0, sigma.ref = NULL, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec[1,1,drop=F])

pr.he.cos <- cosTransformPricer(strikeMat = strikeMatHeston[,,1,drop=F], mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec[1,1,drop=F])

pr.he.gl.3 <- gaussLaguerrePricer(strikeMat = strikeMatHeston, mkt = mkt, N = 64, alpha = 0, sigma.ref = NULL, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec)

pr.he.cos.3 <- cosTransformPricer(strikeMat = strikeMatHeston, mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 1, charFun = affineCF, params.Q = parListHeston$Q, v.0 = stateVec)

stateVec2 <- cbind(c(3.7,1,1.3), c(1,1.5,0.2))

pr.he2.gl.3 <- gaussLaguerrePricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 64, alpha = 0, sigma.ref = NULL, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev2$Q, v.0 = stateVec2)

pr.he2.cos.3 <- cosTransformPricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev2$Q, v.0 = stateVec2)


pr.he3.gl.3 <- gaussLaguerrePricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 64, alpha = 0, sigma.ref = NULL, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev3$Q, v.0 = stateVec2)

pr.he3.cos.3 <- cosTransformPricer(strikeMat = 3*strikeMatHeston, mkt = mkt, N = 220, intLim = c(-3,3), payCoeffFoo = callCosCoeffs, N.factors = 2, charFun = affineCF, params.Q = parListHestonLev3$Q, v.0 = stateVec2)


### test with put pricing ###

library(transformOptionPricer)

bscf <- function(u,sig,t.vec,N.factors,mkt){
  tau <- mkt$t
  return(exp(-1/2*sig^2*u*tau + 1/2*sig^2*u^2*tau))
}

library(RQuantLib)


strikeMat <- t(matrix(log(c(0.9,1.0,1.1)),3,1))

mkt <- data.frame(p=0,q=0.0,r=0.0,t=c(1/2))

sig <- 0.3

intLim <- mkt$r * mkt$t + c(-10,10) * sqrt(sig^2*mkt$t)

pr <- cosTransformPricer(strikeMat = strikeMat, mkt = mkt, N = 50, intLim = intLim, payCoeffFoo = putCosCoeffs, N.factors = 0, charFun = bscf, sig = sig)

pr.true <- list()
pr.true$call <- sapply(exp(as.numeric(strikeMat)), function(x) EuropeanOption(type = 'call', underlying = 1, strike = x, dividendYield = 0, riskFreeRate = 0, maturity = mkt$t, volatility = sig)$value)
pr.true$put <- sapply(exp(as.numeric(strikeMat)), function(x) EuropeanOption(type = 'put', underlying = 1, strike = x, dividendYield = 0, riskFreeRate = 0, maturity = mkt$t, volatility = sig)$value)


dotFoo <- function(...){
  dots.arg <- as.list(substitute(list(...)))[-1L]
  return(1)
}

pr.gl <- gaussLaguerrePricer(strikeMat = strikeMat, mkt = mkt, N = 64, alpha = 0, sigma.ref = 0.5, N.factors = 0, charFun = bscf, preCalc = NULL, sig = sig)