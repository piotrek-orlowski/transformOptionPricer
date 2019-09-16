#' Option pricing and Greeks in affine models
#'
#' @param strike_mat T x K x 1 array
#' @param mkt T x 1 data.frame
#' @param type call, put or otm options (character)
#' @param N integer num integration point see Fang and Oosterle 2008
#' @param intLim integration limits see Fang and Oosterle 2008
#' @param payCoeffFoo ibid
#' @param params Q parameter structure as in jumpDiffusionODEs in affineModelR
#' @param N.factors integer, number of factors in SV model
#' @param v.0 matrix of states in pricing, 1 x N
#' @param S.0 numeric, stock price at which options are evaluated
#' @param ode_sol array N x T x (1 + N.factors) of affine CF coefficients
#' @param ... further arguments to affineCF
#'
#' @return
#' @export
#'
#' @details This works for a single time period!
#' 
option_price_with_greeks <- function(strike_mat
                                     , mkt
                                     , type = c("call", "put", "otm")
                                     , N = 192L
                                     , intLim
                                     , payCoeffFoo = callCosCoeffs
                                     , params # Q-measure list for affineModelR
                                     , N.factors = 1L
                                     , v.0
                                     , S.0
                                     , ode_sol = NULL
                                     , ...){
  # Mon Sep 16 10:39:53 2019 ------------------------------
  # pricing problem dimensions
  T <- length(mkt$t)
  stopifnot(length(mkt$t) == dim(strike_mat)[1])
  S <- 9
  K <- dim(strike_mat)[2]
  
  # This is for affine models only Tue Jun 04 13:17:59 2019 ------------------------------
  charFun <- affineCF
  
  # set dv and ds for derivatives Tue Jun 04 13:25:37 2019 ------------------------------
  dv <- 1e-3 * as.numeric(v.0) / N.factors
  ds <- 1e-4 * S.0
  
  v_fwd <- v.0 + dv
  v_bwd <- v.0 - dv
  
  s_fwd <- S.0 + ds
  s_bwd <- S.0 - ds
  
  # expand strike_mat for the calculation of Greeks Tue Jun 04 13:18:08 2019 ------------------------------
  # strike_mat is T x K x 1
  # In order to calculate the Greeks, we will have to expand it to T x K x 9 and adjust strikes in a number of positions
  # strike_mat[,,1]: pricing at S.0, v.0
  # strike_mat[,,2]: pricing at S.0, v.0 - dv
  # strike_mat[,,3]: pricing at S.0, v.0 + dv
  # strike_mat[,,4]: pricing at S.0 + ds, v.0
  # strike_mat[,,5]: pricing at S.0 - ds, v.0
  # strike_mat[,,6]: pricing at S.0 - ds, v.0 - dv
  # strike_mat[,,7]: pricing at S.0 + ds, v.0 + dv
  # strike_mat[,,8]: pricing at S.0 - ds, v.0 + dv
  # strike_mat[,,9]: pricing at S.0 + ds, v.0 - dv
  
  # Initialise array of zeros
  strike_mat_expanded <- array(0, c(dim(strike_mat)[c(1L,2L)], 9L))
  
  # Everything is based on the initial strikes, so fill the array with them.
  # For slices (1,2,3) there will be no change to the strike_mat
  strike_mat_expanded <- lapply(1:9, function(x) strike_mat)
  strike_mat_expanded <- do.call(abind, strike_mat_expanded)
  
  # For slices (4,7,9) the stock price is increased by ds
  strike_mat_expanded[,,c(4L,7L,9L)] <- log(exp(strike_mat_expanded[,,c(4L,7L,9L)]) * S.0 / s_fwd)
  
  # For slicess (5,6,8) the stock price is decreased by ds
  strike_mat_expanded[,,c(5L,6L,8L)] <- log(exp(strike_mat_expanded[,,c(5L,6L,8L)]) * S.0 / s_bwd)
  
  # expand v.0 for the calculation of Greeks Tue Jun 04 13:18:40 2019 ------------------------------
  v.0_expanded <- do.call(rbind, list(v.0, v_bwd, v_fwd, v.0, v.0, v_bwd, v_fwd, v_fwd, v_bwd))
  S.0_expanded <- c(S.0, S.0, S.0, s_fwd, s_bwd, s_bwd, s_fwd, s_bwd, s_fwd)
  
  # Mon Sep 16 10:30:47 2019 ------------------------------
  # prepare a pre_calc array for cosTransformPricer
  if(is.null(ode_sol)){
    pre_calc <- NULL 
  } else {
    pre_calc <- affineCFevalCpp(coeffs = ode_sol
                                , stateMat = v.0_expanded
                                , retLog = FALSE)
    
    pre_calc <- array(pre_calc, dim = c(N,T,S,K))
    pre_calc <- aperm(pre_calc, perm = c(1,2,4,3))
  }
  
  
  # calculate all option prices Tue Jun 04 13:34:32 2019 ------------------------------
  price_array <- cosTransformPricer(strikeMat = strike_mat_expanded
                                    , mkt = mkt
                                    , N = N
                                    , intLim = intLim
                                    , payCoeffFoo = payCoeffFoo
                                    , N.factors = N.factors
                                    , charFun = charFun
                                    , params.Q = params
                                    , v.0 = v.0_expanded
                                    , preCalc = pre_calc
                                    , ...)[[type[1L]]]
  
  # the array contains relative option prices, they have to be multiplied by the appropriate stock prices Tue Jun 04 13:40:14 2019 ------------------------------
  price_array <- lapply(1L:9L, function(ind){
    price_array[,,ind,drop=FALSE] * S.0_expanded[ind]
  })
  
  price_array <- do.call(abind, price_array)
  
  # all calculations follow the Finite Difference article on Wikipedia as of timestamp Tue Jun 04 13:50:33 2019 ------------------------------
  
  # delta calculation: central difference Tue Jun 04 13:45:16 2019 ------------------------------
  delta <- price_array[,,4L,drop=FALSE] - price_array[,,5L,drop=FALSE]
  delta <- delta / (2.0 * ds)
  
  # vega calculation: central difference Tue Jun 04 13:45:31 2019 ------------------------------
  vega <- price_array[,,3L,drop=FALSE] - price_array[,,2L,drop=FALSE]
  vega <- vega / (2.0 * sum(dv))
  
  # gamma calculation Tue Jun 04 13:45:39 2019 ------------------------------
  gamma <- price_array[,,4L,drop=FALSE] + price_array[,,5L,drop=FALSE] - 2.0 * price_array[,,1L,drop=FALSE]
  gamma <- gamma / (ds^2.0)
  
  # volga calculation Tue Jun 04 13:45:48 2019 ------------------------------
  volga <- price_array[,,3L,drop=FALSE] + price_array[,,2L,drop=FALSE] - 2.0 * price_array[,,1L,drop=FALSE]
  volga <- volga / (sum(dv)^2.0)
  
  # vanna calculation Tue Jun 04 13:45:56 2019 ------------------------------
  vanna <- price_array[,,7L,drop=FALSE] - price_array[,,4L,drop=FALSE] - price_array[,,3L,drop=FALSE] + 2.0 * price_array[,,1L,drop=FALSE] - price_array[,,5L,drop=FALSE] - price_array[,,2L,drop=FALSE] + price_array[,,6L,drop=FALSE]
  vanna <- vanna / (2.0 * sum(dv) * ds)
  
  return(list(option_price = price_array[,,1L,drop=FALSE]
              , delta = delta
              , vega = vega
              , gamma = gamma
              , volga = volga
              , vanna = vanna))
}