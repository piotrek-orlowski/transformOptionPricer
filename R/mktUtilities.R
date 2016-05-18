#' Utilities for creating mkt lists for simulating option prices
#' @description List with fixed expiry (changing maturity day by day)
#' @rdname makeFixedExpiryMktList
#' @param num.days Length of list, how many days.
#' @param t vector of starting maturities
#' @param expiry.switch how often a new series of options is issued
#' @param r interest rate
#' @param q dividend yield
#' @param p stock price
#' @export
#' @return A list mkt lists

makeFixedExpiryMktList <- function(num.days, t = 37/365, expiry.switch = 30, r = 0, p = 1, q = 0){
  base.maturities <- t
  # we will start at the "issue" of the smallest basic maturity date
  issue.mat <- min(base.maturities)
  mkt.list <- list()
  for(ind in 1:num.days){
    t.loc <- (t*365 - (ind %% expiry.switch))/365
    mkt.list[[ind]] <- data.frame(r = r, p = p, q = q, day = ind-1, t = t.loc)
  }
  return(mkt.list)
}

#' @description List with fixed maturity each day
#' @rdname makeFixedExpiryMktList
#' @param num.days Length of list, how many days.
#' @param t vector of starting maturities
#' @param r interest rate
#' @param q dividend yield
#' @param p stock price
#' @export
#' @return A list containing lists svFast,svSlow and jmp

makeFixedMaturityMktList <- function(num.days, r = 0, t = 30/365, q = 0, p = 1){
  simple.mkt.list <- alply(.data = 1:num.days, .margins = 1, .fun = function(aa){return(list(day = aa-1, p = 1, q = q, r = r, t = t))})
  return(simple.mkt.list)
}