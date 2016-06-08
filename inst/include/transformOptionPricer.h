#ifndef _GL_PRICER_H_
#define _GL_PRICER_H_
#include <RcppArmadillo.h>

arma::cube glPricer_cpp(const arma::cube& strikeMat, const arma::mat& mkt, const arma::vec& glWts, const arma::vec& glNodes, const arma::cx_cube& cfVals, int Nfactors, double alpha = 0, double sigmaRef = 1.0);

#endif