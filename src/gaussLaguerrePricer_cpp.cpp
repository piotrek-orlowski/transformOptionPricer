#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//'@export
// [[Rcpp::export]]
arma::cube glPricer_cpp(const arma::cube& strikeMat, const arma::mat& mkt, const arma::vec& glWts, const arma::vec& glNodes, arma::cx_cube& cfVals, int Nfactors, double alpha = 0, double sigmaRef = 1.0){
  
  // define 1i imaginary
  std::complex<double> i1(0,1);
  
  // capture dimensions
  arma::uvec strikeMat_dim(3);
  strikeMat_dim(0) = strikeMat.n_rows;
  strikeMat_dim(1) = strikeMat.n_cols;
  strikeMat_dim(2) = strikeMat.n_slices;
  
  arma::uvec cfVals_dim(3);
  cfVals_dim(0) = cfVals.n_rows;
  cfVals_dim(1) = cfVals.n_cols;
  cfVals_dim(2) = cfVals.n_slices;
  
  // add discounting to characteristic function coefficients
  arma::cx_cube driftMat(cfVals_dim(0), cfVals_dim(1), cfVals_dim(2));
  driftMat.fill(0.0);
  arma::mat dscnting = (mkt.col(0) % (mkt.col(1)-mkt.col(2))).t();
  for(int ss = 0; ss < cfVals_dim(2); ss++){
    cfVals.slice(ss) = cfVals.slice(ss) % arma::exp(i1 * (glNodes * dscnting));
  }
  
  arma::cx_cube cfVals_reshaped(cfVals_dim(1), cfVals_dim(2), cfVals_dim(0));
  for(unsigned int dim1=0; dim1 < cfVals_dim(1); dim1++){
    for(unsigned int dim2=0; dim2 < cfVals_dim(2); dim2++){
      for(unsigned int dim0=0; dim0 < cfVals_dim(0); dim0++){
        cfVals_reshaped(dim1, dim2, dim0) = cfVals(dim0, dim1, dim2);
      }
    }
  }
  
  int T = strikeMat_dim(0);
  int K = strikeMat_dim(1);
  int S;
  if(strikeMat_dim.size() > 2L){
    S = strikeMat_dim(2); 
  } else {
    S = 1L;
  }
  int N = cfVals_dim(0);
  
  // log-forward prices from mkt
  arma::vec logFvec = mkt.col(0) % (mkt.col(1) - mkt.col(2));
  arma::vec Fvec = arma::exp(logFvec);
  
  // calculate field with values of the BS characteristic function, cfBS
  arma::mat M1(T,S,arma::fill::zeros);
  arma::mat M2(T,S,arma::fill::zeros);
  
  M1 = logFvec * arma::ones<arma::rowvec>(S);
  M2 = mkt.col(0) * (sigmaRef * arma::ones<arma::rowvec>(S));
  
  arma::mat M3_filler = arma::repmat(arma::reshape(2.0*M1-M2,T*S,1),N,1) % arma::reshape(arma::repmat(glNodes,1,T*S).t(),T*S*N,1);
  const arma::cube M3(M3_filler.memptr(), T, S, N, false, true);
  
  arma::mat M4_filler = arma::repmat(arma::reshape(-M2,T*S,1),N,1) % arma::reshape(arma::repmat(arma::pow(glNodes,2.0),1,T*S).t(),T*S*N,1);
  const arma::cube M4(M4_filler.memptr(),T,S,N, false, true);
  
  arma::cx_cube cfBStemp = arma::exp(0.5* i1 * M3 + 0.5*M4);
  arma::field<arma::cx_cube> cfBS(K);
  for(int kk = 0; kk < K; kk++){
    cfBS(kk) = cfBStemp;
  }
  
  // Calculate BS prices at reference volatility
  arma::mat sLittle_filler = arma::repmat(arma::reshape(M2,T*S,1),K,1); // s <- aperm(array(rep(as.vector(M2),K),c(T,S,K)),c(1,3,2))
  arma::cube sLittle(sLittle_filler.memptr(),T,K,S,false,true);
  
  arma::mat mLittle_filler = repmat(logFvec,S*K,1);
  arma::cube mLittle(mLittle_filler.memptr(),T,K,S,false,true);
  
  arma::cube d1 = arma::pow(sLittle,-0.5) % (mLittle - strikeMat + 0.5*sLittle);
  arma::cube d2 = d1 - arma::pow(sLittle,0.5);
  
  arma::mat df_filler = arma::repmat(arma::exp(- mkt.col(0) % mkt.col(1)),K*S,1);
  arma::cube df(df_filler.memptr(),T,K,S,false,true);
  
  arma::mat Fcube_filler = arma::repmat(Fvec,K*S,1);
  arma::cube Fcube(Fcube_filler.memptr(),T,K,S,false,true);
  
  arma::cube ones_cube(T,K,S,arma::fill::ones);
  
  arma::cube pBS = df % (Fcube % (0.5*(ones_cube + arma::erf(d1/sqrt(2.0)))) - arma::exp(strikeMat) % (0.5 * (ones_cube + arma::erf(d2/sqrt(2.0)))));
  
  // create field to hold affine CF
  arma::field<arma::cx_cube> cf(K);
  for(int kk = 0; kk < K; kk++){
    cf(kk) = cfVals_reshaped;
  }

  // difference of CFs
  arma::mat uFiller_alpha = arma::reshape(arma::repmat(arma::pow(glNodes, 1.0 + alpha).t(),T*S,1),T*S*N,1);
  arma::mat uFiller = arma::reshape(arma::repmat(glNodes.t(),T*S,1),T*S*N,1);
  arma::cube uScale_alpha(uFiller_alpha.memptr(),N,T,S,true,true);
  arma::cube uScale_alpha_reshaped = arma::reshape(uScale_alpha,T,S,N);
  arma::cube uScale(uFiller.memptr(),N,T,S,true,true);
  arma::cube uScale_reshaped = arma::reshape(uScale,T,S,N);
  
  arma::field<arma::cx_cube> dCF(K);
  for(int kk = 0; kk < K; kk++){
    dCF(kk) = (cf(kk) - cfBS(kk)) / ((i1*uScale_alpha_reshaped) % (arma::ones(T,S,N) - i1*uScale_reshaped));
  }
  
  arma::field<arma::cx_cube> cExp(K);
  for(int kk = 0; kk < K; kk++){
    arma::mat locStrikeTube = strikeMat.tube(0,kk,T-1,kk);
    arma::mat locStrikeMat = arma::repmat(arma::reshape(locStrikeTube,T*S,1),N,1);
    arma::cube locStrikeCube(locStrikeMat.memptr(),T,S,N);
    cExp(kk) = -i1 * locStrikeCube;
    cExp(kk) = cExp(kk) % uScale_reshaped;
    cExp(kk) = arma::exp(cExp(kk));
  }
  
  arma::mat wnFiller = arma::reshape(arma::repmat(glWts.t() % arma::exp(glNodes).t(),T*S,1),T*S*N,1);
  arma::cube wn(wnFiller.memptr(),N,T,S,true,true);
  arma::cube wn_reshaped = arma::reshape(wn,T,S,N);
  
  arma::field<arma::cx_cube> dC(K);
  for(int kk = 0; kk < K; kk++){
    dC(kk) = dCF(kk) % cExp(kk) % wn_reshaped;
  }
  
  arma::mat df1, df2;
  
  arma::field<arma::mat> call(K), put(K);
  for(int kk = 0; kk < K; kk++){
    arma::cx_mat tmp(T,S,arma::fill::zeros);
    for(int nn = 0; nn < N; nn++){
      tmp += dC(kk).slice(nn);
    }
    df1 = arma::repmat(arma::exp(- mkt.col(0) % mkt.col(1)),1,S);
    df2 = arma::repmat(arma::exp(- mkt.col(0) % mkt.col(2)),1,S);
    arma::mat locStrikeMat = strikeMat.tube(0,kk,T-1,kk);
    arma::mat locPBS = pBS.tube(0,kk,T-1,kk);
    call(kk) = - arma::exp(locStrikeMat) % arma::real(tmp)/arma::datum::pi % df1 + locPBS;
    put(kk) = call(kk) + arma::exp(locStrikeMat) % df1 - df2; // call + exp(strikeMat)*df1 - 1*df2
  }
  
  arma::cube otm(T,K,S,arma::fill::zeros);
  for(int kk = 0; kk < K; kk++){
    arma::mat locStrikeMat = strikeMat.tube(0,kk,T-1,kk);
    arma::mat locForwardMat(locStrikeMat.n_rows, locStrikeMat.n_cols);
    for(int tt = 0; tt < T; tt++){
      locForwardMat.row(tt).fill(logFvec(tt));
    }
    arma::umat strikeLessThanForward = arma::find(locStrikeMat <= locForwardMat);
    arma::mat tmpPut = put(kk);
    arma::mat tmp = call(kk);
    tmp(strikeLessThanForward) = tmpPut(strikeLessThanForward);
    otm.tube(0,kk,T-1,kk) = tmp;
  }
  
  return otm;
  // return List::create(Named("otm") = otm, Named("call") = call, Named("put") = put);
}