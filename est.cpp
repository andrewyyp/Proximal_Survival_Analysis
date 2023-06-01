# include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



// [[Rcpp::export]]
SEXP hbridge(double eta,
             NumericVector time,
             NumericVector event,
             NumericVector stime, 
             NumericVector ctime, 
             NumericMatrix Z,
             NumericMatrix W, 
             NumericVector weights){
  
  // variables for estimation
  int N = time.size(), K_s = stime.size(), K_c = ctime.size(), dim_Z = Z.cols(), dim_W = W.cols();
  NumericVector B_numer(dim_W);
  NumericMatrix B(K_s, dim_W), dB(K_s, dim_W), dB_int(N, K_s), B_int(N, K_s + 1), B_denom(dim_Z, dim_W), B_pred_int(N, K_c);
  NumericMatrix dN(N, K_s),  risk_s(N, K_s);
  double del = 0.001;
  
  
  for (int j = K_s - 1; j >= 0; j--) {
    for (int dim_w = 0; dim_w < dim_W; dim_w++) {
      B_numer[dim_w] = 0;
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        B_denom(dim_w, dim_z) = 0;
      }
    }
    for (int i = 0; i < N; i++) {
      // estimation
      dN(i, j) = (time[i] == stime[j]) ? 1:0;
      dN(i, j) *= event[i];
      risk_s(i, j) = (time[i] >= stime[j]) ? 1:0;
      if (stime[j] > eta) 
        continue;
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        B_numer[dim_z] += weights[i] * Z(i, dim_z) * exp(B_int(i, j + 1)) * dN(i, j);
        for (int dim_w = 0; dim_w < dim_W; dim_w++) {
          B_denom(dim_z, dim_w) += weights[i] * Z(i, dim_z) * W(i, dim_w) * risk_s(i, j) * exp(B_int(i, j + 1));
        }
      }
    }
    
    
    if (stime[j] <= eta) {
        // solve linear system using RcppArmadillo
      try
      {
        arma::mat X(B_denom.begin(), dim_W, dim_Z, false);
        arma::mat tmp_pinv = arma::pinv(X, del);
        arma::mat tmp_res = tmp_pinv * as<arma::vec>(B_numer);
        for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
          dB(j, dim_z) = tmp_res(dim_z, 0);
        }
      } catch (...) {
        
      }
      
      
      
      // estimation
      for (int i = 0; i < N; i++) {
        for (int dim_w = 0; dim_w < dim_W; dim_w++) {
          dB_int(i, j) += W(i, dim_w) * dB(j, dim_w);
        }
        B_int(i, j) = B_int(i, j + 1) - dB_int(i, j) ;
      }
      
      
      for (int dim_w = 0; dim_w < dim_W; dim_w++) {
        if (j == K_s - 1) {
          B(j, dim_w) = 0 - dB(j, dim_w);
        }
        else {
          B(j, dim_w) = B(j + 1, dim_w) - dB(j, dim_w);
        }
      }
    }
  }

  int j_s = 0;
  int j_c = 0;
  while(j_c < K_c) {
    if (j_s < K_s) {
      if (ctime[j_c] < stime[j_s + 1])
        for (int i = 0; i < N; i++) {
          B_pred_int(i, j_c) = B_int(i, j_s);
        }
        else {
          ++j_s;
          continue;
        }
    } else {
      for (int i = 0; i < N; i++) {
        B_pred_int(i, j_c) = B_int(i, j_s);
      }
    }
    ++j_c;
  }
  
  return List::create(
    Named("dN_T") = dN,
    Named("risk_s") = risk_s,
    Named("dB") = dB,
    Named("B") = B,
    Named("dB_int") = dB_int,
    Named("B_int") = B_int,
    Named("B_pred_int") = B_pred_int);
}


// [[Rcpp::export]]
SEXP qbridge(NumericVector time,
              NumericVector event,
              NumericVector stime, 
              NumericVector ctime, 
              NumericMatrix Z,
              NumericMatrix W, 
              NumericVector weights){
  
  // variables for estimation
  int N = time.size(), K_s = stime.size(), K_c = ctime.size(), dim_Z = Z.cols(), dim_W = W.cols();
  NumericVector A_numer(dim_W);
  NumericMatrix A(K_c, dim_W), dA(K_c, dim_W), A_int(N, K_c + 1), dA_int(N, K_c), A_denom(dim_Z, dim_W), A_pred_int(N, K_s);
  NumericMatrix dN(N, K_c),  risk_c(N, K_c);
  double del = 0.001;

  
  for (int j = 0; j < K_c; j++) {
    for (int dim_w = 0; dim_w < dim_W; dim_w++) {
      A_numer[dim_w] = 0;
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        A_denom(dim_w, dim_z) = 0;
      }
    }
    for (int i = 0; i < N; i++) {
      // estimation
      dN(i, j) = (time[i] == ctime[j]) ? 1:0;
      dN(i, j) *= 1 - event[i];
      risk_c(i, j) = (time[i] >= ctime[j]) ? 1:0;
      for (int dim_w = 0; dim_w < dim_W; dim_w++) {
        A_numer[dim_w] += weights[i] * W(i, dim_w) * exp(A_int(i, j)) * dN(i, j);
        for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
          A_denom(dim_w, dim_z) += weights[i] * W(i, dim_w) * Z(i, dim_z) * risk_c(i, j) * exp(A_int(i, j));
        }
      }
    }
    


    // solve linear system using RcppArmadillo
    try
    {
      arma::mat X(A_denom.begin(), dim_Z, dim_W, false);
      arma::mat tmp_pinv = arma::pinv(X, del);
      arma::mat tmp_res = tmp_pinv * as<arma::vec>(A_numer);
      for (int dim_w = 0; dim_w < dim_W; dim_w++) {
        dA(j, dim_w) = tmp_res(dim_w, 0);
      }
    } catch (...) {
      
    }
      
      

    
    // estimation
    for (int i = 0; i < N; i++) {
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        dA_int(i, j) += Z(i, dim_z) * dA(j, dim_z);
      }
      A_int(i, j + 1) += A_int(i, j) + dA_int(i, j);
    }
    

    for (int dim_w = 0; dim_w < dim_W; dim_w++) {
      if (j == 0) {
        A(j, dim_w) = 0 + dA(j, dim_w);
      }
      else{
        A(j, dim_w) = A(j - 1, dim_w) + dA(j, dim_w);
      }
    }
  }
  int j_s = 0;
  int j_c = 0;
  while(j_s < K_s) {
    if (j_c < K_c) {
      if (stime[j_s] < ctime[j_c + 1]) {
        for (int i = 0; i < N; i++) {
          A_pred_int(i, j_s) = A_int(i, j_c);
        }
      } else {
        ++j_c;
        continue;
      }
    } else {
      for (int i = 0; i < N; i++) {
        A_pred_int(i, j_s) = A_int(i, j_c);
      }
    }
    ++j_s;
  }
  
  
  return List::create(
    Named("dN_C") = dN,
    Named("risk_c") = risk_c,
    Named("dA") = dA,
    Named("A") = A,
    Named("A_int") = A_int,
    Named("dA_int") = dA_int,
    Named("A_numer") = A_numer,
    Named("A_denom") = A_denom,
    Named("A_pred_int") = A_pred_int);
}