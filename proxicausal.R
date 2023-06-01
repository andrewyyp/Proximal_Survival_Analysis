PEB_estimator <- function(B_int, weights) {
  return(mean(weights * exp(B_int[, 1])))
}


PCB_estimator <- function(eta, time, tau, A_pred_int, weights) {
  time_order <- order(time)
  j <- 1
  numer <- 0
  denom <- 0
  for (i in time_order) {
    if(tau[i] == 1) { 
      numer <- numer + weights[i] * exp(A_pred_int[i, j]) * (time[i] > eta)
      denom <- denom + weights[i] * exp(A_pred_int[i, j]) 
      if(time[i] < eta) {
        j <- j + 1
      }
    }
  }
  return(numer/denom)
}


# PDR_estimator <- function(K, dN_T, risk_s, B_int, dB_int, A_pred_int, weights) {
#   numer <- sum(weights * exp(B_int[, 1])) + sum(t(weights * risk_s * exp(A_pred_int[, 1:K]) * exp(B_int[, 1:K]) * (dN_T - dB_int)) * (stime <= eta))
#   denom <- sum(weights) + sum(t(weights * risk_s * exp(A_pred_int[, 1:K]) * dN_T) * (time <= eta))
#   return(numer/denom)
# }

PDR_estimator <- function(K_c, dN_C, risk_c, A_int, dA_int, A_pred_int, B_pred_int, weights) {
  time_order <- order(time)
  j <- 1
  numer <- 0
  denom <- 0
  for (i in time_order) {
    if(tau[i] == 1) {
      numer <- numer + weights[i] * exp(A_pred_int[i, j]) * (time[i] > eta)
      denom <- denom + weights[i] * exp(A_pred_int[i, j])
      if(time[i] < eta) {
        j <- j + 1
      }
    }
  }
  numer <- numer - sum(t(t(weights * risk_c * exp(B_pred_int[, 1:K_c]) * exp(A_int[, 1:K_c]) * (dA_int - dN_C)) * (ctime <= eta)))
  denom <- denom - sum(t(t(weights * risk_c * exp(A_int[, 1:K_c]) * (dA_int - dN_C)) * (ctime <= eta)))
  return(numer/denom)
}


IPCW_estimator <- function(eta, time, tau, IPCW, weights) {
  numer <- sum(weights * tau * IPCW * (time > eta))
  denom <- sum(weights * tau * IPCW)
  return(numer/denom)
}


transform_W <- function(x) {
  return(1 / (1 + exp(-3 * (x - 1.5))) - 0.5)
}

transform_Z <- function(x) {
  return(1 / (1 + exp(-3 * (x - 1.5))) - 0.5)
}
