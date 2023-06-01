DR_df <- df
#PDR_est <- PDR_estimator(K, dN_T = hfit$dN_T, risk_s = hfit$risk_s, B_int = hfit$B_int, dB_int = hfit$dB_int, A_pred_int = qfit$A_pred_int)

h_Z <- with(DR_df, cbind(1, W, Z, X))
h_W <- with(DR_df, cbind(1, W, Z, X))
hfit <- hbridge(eta = eta, time = time, event = event, stime = stime, ctime = ctime, Z = h_Z, W = h_W, weights = weights)

q_Z <- with(DR_df, cbind(1, W, Z, X))
q_W <- with(DR_df, cbind(1, W, Z, X))
qfit <- qbridge(time = time, event = event, stime = stime, ctime = ctime, Z = q_Z, W = q_W, weights = weights)


DR_est <- PDR_estimator(K_c, dN_C = qfit$dN_C, risk_c = qfit$risk_c, 
                        A_int = qfit$A_int, dA_int = qfit$dA_int, 
                        A_pred_int = qfit$A_pred_int, B_pred_int = hfit$B_pred_int, 
                        weights = weights)
