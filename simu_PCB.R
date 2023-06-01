PCB_df <- df
q_Z <- with(PCB_df, cbind(1, Z, X))
q_W <- with(PCB_df, cbind(1, W, X))
qfit <- qbridge(time = time, event = event, stime = stime, ctime = ctime, Z = q_Z, W = q_W, weights = weights)
PCB_est <- PCB_estimator(eta, time, tau, qfit$A_pred_int, weights = weights)
