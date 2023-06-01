IPCW_df <- df


IPCW <- with(IPCW_df, exp(cbind(1, X, U0, U1) %*% para_set$mu_C * pmin(time, eta)))

IPCW_est <- IPCW_estimator(eta, time, tau, IPCW, weights = weights)
