PEB_df <- df
h_Z <- with(PEB_df, cbind(1, Z, X))
h_W <- with(PEB_df, cbind(1, W, X))
hfit <- hbridge(eta = eta, time = time, event = event, stime = stime, ctime = ctime, Z = h_Z, W = h_W, weights = weights)
PEB_est <- PEB_estimator(hfit$B_int, weights = weights)