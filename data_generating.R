para_set <- list(mu_X = 1.5,
                 sigma_X = 0.35,
                 mu_U = 1.5,
                 sigma_U = 0.35,
                 mu_Z = c(0.2, 0.3, 0.8),
                 sigma_Z = 0.25,
                 mu_W = c(0.2, 0.5, 0.8),
                 sigma_W = 0.25,
                 mu_C = c(0.1, 0.35, 0.9),
                 mu_TT = c(0.25, 0.3, 0.7)
)



data_gen <- function(N, para_set) {
  # generate X0, U0
  X <- pmax(0, para_set$mu_X + rnorm(N, 0, para_set$sigma_X))
  U <- pmax(0, para_set$mu_U + rnorm(N, 0, para_set$sigma_U))

  
  # generate Z0
  Z <- cbind(1, X, U) %*% para_set$mu_Z + rnorm(N, 0, para_set$sigma_Z)
  
  # generate W0
  W <- cbind(1, X, U) %*% para_set$mu_W + rnorm(N, 0, para_set$sigma_W)

  C <- rexp(N, cbind(1, X, U) %*% para_set$mu_C)
  C <- pmin(C, 3)
  #generate Y
  TT <- rexp(N, cbind(1, X, U) %*% para_set$mu_TT)
  

  df <- data.frame(X, U, Z, W, C, TT, CTT = pmin(TT, C), Delta = TT <= C)
  return(df)
}
