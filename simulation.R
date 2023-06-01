rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Rcpp)
library(RcppArmadillo)
source("data_generating.R")
source("proxicausal.R")

Rcpp::sourceCpp("est.cpp")

eta <- 0.5
N <- 10000000
df <- data_gen(N, para_set)
truth <- mean(df$TT > eta)


N <- 500
rep = 1
rep_num <- 1000
B <- 1

PEB_est_result <- c()
PCB_est_result <- c()
PDR_est_result <- c()
DR_est_result <- c()

PEB_sd_result <- c()
PCB_sd_result <- c()
PDR_sd_result <- c()
DR_sd_result <- c()

PEB_covering <- c()
PCB_covering <- c()
PDR_covering <- c()
DR_covering <- c()
while(rep <= rep_num) {
  df <- data_gen(N, para_set)
  time <- df$CTT
  event <- df$Delta
  stime <- sort(unique(time[event == 1]))
  ctime <- sort(unique(time[event == 0]))
  K <- length(stime)
  K_c <- length(ctime)
  tau <- df$CTT > eta | (df$CTT < eta & df$Delta == 1)
  boot_PEB_est_result <- c()
  boot_PCB_est_result <- c()
  boot_PDR_est_result <- c()
  boot_DR_est_result <- c()
  for (b in 1:B) {
    weights <- rexp(N)
    weights <- weights/mean(weights)
    source("simu_PEB.R")
    boot_PEB_est_result <- c(boot_PEB_est_result, PEB_est)
    
    source("simu_PCB.R")
    boot_PCB_est_result <- c(boot_PCB_est_result, PCB_est)
 
    source("simu_PDR.R")
    boot_PDR_est_result <- c(boot_PDR_est_result, PDR_est)
    
    source("simu_DR.R")
    boot_DR_est_result <- c(boot_DR_est_result, DR_est)
  }
  PEB_sd_result <- c(PEB_sd_result, sd(boot_PEB_est_result))
  PCB_sd_result <- c(PCB_sd_result, sd(boot_PCB_est_result))
  PDR_sd_result <- c(PDR_sd_result, sd(boot_PDR_est_result))
  DR_sd_result <- c(DR_sd_result, sd(boot_DR_est_result))
  
  
  weights <- rep(1, N)
  
  source("simu_PEB.R")
  PEB_est_result <- c(PEB_est_result, PEB_est)

  source("simu_PCB.R")
  PCB_est_result <- c(PCB_est_result, PCB_est)

  source("simu_PDR.R")
  PDR_est_result <- c(PDR_est_result, PDR_est)
  
  source("simu_DR.R")
  DR_est_result <- c(DR_est_result, DR_est)


  PEB_covering <- c(PEB_covering, truth <= PEB_est + 1.96 * sd(boot_PEB_est_result) & truth >= PEB_est - 1.96 * sd(boot_PEB_est_result))
  PCB_covering <- c(PCB_covering, truth <= PCB_est + 1.96 * sd(boot_PCB_est_result) & truth >= PCB_est - 1.96 * sd(boot_PCB_est_result))
  PDR_covering <- c(PDR_covering, truth <= PDR_est + 1.96 * sd(boot_PDR_est_result) & truth >= PDR_est - 1.96 * sd(boot_PDR_est_result))
  DR_covering <- c(DR_covering, truth <= DR_est + 1.96 * sd(boot_DR_est_result) & truth >= DR_est - 1.96 * sd(boot_DR_est_result))
  
  print(rep)
  rep = rep + 1
}




 
running_result <- list(PEB_est_result = PEB_est_result,
                       PCB_est_result = PCB_est_result,
                       PDR_est_result = PDR_est_result,
                        DR_est_result = DR_est_result,
                       PEB_sd_result = PEB_sd_result,
                       PCB_sd_result = PCB_sd_result,
                       PDR_sd_result = PDR_sd_result,
                        DR_sd_result = DR_sd_result,
                       PEB_covering = PEB_covering,
                       PCB_covering = PCB_covering,
                       PDR_covering = PDR_covering,
                        DR_covering = DR_covering
)
 
stats_result <- with(running_result, 
                     list(
                     PEB_bias = mean(PEB_est_result) - truth,
                      PCB_bias = mean(PCB_est_result) - truth,
                      PDR_bias = mean(PDR_est_result) - truth, 
                     DR_bias = mean(DR_est_result) - truth,
                     PEB_see = sd(PEB_est_result),
                      PCB_see = sd(PCB_est_result),
                      PDR_see = sd(PDR_est_result), 
                     DR_see = sd(DR_est_result),
                     PEB_sd = mean(PEB_sd_result),
                      PCB_sd = mean(PCB_sd_result),
                      PDR_sd = mean(PDR_sd_result), 
                     DR_sd = mean(DR_sd_result),
                     PEB_covering = mean(PEB_covering),
                      PCB_covering = mean(PCB_covering),
                      PDR_covering = mean(PDR_covering),
                     DR_covering = mean(DR_covering)
                     ))
stats_result

