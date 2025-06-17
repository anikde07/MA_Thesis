args       = commandArgs(trailingOnly = TRUE)
N          = as.numeric(args[1])
TM         = as.numeric(args[2])
MC         = 10000

print(paste0("KSS only. N = ", N, " MC = ", MC, " T = ", TM))

library("parallel")
library("phtt")
library("ffscb")
library("pspline")

st = Sys.time()

alpha_confint = 0.05
z_value       = qnorm(1 - (alpha_confint / 2))
beta_true     = -5
grid          = make_grid(TM, rangevals = c(0, 1))
mu            = meanf_poly(grid)
cov.m = make_cov_m(cov.f = ffscb::covf_st_matern, grid = grid, cov.f.params = c(3 / 2, 0.5)) 

## Parallelization cluster ##

num_cores = detectCores() 
cl        = makeCluster(num_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, c(library("phtt"), library("ffscb"), library("pspline")))
######################
dml <- parLapply(cl, 1:MC, function(v) {
  
  gc()

  FE_fcts_mat   = make_sample(mean.v = rep(0, TM), cov.m = cov.m, N = N, dist = "rnorm")
  eps_y_mat     = matrix(rnorm(N * TM, 0, 0.2), nrow = TM, ncol = N)
  eps_p_mat     = matrix(rnorm(N * TM, 0, 0.1), nrow = TM, ncol = N)
  
  probabilities = pmax(pmin(plogis(FE_fcts_mat+ eps_p_mat), 0.99), 0.1)
  
  D_mat         = matrix(rbinom(N * TM, size = 1, prob = as.vector(probabilities)), nrow = TM, ncol = N)
  Y_mat         = mu + beta_true * D_mat + FE_fcts_mat + eps_y_mat
  
  ###########  Estimating  #################
  KSS_obj = KSS(Y_mat ~ -1 + D_mat, additive.effects = "time")
  beta_hat = KSS_obj$slope.para
  sigma_hat_std   = summary(KSS_obj)$coefficients[2]
  ##########################################

  ci_se = c(beta_hat - (z_value * sigma_hat_std), beta_hat + (z_value * sigma_hat_std))

  list(CI_se = ci_se, CI_sigma = ci_sigma, theta = beta_hat, sigma = sigma_hat_std)
})

stopCluster(cl)

confint_matrix_se      = do.call(cbind, lapply(dml, `[[`, "CI_se"))
theta               = sapply(dml, `[[`, "theta")
sigma               = sapply(dml, `[[`, "sigma")
coverage_proportion_se = mean((confint_matrix_se[1, ] <= beta_true & beta_true <= confint_matrix_se[2, ]))

et   = Sys.time()
time = (et - st) 

print(paste0("Time taken = ", time))
print(paste0("Coverage Proportion SE = ", coverage_proportion_se))

write.table(confint_matrix_se, paste0("Output/KSS/Confidence_intervals_SE for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(theta, paste0("Output/KSS/Theta_Hat for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sigma, paste0("Output/KSS/Sigma_Hat for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
