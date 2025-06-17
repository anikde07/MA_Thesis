args       = commandArgs(trailingOnly = TRUE)
N          = as.numeric(args[1])
TM         = as.numeric(args[2])
MC         = 10000

print(paste0("PLM only. N = ", N, " MC = ", MC, " T = ", TM))

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
clusterEvalQ(cl, c(library("phtt"), library("ffscb"), library("pspline"), library("plm")))

dml <- parLapply(cl, 1:MC, function(v) {
  
  gc()

  FE_fcts_mat   = make_sample(mean.v = rep(0, TM), cov.m = cov.m, N = N, dist = "rnorm")
  eps_y_mat     = matrix(rnorm(N * TM, 0, 0.2), nrow = TM, ncol = N)
  eps_p_mat     = matrix(rnorm(N * TM, 0, 0.1), nrow = TM, ncol = N)
  
  probabilities = pmax(pmin(plogis(FE_fcts_mat+ eps_p_mat), 0.99), 0.1)
  
  D_mat         = matrix(rbinom(N * TM, size = 1, prob = as.vector(probabilities)), nrow = TM, ncol = N)
  Y_mat         = mu + beta_true * D_mat + FE_fcts_mat + eps_y_mat

  df            = data.frame(time = rep(1:TM, each = N), individual = rep(1:N, times = TM),
                  Y = as.vector(Y_mat), D = as.vector(D_mat))

  ###########  Estimating  #################
  model_plm <- plm(Y ~ D, data = df, index = c("individual", "time"), model = "within", effect = "twoways")
  beta_hat = summary(model_plm)$coefficients[1]
  se_hat   = summary(model_plm)$coefficients[2]
  ##########################################

  ci = c(beta_hat - (z_value * se_hat), beta_hat + (z_value * se_hat))
  
  list(CI = ci, theta = beta_hat, sigma = se_hat)
})

stopCluster(cl)

confint_matrix      = do.call(cbind, lapply(dml, `[[`, "CI"))
theta               = sapply(dml, `[[`, "theta")
sigma               = sapply(dml, `[[`, "sigma")
coverage_proportion = mean((confint_matrix[1, ] <= beta_true & beta_true <= confint_matrix[2, ]))

et   = Sys.time()
time = (et - st) 

print(paste0("Time taken = ", time))
print(paste0("Coverage Proportion = ", coverage_proportion))

write.table(confint_matrix, paste0("Output/Twoways/Confidence_intervals for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(theta, paste0("Output/Twoways/Theta_Hat for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sigma, paste0("Output/Twoways/Sigma_Hat for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
