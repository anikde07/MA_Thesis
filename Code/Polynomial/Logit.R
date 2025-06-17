args       = commandArgs(trailingOnly = TRUE)
N          = as.numeric(args[1])
TM         = as.numeric(args[2])
K          = 5
MC         = 10000

print(paste0("Logit with KSS being used to estimate X. N = ", N, " K = ", K, " MC = ", MC, " T = ", TM))

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

bernstein_basis <- function(n, x) 
{
  sapply(0:n, function(i) choose(n, i) * x^i * (1 - x)^(n - i))
}

## Parallelization cluster ##

num_cores = detectCores()
cl        = makeCluster(num_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, c(library("phtt"), library("ffscb"), library("pspline"), library("caret"), library("foreach"), library("plm")))

######################

dml <- parLapply(cl, 1:MC, function(v) {
  
  gc()
  
  FE_fcts_mat = matrix(NA, nrow = TM, ncol = N)
  a_values <- rnorm(N, mean = 0, sd = 1)
  b_values <- rnorm(N, mean = 0, sd = 1)
  c_values <- rnorm(N, mean = 0, sd = 1)

  for (i in 1:N) 
  {
    coefficients    <- c(a_values[i], b_values[i], c_values[i])
    degree          <- length(coefficients) - 1  
    basis           <- bernstein_basis(degree, grid)
    FE_fcts_mat[,i] <- basis %*% coefficients  
  }

  eps_y_mat     = matrix(rnorm(N * TM, 0, 0.2), nrow = TM, ncol = N)
  eps_p_mat     = matrix(rnorm(N * TM, 0, 0.1), nrow = TM, ncol = N)
  
  probabilities = pmax(pmin(plogis(FE_fcts_mat + eps_p_mat), 0.99), 0.1)
  
  D_mat         = matrix(rbinom(N * TM, size = 1, prob = as.vector(probabilities)), nrow = TM, ncol = N)
  Y_mat         = mu + beta_true * D_mat + FE_fcts_mat + eps_y_mat
  
  ####### Estimating \hat{X}_it ############
  
  KSS_obj = KSS(Y_mat ~ -1 + D_mat, additive.effects = "time")
  X_hat_mat   = KSS_obj$unob.fact.stru

  ####### Cross-Fitting ##############
  stratification_vector = sample(rep(1:K, length.out = N))
  skf                   = createFolds(stratification_vector, k = K, list = TRUE)

  cross_fit <- foreach(i = 1:K, .combine = 'list', .multicombine = TRUE, .maxcombine = K) %do% {
    
    train_indices = unlist(skf[-i])
    eval_indices  = unlist(skf[i])
    
    y_train = Y_mat[, train_indices]
    d_train = D_mat[, train_indices]
    x_train = X_hat_mat[, train_indices]
    
    y_eval  = Y_mat[, eval_indices]
    d_eval  = D_mat[, eval_indices]
    x_eval  = X_hat_mat[, eval_indices]
    
    y_train_sub           = y_train - x_train
    y_train_sub_minus_bar = y_train_sub - rowMeans(y_train_sub)
    d_train_minus_bar     = d_train - rowMeans(d_train)
    
    y_train_sub_minus_bar = matrix(y_train_sub_minus_bar, dimnames = list(t(outer(colnames(y_train_sub_minus_bar), rownames(y_train_sub_minus_bar), FUN = paste)), NULL))
    d_train_minus_bar     = matrix(d_train_minus_bar, dimnames = list(t(outer(colnames(d_train_minus_bar), rownames(d_train_minus_bar), FUN = paste)), NULL))
    
    ## Training g_0 ##
    
    model    = lm(y_train_sub_minus_bar ~ d_train_minus_bar)
    beta_hat = as.numeric(coef(model)[2])
    mu_hat   = rowMeans(y_train_sub) - beta_hat * rowMeans(d_train)
    
    ### Evaluating \hat{g_0}_(0,.) and \hat{g_0}_(1,.) at x_eval ### 
    
    g_0_hat_d_0 = mu_hat + x_eval
    g_0_hat_d_1 = mu_hat + beta_hat + x_eval
    
    ### Training \hat{m_0} ###
    
    m_0_hat = NULL
    
    for (j in 1:TM)
    {
      model_j  = glm(d ~ x, data = data.frame(d = d_train[j,], x = x_train[j,]), family = binomial)
      
      m_0t_hat = predict(model_j, newdata = data.frame(x = x_eval[j,]), type = "response")
      m_0t_hat = pmin(pmax(m_0t_hat, 0.01), 0.99)
      
      m_0_hat  = rbind(m_0_hat, m_0t_hat)
    }
    
    # Compute scores
    score_i_t_k = (g_0_hat_d_1 - g_0_hat_d_0) + ((d_eval * (y_eval - g_0_hat_d_1)) / m_0_hat) - (((1 - d_eval) * (y_eval - g_0_hat_d_0)) / (1 - m_0_hat))

    list(score_i_t_k = score_i_t_k)
  }

  scores_list_of_matrices = lapply(cross_fit, `[[`, "score_i_t_k")
  
  theta_t_k_hat_matrix = do.call(cbind, lapply(scores_list_of_matrices, rowMeans))
  theta_t_hat = rowMeans(theta_t_k_hat_matrix)
  theta_hat = mean(theta_t_hat)
  
  sigma_squared_list_of_matricies <- lapply(scores_list_of_matrices, function(score_matrix) {
    (score_matrix - theta_hat)^2
  })
  sigma_squared_t_k_hat_matrix = do.call(cbind, lapply(sigma_squared_list_of_matricies, rowMeans))
  sigma_squared_t_hat = rowMeans(sigma_squared_t_k_hat_matrix)
  sigma_squared_hat = mean(sigma_squared_t_hat)
  
  ci = c(theta_hat - (z_value * sqrt(sigma_squared_hat/(TM*N))), theta_hat + (z_value * sqrt(sigma_squared_hat/(TM*N))))
  
  list(CI = ci, theta = theta_hat, sigma_squared = sigma_squared_hat)
})

stopCluster(cl)

confint_matrix      = do.call(cbind, lapply(dml, `[[`, "CI"))
theta               = sapply(dml, `[[`, "theta")
sigma_squared       = sapply(dml, `[[`, "sigma_squared")
coverage_proportion = mean((confint_matrix[1, ] <= beta_true & beta_true <= confint_matrix[2, ]))

et   = Sys.time()
time = (et - st)

print(paste0("Time taken = ", time))
print(paste0("Coverage Proportion = ", coverage_proportion))

write.table(confint_matrix, paste0("Output/Logit/Confidence_intervals for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(theta, paste0("Output/Logit/Theta_Hat for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
write.table(sigma_squared, paste0("Output/Logit/Sigma_Squared_Hat for N = ", N, "T = ", TM, ".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
