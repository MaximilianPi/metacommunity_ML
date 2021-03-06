library(sjSDM)
library(xgboost)
library(ranger)
source("code/simulation.R")
source("code/random_forest.R")
source("code/brt.R")

set.seed(42)

# Comparison between MAR, BRT, and RF
extract_results = function(data) {
  MAR_p = sjSDM(data$count_prepared$Y, data$count_prepared$X, sampling = 100L, family = poisson(), device = sample(c(0,1), 1) )
  BRT_p = community_BRT_xg(data$count_prepared$X, data$count_prepared$Y, E = nrow(data$data_count$W), response = "poisson")
  RF_p = community_RF(data$count_prepared$X, data$count_prepared$Y, E = nrow(data$data_count$W), response = "poisson", impurity = "permutation")
  
  MAR_b = sjSDM(data$binary_prepared$Y, data$binary_prepared$X, sampling = 100L, family = binomial(), device = sample(c(0,1), 1) )
  BRT_b = community_BRT_xg(data$binary_prepared$X, data$binary_prepared$Y, E = nrow(data$data_binary$W), response = "binomial")
  RF_b = community_RF(data$binary_prepared$X, data$binary_prepared$Y, E = nrow(data$data_binary$W), response = "binomial", impurity = "permutation")
  
  MAR_A_p = t(coef(MAR_p)[[1]])[(nrow(data$data_count$W)+2):(nrow(t(coef(MAR_p)[[1]]))), ]
  MAR_A_b = t(coef(MAR_b)[[1]])[(nrow(data$data_binary$W)+2):(nrow(t(coef(MAR_b)[[1]]))), ]
  BRT_A_p = t(BRT_p$A)
  BRT_A_b = t(BRT_b$A)
  RF_A_p = t(RF_p$A)
  RF_A_b = t(RF_b$A)
  
  rmse_f = function(true, obs) sqrt(mean((true-obs)**2))
  rmse_p = 
    sapply(list(MAR_A_p, 
                BRT_A_p, 
                RF_A_p, 
                matrix(runif(nrow(MAR_A_p)**2, 0, 0.8 ),nrow(MAR_A_p), nrow(MAR_A_p) ),
                diag(0.0, nrow(MAR_A_b))
    ), function(Ap) rmse_f(abs(data$data_count$A), abs(Ap)))
  correlation_p = 
    sapply(list(MAR_A_p, 
                BRT_A_p, 
                RF_A_p, 
                matrix(runif(nrow(MAR_A_p)**2, 0, 0.8 ),nrow(MAR_A_p), nrow(MAR_A_p) ),
                diag(0.0, nrow(MAR_A_b))
    ), function(Ap) cor(as.vector(abs(data$data_binary$A)), as.vector(abs(Ap)), method = "spearman"))
  
  rmse_b = 
    sapply(list(MAR_A_b, 
                BRT_A_b, 
                RF_A_b, 
                matrix(runif(nrow(MAR_A_p)**2, 0, 0.8 ),nrow(MAR_A_p), nrow(MAR_A_p) ),
                diag(0.0, nrow(MAR_A_b))
    ), function(Ap) rmse_f(abs(data$data_binary$A), abs(Ap)))
  
  correlation_b = 
    sapply(list(MAR_A_b, 
                BRT_A_b, 
                RF_A_b, 
                matrix(runif(nrow(MAR_A_p)**2, 0, 0.8 ),nrow(MAR_A_p), nrow(MAR_A_p) ),
                diag(0.0, nrow(MAR_A_b))
    ), function(Ap) cor(abs(as.vector(data$data_binary$A)), as.vector(abs(Ap)), method = "spearman"))
  gc()
  sjSDM:::pkg.env$torch$cuda$empty_cache()
  return(cbind(rmse_p, rmse_b, correlation_p, correlation_b))
}


data_sets = 
  lapply(c(5, 10, 20), function(species) {
    N = species
    
    A1 = diag(0.5,N)
    A2 = matrix(runif(N**2, -0.3, 0.3), N, N)
    diag(A2) = 0.0
    A3 = matrix(runif(N**2, -0.2, 0.2), N, N)
    diag(A3) = runif(N, 0.25, 0.5)
    A4 = diag(0.0, N)
    for(i in 1:(N-1)) {
      A4[i+1, i] = runif(1, 0, 0.7)
      A4[i, i+1] = runif(1, 0, 0.7)
    }
    A5 = diag(0.0, N)
    
    data_set = 
      lapply(1:20, function(i) {
        data = 
          lapply(list(A1, A2, A3, A4, A5), function(A) {
            data_count = data = simulate(scenario ="spatio-temporal", SP = N, count = TRUE, N = 50, time = 100, E = 2,A = A)
            data_binary = data = simulate(scenario ="spatio-temporal", SP = N, count = FALSE, N = 50, time = 100, E = 2, A = A)
            count_prepared = prepare_data(data_count)
            binary_prepared = prepare_data(data_binary)
            return(list(data_count = data_count, data_binary = data_binary, count_prepared = count_prepared, binary_prepared = binary_prepared))
          })
        return(data)
      })
  
  })

lapply(data_sets, function(d) lapply(d, function(dd) sapply(dd, function(ddd) sum(ddd$data_count$Y[[100]]) ) ))

cl = parallel::makeCluster(10L)
parallel::clusterExport(cl, varlist = list("extract_results", "community_RF", "community_BRT_xg"), envir = environment())
parallel::clusterEvalQ(cl, {library(sjSDM);library(xgboost);library(ranger)})

for(sp in 1:4) {
  sp_data = data_sets[[sp]]
  for(a in 1:5) {
    parallel::clusterExport(cl, list("sp_data", "a", "sp"), envir = environment())
    res = abind::abind(parallel::parLapply(cl, 1:20, function(i)  extract_results(sp_data[[i]][[a]] ) ), along = 0)
    saveRDS(res, paste0("results/","a_", a,"_sp_", sp,"_.RDS"))
  }
}

