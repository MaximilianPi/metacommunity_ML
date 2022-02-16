#' Simulation function for different meta-community scenarios
#' 
#' @param N number of plots
#' @param SP number of (initial) species
#' @param E number of environmental factors
#' @param time number of time points
#' @param scenario which scenario
#' @param upperA upper boundary of A entries
#' @param diagA diagonal values for interactions within species
#' 
#' @details 
#' 
#' Simple MAR models to generate benchmark data for the ML approach
#' 

simulate  = function(N = 100, 
                     SP = 5, 
                     E = 3, 
                     time = 100, 
                     scenario = c("temporal", "spatial", "spatio-temporal"), 
                     upperA = 0.5,
                     diagA = 0.0,
                     count = FALSE, 
                     A = NULL
                     ) {
  
  scenario = match.arg(scenario)
  
  # Interaction Matrix
  if(is.null(A)) {
    A = matrix(runif(SP*SP, 0.0, upperA), SP, SP)
    W = matrix(rnorm(SP*E), E, SP)
    diag(A) = diagA
  }
  
  if(!count) {
    link = function(V) 1/(1+exp(-V))
    
    sample_f = function(P) apply(P, 1:2, function(v) rbinom(1, 1, link(v) ))
  } else {
    link = function(V) exp(V)
    
    sample_f = function(P) apply(P, 1:2, function(v) rpois(1, link(v) ))
  }
  
  if(scenario == "temporal") {
    # one plot? I think we agreed on this, right?
    N = 1
    X = matrix(runif(N*E, -1 ,1), N, E)
    W = matrix(runif(SP*E,-1, 1), E, SP)
  
    Y = sample_f(X%*% W)
    comms = vector("list", time)
    Xs = vector("list", time)
    comms[[1]] = Y
    Xs[[1]] = X
    for(i in 2:time) {
      comms[[i]] =  ( (comms[[i-1]]%*%A) + X%*% W )  + mvtnorm::rmvnorm(1, rep(0, SP), sigma = diag(0.3, SP))
      Xs[[i]] = X
    } 
    comms = lapply(comms, link)
  }
   if(scenario == "spatio-temporal") {
     X = matrix(runif(N*E, -1 ,1), N, E)
     W = matrix(runif(SP*E, -1, 1), E, SP)
     
     Y = sample_f(X%*% W)
     comms = vector("list", time)
     Xs = vector("list", time)
     comms[[1]] = Y
     Xs[[1]] = X
     for(i in 2:time) {
       comms[[i]] =  ( (comms[[i-1]]%*%A) + X%*% W  + mvtnorm::rmvnorm(N, rep(0, SP), sigma = diag(0.3, SP))) 
       #plot(sapply(comms, function(com) com[1,1])[20:100], type = "l")
       Xs[[i]] = X
     }
     comms = lapply(comms, link)
   }
    
  
    if(scenario == "spatio") {
      time = 1
      X = matrix(runif(N*E), N, E)
      W = matrix(rnorm(SP*E), E, SP)
      Y = sample_f(X%*% W)
      comms = vector("list", time)
      Xs = vector("list", time)
      comms[[1]] = Y
      Xs[[1]] = X
    }  

  return(list(A = A, W = W, X = Xs, Y = comms, time = time, N = N, scenario=scenario))
}





prepare_data = function(data, lag = 1, scale_d=TRUE) {
  time_steps = length(data$X)
  sites = nrow(data$X[[1]])
  XX = do.call(rbind, data$X)
  YY = do.call(rbind, data$Y)
  lag_indices = (1: (lag*sites) )#+(lag-1)*sites
  X_corrected = cbind(XX[-lag_indices,], YY[-rev(nrow(YY)-lag_indices+1),])
  colnames(X_corrected) = c(paste0("E_", 1:ncol(XX)), paste0("Y_", 1:ncol(YY)))
  YY = YY[-lag_indices,]
  XX = X_corrected
  if(scale_d) XX = scale(XX)
  XX = cbind(XX[,1:ncol(data$X[[1]])], apply(XX[,(ncol(data$X[[1]])+1):ncol(XX)], 2, function(x) (x+min(x))/(max(x+min(x)))))
  return(list(X = XX, Y = round(YY)))
}
