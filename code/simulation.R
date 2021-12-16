#' Simulation function for different meta-community scenarios
#' 
#' @param N number of plots
#' @param SP number of (initial) species
#' @param E number of environmental factors
#' @param time number of time points
#' @param scenario which scenario
#' 
#' @details 
#' 
#' Simple MAR models to generate benchmark data for the ML approach
#' 

simulate  = function(N = 100, SP = 5, E = 3, time = 100, scenario = c("temporal", "spatial", "spatio-temporal")) {
  
  scenario = match.arg(scenario)
  
  # Interaction Matrix
  A = matrix(runif(SP*SP, -0.5, 0.5), SP, SP)
  W = matrix(rnorm(SP*E), E, SP)
  diag(A) = 1.0
  
  if(scenario == "temporal") {
    # one plot? I think we agreed on this, right?
    N = 1
    X = matrix(runif(N*E), N, E)
    W = matrix(rnorm(SP*E), E, SP)
  
    Y = sample_f(X%*% W)
    comms = vector("list", time)
    Xs = vector("list", time)
    comms[[1]] = Y
    Xs[[1]] = X
    for(i in 2:t) {
      comms[[i]] =  sample_f( comms[[i-1]]%*%A + X%*% W  )
      Xs[[i]] = X
    } 
  }
   if(scenario == "spatio-temporal") {
     X = matrix(runif(N*E), N, E)
     W = matrix(rnorm(SP*E), E, SP)
     
     Y = sample_f(X%*% W)
     comms = vector("list", time)
     Xs = vector("list", time)
     comms[[1]] = Y
     Xs[[1]] = X
     for(i in 2:t) {
       comms[[i]] =  sample_f( comms[[i-1]]%*%A + X%*% W  )
       Xs[[i]] = X
     }
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



link = function(V) 1/(1+exp(-V))

sample_f = function(P) apply(P, 1:2, function(v) rbinom(1, 1, link(v) ))