library(ranger)

#' Community Random Forest 
#' 
#' @param XX predictors
#' @param YY N*SP species occurrence matrix
#' @param E number of environmental predictors (without intercept)
#' 
#' @import ranger
community_RF = function(XX, YY, E = 3) {
  SP = ncol(YY)
  models = vector("list", SP)
  for(i in 1:SP) {
    models[[i]] = ranger(x = cbind(1, X_corrected), y = YY[,i,drop=FALSE], 
                         importance = "impurity_corrected", probability = TRUE)
  }
  
  ## naive ###
  A = matrix(NA, SP, SP)
  for(i in 1:SP) {
    A[i, ] = models[[i]]$variable.importance[(E+2):length(models[[i]]$variable.importance)]
  }
  
  W = matrix(NA, SP, E+1)
  for(i in 1:SP) {
    W[i, ] = models[[i]]$variable.importance[1:(E+1)]
  }
  
  result = list(models = models, A_impurity_corrected = A, W_impurity_corrected = W)
  class(result) = c("community_rf")
  
  return(result)
}
coef.community_rf = function(object, ...) object$W_impurity_corrected

predict.community_rf = function(object, data, ...) {
  preds = list()
  for(i in 1:length(object$models)) {
    preds[[i]] = predict(object$models[[i]], data = data)$predictions
  }
  return(do.call(cbind, preds))
}
