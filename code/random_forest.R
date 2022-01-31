library(ranger)

#' Community Random Forest 
#' 
#' @param XX predictors
#' @param YY N*SP species occurrence matrix
#' @param E number of environmental predictors (without intercept)
#' @param impurity which impurity
#' 
#' @import ranger
community_RF = function(XX, YY, E = 3, impurity = "impurity") {
  SP = ncol(YY)
  models = null = preds = vector("list", SP)
  XX_n = apply(XX, 2, sample)
  YY_n = apply(YY, 2, sample)
  for(i in 1:SP) {
    models[[i]] = ranger(x = cbind(1, XX), y = YY[,i,drop=FALSE], 
                         importance = impurity, probability = TRUE)
    
    preds[[i]] = predict(models[[i]], data = cbind(1, XX))$predictions[,1] # return probabilities for class 1
    
    # maybe useful?
    null[[i]]   = ranger(x = cbind(1, XX_n), y = YY_n[,i,drop=FALSE], 
                         importance = impurity, probability = TRUE)
  }
  
  ## naive ###
  # biotic importance matrix
  A = matrix(NA, SP, SP)
  for(i in 1:SP) {
    A[i, ] = models[[i]]$variable.importance[(E+2):length(models[[i]]$variable.importance)] 
  }
  
  # spatial importance matrix
  W = matrix(NA, SP, E+1)
  for(i in 1:SP) {
    W[i, ] = models[[i]]$variable.importance[1:(E+1)] 
  }
  
  # species-species association matrix
  preds = do.call(cbind, preds)
  Sigma = cov(YY - preds)
  
  
  result = list(models = models, 
                A_impurity_corrected = A, 
                W_impurity_corrected = W, 
                null = null, 
                Sigma = Sigma, 
                Pred = preds)
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

