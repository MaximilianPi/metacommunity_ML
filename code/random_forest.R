library(ranger)
library(xgboost)

#' Community Random Forest 
#' 
#' @param XX predictors
#' @param YY N*SP species occurrence matrix
#' @param E number of environmental predictors (without intercept)
#' @param impurity which impurity
#' @param response poisson or binomial
#' 
#' @import ranger
community_RF = function(XX, YY, E = 3, impurity = "impurity", response = "binomial") {
  SP = ncol(YY)
  models = null = preds = vector("list", SP)

  x_data = cbind(rnorm(nrow(XX)), XX)
  for(i in 1:SP) {
    models[[i]] = ranger(x = x_data, y = YY[,i,drop=FALSE], 
                         importance = impurity, 
                         probability = ifelse(response == "binomial", TRUE, FALSE), 
                         classification = ifelse(response == "binomial", TRUE, FALSE),
                         num.threads = 2L)
    
    if(response == "binomial") {
      preds[[i]] = predict(models[[i]], data = x_data)$predictions[,1] # return probabilities for class 1
    } else {
      
    preds[[i]] = predict(models[[i]], data = x_data)$predictions
    }

  }
  
  ## naive ###
  # biotic importance matrix
  A = matrix(NA, SP, SP)
  for(i in 1:SP) {
    imp = (models[[i]]$variable.importance-models[[i]]$variable.importance[1])[-1]
    A[i, ] = imp[(E+1):length(imp)] 
  }
  
  # spatial importance matrix
  W = matrix(NA, SP, E)
  for(i in 1:SP) {
    imp = (models[[i]]$variable.importance-models[[i]]$variable.importance[1])[-1]
    W[i, ] = imp[1:(E)] 
  }
  
  # species-species association matrix
  preds = do.call(cbind, preds)
  Sigma = cov(YY - preds)
  
  
  result = list(models = models, 
                A = A, 
                W = W, 
                null = null, 
                Sigma = Sigma, 
                Pred = preds)
  class(result) = c("community_rf")
  
  return(result)
}
coef.community_rf = function(object, ...) object$W

predict.community_rf = function(object, data, ...) {
  preds = list()
  for(i in 1:length(object$models)) {
    preds[[i]] = predict(object$models[[i]], data = data)$predictions
  }
  return(do.call(cbind, preds))
}
