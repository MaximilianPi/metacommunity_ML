library(gbm)
library(xgboost)

#' Community BRT
#' 
#' @param XX predictors
#' @param YY N*SP species occurrence matrix
#' @param E number of environmental predictors (without intercept)
#' @param impurity which impurity
#' @param response binomial or poisson
#' 
#' @import ranger
community_BRT = function(XX, YY, E = 3, response = "binomial", ...) {
  SP = ncol(YY)
  models =  preds = vector("list", SP)
  dist = ifelse(response == "binomial", "bernoulli", "poisson")
  for(i in 1:SP) {
    df = data.frame(XX, Y_target = YY[,i,drop=FALSE])
    models[[i]]  = gbm::gbm(Y_target~., data = df, distribution = dist, 
                            interaction.depth = 3, cv = 2)
    preds[[i]] = predict(models[[i]]) # return probabilities for class 1
    
  }
  
  ## naive ###
  # biotic importance matrix
  A = matrix(NA, SP, SP)
  for(i in 1:SP) {
    imps = gbm::relative.influence(models[[i]]) #, n.trees = gbm.perf(models[[i]], method = "OOB", plot.it = FALSE))
    A[i, ] = imps[(E+1):length(imps)] 
  }
  
  # spatial importance matrix
  W = matrix(NA, SP, E)
  for(i in 1:SP) {
    imps = gbm::relative.influence(models[[i]]) #, n.trees = gbm.perf(models[[i]], method = "OOB", plot.it = FALSE))
    W[i, ] = imps[1:E] 
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
  class(result) = c("community_brt")
  
  return(result)
}



#' Community BRT
#'
#' @param XX predictors
#' @param YY N*SP species occurrence matrix
#' @param E number of environmental predictors (without intercept)
#' @param impurity which impurity
#'
#' @import ranger
community_BRT_xg = function(XX, YY, E = 3, response = "binomial",...) {
  SP = ncol(YY)
  models = null = preds = vector("list", SP)
  dist = ifelse(response == "binomial", "binary:logistic", "count:poisson" )
  for(i in 1:SP) {
    df = data.frame(XX, Y_target = YY[,i,drop=FALSE])
    df = xgboost::xgb.DMatrix(XX, label = YY[,i,drop=FALSE])
    models[[i]]  = xgboost::xgboost(df, nrounds = 100, objective = dist, max_depth = 3, early_stopping_rounds = 1)
    preds[[i]] = predict(models[[i]], df) # return probabilities for class 1

  }

  ## naive ###
  # biotic importance matrix
  A = matrix(NA, SP, SP)
  for(i in 1:SP) {
    imps = xgboost::xgb.importance(model = models[[i]])
    importances = imps$Gain
    names(importances) = imps$Feature
    imps = importances[order(names(importances))]
    A[i, ] = imps[(E+1):length(importances)]
  }

  # spatial importance matrix
  W = matrix(NA, SP, E)
  for(i in 1:SP) {
    imps = xgboost::xgb.importance(model = models[[i]])
    importances = imps$Gain
    names(importances) = imps$Feature
    imps = importances[order(names(importances))]
    W[i, ] = imps[1:E]
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
  class(result) = c("community_brt_xg")

  return(result)
}
