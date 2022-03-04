library(xgboost)


#' Community BRT
#'
#' @param XX predictors
#' @param YY N*SP species occurrence matrix
#' @param response response type, binomial or poisson
#' @param MAR fit multivariate autoregressive or not
#' @param nrounds number of trees
#'
#' @import ranger
community_BRT_xg = function(XX, YY,response = "binomial", MAR=TRUE, nrounds = 50, ...) {
  XX = as.matrix(XX)
  YY = as.matrix(YY)
  library(xgboost)
  
  SP = ncol(YY)
  models = preds = vector("list", SP)
  dist = ifelse(response == "binomial", "binary:logistic", "count:poisson" )
  
  for(i in 1:SP) {
    #df = data.frame(XX, Y_target = YY[,i,drop=FALSE])
    df = xgboost::xgb.DMatrix(XX, label = YY[,i,drop=FALSE])
    models[[i]]  = xgboost::xgboost(df, nrounds = nrounds, objective = dist, max_depth = 5, early_stopping_rounds = 1, nthread = 2)
    preds[[i]] = predict(models[[i]], df) # return probabilities for class 1
    
  }
  E=0
  if(MAR) {
    E = ncol(XX) - ncol(YY)
  }
  
  
  A = matrix(NA, SP, SP)
  W = matrix(NA, SP, E)
  if(MAR) {
    ## naive ###
    # biotic importance matrix
    for(i in 1:SP) {
      imps = tryCatch({xgboost::xgb.importance(model = models[[i]])}, error = function(e) e)
      if(!inherits(imps, "error")) {
        importances = imps$Gain
        names(importances) = imps$Feature
        imps = importances[colnames(XX)]
        if(correct) {
          imps = (imps-imps[1])[-1]
          nns = colnames(df)[-1]
        } else {
          nns = colnames(df)
        }
        all_in = (nns %in% names(imps))
        if(any(!all_in)){
          fill = rep(0, length(nns[!all_in]))
          names(fill) = nns[!all_in]
          imps = c(imps, fill)
        }
        imps = imps[colnames(XX)]
        A[i, ] = imps[(E+1):length(imps)] - min(imps)
        W[i, ] = imps[1:E] - min(imps)
      }
    }
  } else {
    for(i in 1:SP) {
      imps = tryCatch({xgboost::xgb.importance(model = models[[i]])}, error = function(e) e)
      if(!inherits(imps, "error")) {
        importances = imps$Gain
        names(importances) = imps$Feature
        imps = importances[colnames(XX)]
        if(correct) {
          imps = (imps-imps[1])[-1]
          nns = colnames(df)[-1]
        } else {
          nns = colnames(df)
        }
        all_in = (nns %in% names(imps))
        if(any(!all_in)){
          fill = rep(0, length(nns[!all_in]))
          names(fill) = nns[!all_in]
          imps = c(imps, fill)
        }
        imps = imps[colnames(XX)]
        W[i, ] = imps - min(imps)
      }
    }
  }
  
  
  # species-species association matrix
  preds = do.call(cbind, preds)
  Sigma = cov(YY - preds)
  
  
  result = list(models = models,
                A = A,
                W = W,
                Sigma = Sigma,
                Pred = preds)
  class(result) = c("community_brt_xg")
  
  return(result)
}
