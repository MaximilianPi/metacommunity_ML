library(sjSDM)
source("code/simulation.R")
source("code/random_forest.R")
data = simulate(scenario ="spatio-temporal", diagA = 0,SP = 5, upperA = 0.5)

### Transform data ###
# TODO create function for the following chunk!

XX = do.call(rbind, data$X)
YY = do.call(rbind, data$Y)
X_corrected = cbind(XX[-c(1:100),], YY[-c( (nrow(YY)-100 + 1):nrow(YY) ),])
colnames(X_corrected) = c(paste0("E_", 1:ncol(XX)), paste0("Y_", 1:ncol(YY)))
YY = YY[-c(1:100),]
XX = X_corrected

### sjSDM ###
m = sjSDM(YY, X_corrected, sampling = 100L)
end = nrow(data$W)+1+ncol(data$A)
fields::image.plot(t(coef(m)[[1]][,5:end]))
fields::image.plot(data$A)
cor(as.vector(t(coef(m)[[1]][,5:end])), as.vector(data$A), method = "spearman")
plot(as.vector(t(coef(m)[[1]][,5:end])), as.vector(data$A))


### RF ###
results_RF = community_RF(XX, YY, 3, impurity = "impurity_corrected")
# the signs don't tell us anything, I think, tbh I have first to read up this 'impurity_corrected' 
# but I guess that we should use a general permutation approach anyways to make it independent of the ML algorithm!

cor(abs(as.vector(results_RF$A_impurity_corrected)), as.vector(data$A), method = "spearman")
cor(abs(as.vector(results_RF$W_impurity_corrected[,-1])), abs(as.vector(data$W)), method = "spearman")




### Post-hoc VI ###
library(DALEX)
library(iml)
VIs = lapply(1:5, function(i) {
  pred <- function(model, newdata)  {
    return(predict(model$models[[i]], data=cbind(1, newdata))$predictions[,2])
  }
  
  explainer_rf = explain.default(results_RF, data = XX, 
                         predict_function = pred, y = YY[,i])
  
  
  VI = DALEX::variable_importance(explainer_rf, type = "raw")
  res = as.data.frame(VI) %>% group_by(variable) %>% summarise(n =  mean(dropout_loss))
  return(res)
})
AB = sapply(1:5, function(i) (VIs[[i]]$n - VIs[[i]]$n[1])[6:10])

cor(abs(as.vector(t(AB))), as.vector(data$A), method = "spearman")


VIs = lapply(1:5, function(i) {
  pred <- function(model, newdata)  {
    return(predict(model$models[[i]], data=cbind(1, as.matrix(newdata)))$predictions[,2])
  }
  
  P = Predictor$new(results_RF, data = data.frame(XX), predict.function = pred, y = YY[,i], type="prob")
  P$task = "classif"
  VI = iml::FeatureImp$new(P, loss = "rmse")
  KK = iml::FeatureEffects$new(P)
  
  VI = DALEX::variable_importance(explainer_rf)
  res = as.data.frame(VI) %>% group_by(variable) %>% summarise(n =  mean(dropout_loss))
  return(res)
})
AB = sapply(1:5, function(i) (VIs[[i]]$n - VIs[[i]]$n[1])[6:10])

cor(abs(as.vector((AB))), as.vector(data$A), method = "spearman")

AB = 
sapply(1:5, function(i) {
  PP = pred(results_RF, XX)
  return(coef(summary(lm(PP~XX)))[5:9,1])
})

