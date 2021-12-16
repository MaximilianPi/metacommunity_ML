source("code/simulation.R")
source("code/random_forest.R")
data = simulate(scenario ="spatio-temporal")

### Transform data ###
# TODO create function for the following chunk!

XX = do.call(rbind, data$X)
YY = do.call(rbind, data$Y)
X_corrected = cbind(XX[-c(1:100),], YY[-c( (nrow(YY)-100 + 1):nrow(YY) ),])
colnames(X_corrected) = c(paste0("E_", 1:ncol(XX)), paste0("Y_", 1:ncol(YY)))
YY = YY[-c(1:100),]
XX = X_corrected

### sjSDM ###
library(sjSDM)
m = sjSDM(YY, X_corrected, sampling = 100L)
fields::image.plot(t(coef(m)[[1]][,5:9]))
fields::image.plot(data$A)
cor(as.vector(t(coef(m)[[1]][,5:9])), as.vector(data$A), method = "spearman")
plot(as.vector(t(coef(m)[[1]][,5:9])), as.vector(data$A))


### RF ###
results_RF = community_RF(XX, YY, 3)
# the signs don't tell us anything, I think, tbh I have first to read up this 'impurity_corrected' 
# but I guess that we should use a general permutation approach anyways to make it independent of the ML algorithm!

cor(as.vector(results_RF$A_impurity_corrected), as.vector(data$A), method = "spearman")
plot(as.vector(results_RF$A_impurity_corrected), as.vector(data$A))
