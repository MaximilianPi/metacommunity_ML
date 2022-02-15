library(sjSDM)
source("code/simulation.R")
source("code/random_forest.R")
source("code/brt.R")
data = simulate(scenario ="spatio-temporal", diagA = 0.1,SP = 5, upperA = 0.3, count = TRUE, N = 50, time = 100, E = 3)

prepared_data = prepare_data(data, lag = 1)


X = scale(prepared_data$X)
X[,1:3] = prepared_data$X[,1:3]
X[,4:8] = apply(prepared_data$X[,4:8], 2, function(x) (x+min(x))/(max(x+min(x))))
Y = round(prepared_data$Y)

### sjSDM ###
m = sjSDM(Y, X, sampling = 100L, family = poisson())
end = nrow(data$W)+1+ncol(data$A)
cor(as.vector(t(coef(m)[[1]][,5:end])), as.vector(data$A), method = "spearman")


### RF ###
results_RF2 = community_RF(X, Y, 3, impurity = "permutation", response = "count")

cor(abs(as.vector(results_RF$A)), as.vector(data$A), method = "spearman")
cor(abs(as.vector(results_RF2$W)), abs(as.vector(t(data$W))), method = "spearman")


### BRT ###
results_BRT = community_BRT(X, Y, 3, response = "count")
cor(abs(as.vector(results_BRT$A)), as.vector(t(data$A)), method = "spearman")
cor(abs(as.vector(results_BRT$W)), abs(as.vector(t(data$W))), method = "spearman")


### BRT 2 ###
results_BRT = community_BRT_xg(X, Y, 3, response = "count", correct = FALSE)
cor(abs(as.vector(results_BRT$A)), as.vector(t(data$A)), method = "spearman")
cor(abs(as.vector(results_BRT$W)), abs(as.vector(t(data$W))), method = "spearman")
