library(sjSDM)
source("code/simulation.R")
source("code/random_forest.R")
source("code/brt.R")

AT = diag(0, 20)
#AT[lower.tri(AT)] = runif(sum(lower.tri(AT)), 0, 0.3)
for(i in 1:19) {
  AT[i+1, i] = runif(1, 0, 0.7)
  AT[i, i+1] = runif(1, 0, 0.7)
}
AT = diag(0, 20)
data = simulate(scenario ="spatio-temporal", 
                diagA = 0.1,
                SP = 20, 
                upperA = 0.3, 
                count = TRUE, 
                N = 50, 
                time = 100, 
                E = 3,
                A = AT)

prepared_data = prepare_data(data, lag = 1)


X = scale(prepared_data$X)
X[,1:3] = prepared_data$X[,1:3]
X[,4:8] = apply(prepared_data$X[,4:8], 2, function(x) (x+min(x))/(max(x+min(x))))
Y = round(prepared_data$Y)

### sjSDM ###
m = sjSDM(Y, X, sampling = 100L, family = binomial())
end = nrow(data$W)+1+ncol(data$A)
cor(as.vector(t(coef(m)[[1]][,5:end])), as.vector(data$A), method = "spearman")


### RF ###
results_RF2 = community_RF(X, Y, 3, impurity = "permutation", response = "count")
cor(abs(as.vector(results_RF2$A)), as.vector(t(data$A)), method = "spearman")
cor(abs(as.vector(results_RF2$W)), abs(as.vector(t(data$W))), method = "spearman")


### BRT ###
results_BRT = community_BRT(X, Y, 3, response = "count")
cor(abs(as.vector(results_BRT$A)), as.vector(t(data$A)), method = "spearman")
cor(abs(as.vector(results_BRT$W)), abs(as.vector(t(data$W))), method = "spearman")


### BRT 2 ###
results_BRT = community_BRT_xg(X, Y, 3, response = "poisson", correct = FALSE)
AA = results_BRT$A
AA[AA<0] = 0
cor((as.vector(AA)), as.vector(t(data$A)), method = "spearman")
cor(abs(as.vector(results_BRT$W)), abs(as.vector(t(data$W))), method = "spearman")


