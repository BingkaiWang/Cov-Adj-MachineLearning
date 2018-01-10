set.seed(100)
library(dplyr)
library(purrr)
library(boot)

load("METS.rdata")
source(file = "Rotnitzky_estimator.R")
METS <- subset(METS, !is.na(METS$weightchange))
y <- METS$weightchange
a <- METS$treatment == "Metformin"
w <- select(METS, age, gender, CGI, tobacco, drug, alcohol, weight_baseline, bmi)
mets <- cbind(y,a,w)

# unadjusted estimator
boot_unadj <- boot(mets, statistic = function(data, i){
  data = data[i,]
  mean(data$y[data$a]) - mean(data$y[!data$a])
}, R = 1000, stype = 'i')
bias_unadj <- mean(boot_unadj$t) - boot_unadj$t0
var_unadj <- var(boot_unadj$t)

# ANCOVA with all variable (LOOCV)
pred_y <- map_dbl(1:nrow(mets), function(i){
  predict(lm(y~., data = mets[-i,]), mets[i,])
})
var_ancova_all <- 4*mean((mets$y - pred_y)^2)/nrow(mets)

# ANCOVA with only baseline score (LOOCV)
pred_y <- map_dbl(1:nrow(mets), function(i){
  predict(lm(y~a+weight_baseline, data = mets[-i,]), mets[i,])
})
var_ancova_baseline <- 4*mean((mets$y - pred_y)^2)/nrow(mets)

# Rotnitzky's estimator
boot_rotnitzky <- boot(mets, statistic = function(data, i){
  data = data[i,]
  beta_1 <- suppressWarnings(DR(Y = data$y, A = data$a, piX = as.matrix(data[,-(1:2)]), phiX = as.matrix(data[,-(1:2)])))$betaDR
  beta_0 <- suppressWarnings(DR(Y = data$y, A = !data$a, piX = as.matrix(data[,-(1:2)]), phiX = as.matrix(data[,-(1:2)])))$betaDR
  beta_1 - beta_0
}, R = 1000, stype = 'i')
bias_rotnitzky <- mean(boot_rotnitzky$t) - boot_rotnitzky$t0
var_rotnitzky <- var(boot_rotnitzky$t)

# result
METS_result <- data.frame(bias = c(bias_unadj,NA,NA, bias_rotnitzky),
                          var = c(var_unadj,var_ancova_baseline, var_ancova_all, var_rotnitzky),
                          RE = c(1,var_unadj/var_ancova_baseline, var_unadj/var_ancova_all, var_unadj/var_rotnitzky))
round(METS_result,2)