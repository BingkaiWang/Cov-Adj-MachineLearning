set.seed(100)
library(dplyr)
library(purrr)
library(boot)
library(SuperLearner)
library(foreach)
library(doParallel)

load("ACDS.rdata")
source(file = "Rotnitzky_estimator.R")
d <- subset(d, (d$arm != "Vitamin E") & (!is.na(d$Y18)))
y <- d$Y18 - d$Y0
a <- d$arm == "Donepezil"
w <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
mci <- cbind(y,a,w)
n1 <- sum(a)
n0 <- sum(!a)

# unadjusted estimator
boot_unadj <- map_dbl(1:10000, function(i){
  mci_boot1 <- sample_n(mci, size = n1, replace = T)
  mci_boot0 <- sample_n(mci, size = n0, replace = T)
  mean(mci_boot1$y) - mean(mci_boot0$y)
})
bias_unadj <- mean(boot_unadj)
var_unadj <- var(boot_unadj)
# boot_unadj1 <- boot(mci[mci$a,], statistic = function(data, i) mean(data[i,]$y), R = 10000, stype = 'i')
# boot_unadj0 <- boot(mci[!mci$a,], statistic = function(data, i) mean(data[i,]$y), R = 10000, stype = 'i')
# bias_unadj <- mean(boot_unadj1$t - boot_unadj0$t) - (boot_unadj1$t0 - boot_unadj0$t0)
# var_unadj <- var(boot_unadj1$t - boot_unadj0$t)

# ANCOVA with all variable (LOOCV)
boot_ancova_all <- map_dbl(1:10000, function(i){
  mci_boot1 <- sample_n(mci, size = n1, replace = T)
  mci_boot1$a <- T
  mci_boot0 <- sample_n(mci, size = n0, replace = T)
  mci_boot0$a <- F
  mci_boot <- rbind(mci_boot1, mci_boot0)
  lm(y~.,data = mci_boot)$coef[2]
})
bias_ancova_all <- mean(boot_ancova_all)
var_ancova_all <- var(boot_ancova_all)
# pred_y <- map_dbl(1:nrow(mci), function(i){
#   predict(lm(y~., data = mci[-i,]), mci[i,])
# })
# var_ancova_all <- 4*mean((mci$y - pred_y)^2)/nrow(mci)

# ANCOVA with only baseline score (LOOCV)
boot_ancova_baseline <- map_dbl(1:10000, function(i){
  mci_boot1 <- sample_n(mci, size = n1, replace = T)
  mci_boot1$a <- T
  mci_boot0 <- sample_n(mci, size = n0, replace = T)
  mci_boot0$a <- F
  mci_boot <- rbind(mci_boot1, mci_boot0)
  lm(y~a+Y0,data = mci_boot)$coef[2]
})
bias_ancova_baseline <- mean(boot_ancova_baseline)
var_ancova_baseline <- var(boot_ancova_baseline)
# pred_y <- map_dbl(1:nrow(mci), function(i){
#   predict(lm(y~a+Y0, data = mci[-i,]), mci[i,])
# })
# var_ancova_baseline <- 4*mean((mci$y - pred_y)^2)/nrow(mci)

# Rotnitzky's estimator
boot_rotnitzky <- map_dbl(1:10000, function(i){
  mci_boot1 <- sample_n(mci, size = n1, replace = T)
  mci_boot1$a <- T
  mci_boot0 <- sample_n(mci, size = n0, replace = T)
  mci_boot0$a <- F
  mci_boot <- rbind(mci_boot1, mci_boot0)
  beta_1 <- suppressWarnings(DR(Y = mci_boot$y, A = mci_boot$a, piX = as.matrix(mci_boot[,-(1:2)]), phiX = as.matrix(mci_boot[,-(1:2)])))$betaDR
  beta_0 <- suppressWarnings(DR(Y = mci_boot$y, A = !mci_boot$a, piX = as.matrix(mci_boot[,-(1:2)]), phiX = as.matrix(mci_boot[,-(1:2)])))$betaDR
  beta_1 - beta_0
})
# beta_1 <- suppressWarnings(DR(Y = mci$y, A = mci$a, piX = as.matrix(mci[,-(1:2)]), phiX = as.matrix(mci[,-(1:2)])))$betaDR
# beta_0 <- suppressWarnings(DR(Y = mci$y, A = !mci$a, piX = as.matrix(mci[,-(1:2)]), phiX = as.matrix(mci[,-(1:2)])))$betaDR
bias_rotnitzky <- mean(boot_rotnitzky)# - (beta_1 - beta_0)
var_rotnitzky <- var(boot_rotnitzky)

# SuperLearner
SL.library = c('SL.gam', 'SL.randomForest','SL.nnet', 'SL.svm')
cl <- makeCluster(16)
registerDoParallel(cl)
boot_superlearner <- foreach(i = 1:5000, .combine = cbind, .packages = c('SuperLearner', 'dplyr')) %dopar% {
  mci_boot1 <- sample_n(mci, size = n1, replace = T)
  mci_boot1$a <- T
  mci_boot0 <- sample_n(mci, size = n0, replace = T)
  mci_boot0$a <- F
  mci_boot <- rbind(mci_boot1, mci_boot0)
  SLfit_boot <- predict(SuperLearner(Y = mci_boot$y, X = mci_boot[,-(1:2)], SL.library = SL.library))
  ensembled <- mean(mci_boot$y[mci_boot$a]-SLfit_boot$pred[mci_boot$a]) - mean(mci_boot$y[!mci_boot$a]-SLfit_boot$pred[!mci_boot$a])
  c(ensembled=ensembled, mean(mci_boot$y[mci_boot$a])- colMeans(SLfit_boot$library.predict[mci_boot$a,]) - mean(mci_boot$y[!mci_boot$a]) + colMeans(SLfit_boot$library.predict[!mci_boot$a,]))
}
stopCluster(cl)
# SLfit_boot <- predict(SuperLearner(Y = mci$y, X = mci[,-(1:2)], SL.library = SL.library))
# ensembled <- mean(mci$y[mci$a]-SLfit_boot$pred[mci$a]) - mean(mci$y[!mci$a] - SLfit_boot$pred[!mci$a])
bias_superlearner <- colMeans(t(boot_superlearner)) #- c(ensembled, mean(mci$y[mci$a]) - colMeans(SLfit_boot$library.predict[mci$a,]) - mean(mci$y[!mci$a]) + colMeans(SLfit_boot$library.predict[!mci$a,]))
var_superlearner <- map_dbl(data.frame(t(boot_superlearner)), var)

# boot_superlearner <- map(1:10, function(i){
#   mci_boot <- rbind(sample_frac(mci[mci$a,], replace = T), sample_frac(mci[!mci$a,], replace = T))
#   SLfit_boot <- predict(SuperLearner(Y = mci_boot$y, X = mci_boot[,-(1:2)], SL.library = c('SL.gam', 'SL.randomForest','SL.nnet', 'SL.glmnet', 'SL.gbm', 'SL.svm')))
#   ensembled <- mean(mci_boot$y[mci_boot$a]-SLfit_boot$pred[mci_boot$a]) - mean(mci_boot$y[!mci_boot$a]-SLfit_boot$pred[!mci_boot$a])
#   c(ensembled=ensembled, mean(mci_boot$y[mci_boot$a])- colMeans(SLfit_boot$library.predict[mci_boot$a,]) - mean(mci_boot$y[!mci_boot$a]) + colMeans(SLfit_boot$library.predict[!mci_boot$a,]))
# })
# boot_superlearner <- matrix(unlist(boot_superlearner), ncol = 10)
# SLfit_boot <- predict(SuperLearner(Y = mci$y, X = mci[,-(1:2)], SL.library = c('SL.gam', 'SL.randomForest','SL.nnet', 'SL.glmnet', 'SL.gbm', 'SL.svm')))
# ensemble <- mean(mci$y[mci$a]-SLfit_boot$pred[mci$a]) - mean(mci$y[!mci$a] - SLfit_boot$pred[!mci$a])
# bias_superlearner <- colMeans(boot_superlearner) - c(ensemble, mean(mci$y[mci$a]) - colMeans(SLfit_boot$library.predict[mci$a,]) - mean(mci$y[!mci$a]) + colMeans(SLfit_boot$library.predict[!mci$a,]))
# var_superlearner <- map_dbl(boot_superlearner, var)


# result
MCI_result <- data.frame(bias = c(bias_unadj,NA,NA, bias_rotnitzky, bias_superlearner),
                         var = c(var_unadj,var_ancova_baseline, var_ancova_all, var_rotnitzky, var_superlearner),
                         RE = c(1,var_unadj/var_ancova_baseline, var_unadj/var_ancova_all, var_unadj/var_rotnitzky, var_unadj/var_superlearner))
rownames(MCI_result) <- c("unadjusted", "ANCOVA(all)", "ANCOVA(baseline score)", "Rotnitzky", names(bias_superlearner))
round(MCI_result,3)
