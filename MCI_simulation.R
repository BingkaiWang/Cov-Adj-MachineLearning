set.seed(100)
library(dplyr)
library(purrr)
library(boot)

load("ACDS.rdata")
source(file = "Rotnitzky_estimator.R")
d <- subset(d, (d$arm != "Vitamin E") & (!is.na(d$Y18)))
y <- d$Y18 - d$Y0
a <- d$arm == "Donepezil"
w <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
mci <- cbind(y,a,w)

# unadjusted estimator
boot_unadj <- boot(mci, statistic = function(data, i){
  data = data[i,]
  mean(data$y[data$a]) - mean(data$y[!data$a])
}, R = 1000, stype = 'i')
bias_unadj <- mean(boot_unadj$t) - boot_unadj$t0
var_unadj <- var(boot_unadj$t)

# ANCOVA with all variable (LOOCV)
pred_y <- map_dbl(1:nrow(mci), function(i){
  predict(lm(y~., data = mci[-i,]), mci[i,])
})
var_ancova_all <- 4*mean((mci$y - pred_y)^2)/nrow(mci)

# ANCOVA with only baseline score (LOOCV)
pred_y <- map_dbl(1:nrow(mci), function(i){
  predict(lm(y~a+Y0, data = mci[-i,]), mci[i,])
})
var_ancova_baseline <- 4*mean((mci$y - pred_y)^2)/nrow(mci)

# Rotnitzky's estimator
boot_rotnitzky <- boot(mci, statistic = function(data, i){
  data = data[i,]
  beta_1 <- suppressWarnings(DR(Y = data$y, A = data$a, piX = as.matrix(data[,-(1:2)]), phiX = as.matrix(data[,-(1:2)])))$betaDR
  beta_0 <- suppressWarnings(DR(Y = data$y, A = !data$a, piX = as.matrix(data[,-(1:2)]), phiX = as.matrix(data[,-(1:2)])))$betaDR
  beta_1 - beta_0
}, R = 1000, stype = 'i')
bias_rotnitzky <- mean(boot_rotnitzky$t) - boot_rotnitzky$t0
var_rotnitzky <- var(boot_rotnitzky$t)

# result
MCI_result <- data.frame(bias = c(bias_unadj,NA,NA, bias_rotnitzky),
                         var = c(var_unadj,var_ancova_baseline, var_ancova_all, var_rotnitzky),
                         RE = c(1,var_unadj/var_ancova_baseline, var_unadj/var_ancova_all, var_unadj/var_rotnitzky))
rownames(MCI_result) <- c("unadjusted", "ANCOVA(all)", "ANCOVA(baseline score)", "Rotnitzky")
round(MCI_result,2)
