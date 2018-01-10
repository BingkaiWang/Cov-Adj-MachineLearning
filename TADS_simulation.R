set.seed(300)
library(dplyr)
library(purrr)
library(boot)

load("TADS.rdata")
source(file = "Rotnitzky_estimator.R")

sim_tad <- function(data, treatment_arm){
  data <- data[data$treatment %in% treatment_arm, ]
  y <- data$change_score
  a <- data$treatment == treatment_arm[1]
  # w <- select(data, age, gender, CDRS_baseline, CGI, CGAS, RADS, suicide_ideation, depression_episode, comorbidity)
  w <- select(data, age, gender, CDRS_baseline, CGI, CGAS,RADS, depression_episode)
  tads <- cbind(y, a, w)
  # unadjusted estimator
  boot_unadj <- boot(tads, statistic = function(data, i){
    data = data[i,]
    mean(data$y[data$a]) - mean(data$y[!data$a])
  }, R = 1000, stype = 'i')
  bias_unadj <- mean(boot_unadj$t) - boot_unadj$t0
  var_unadj <- var(boot_unadj$t)
  
  # ANCOVA with all variable (LOOCV)
  pred_y <- map_dbl(1:nrow(tads), function(i){
    predict(lm(y~., data = tads[-i,]), tads[i,])
  })
  var_ancova_all <- 4*mean((tads$y - pred_y)^2)/nrow(tads)
  
  # ANCOVA with only baseline score (LOOCV)
  pred_y <- map_dbl(1:nrow(tads), function(i){
    predict(lm(y~a+CDRS_baseline, data = tads[-i,]), tads[i,])
  })
  var_ancova_baseline <- 4*mean((tads$y - pred_y)^2)/nrow(tads)
  
  # Rotnitzky's estimator
  boot_rotnitzky <- boot(tads, statistic = function(data, i){
    data = data[i,]
    beta_1 <- DR(Y = data$y, A = data$a, piX = as.matrix(data[,-(1:2)]), phiX = as.matrix(data[,-(1:2)]))$betaDR
    beta_0 <- suppressWarnings(DR(Y = data$y, A = !data$a, piX = as.matrix(data[,-(1:2)]), phiX = as.matrix(data[,-(1:2)])))$betaDR
    beta_1 - beta_0
  }, R = 1000, stype = 'i')
  bias_rotnitzky <- mean(boot_rotnitzky$t) - boot_rotnitzky$t0
  var_rotnitzky <- var(boot_rotnitzky$t)
  
  # result
  tads_result <- data.frame(bias = c(bias_unadj,NA,NA, bias_rotnitzky),
                            var = c(var_unadj,var_ancova_baseline, var_ancova_all, var_rotnitzky),
                            RE = c(1,var_unadj/var_ancova_baseline, var_unadj/var_ancova_all, var_unadj/var_rotnitzky))
  return(tads_result)
}

sim_tad(data = subset(tad, !is.na(tad$CDRS_12)), treatment_arm = c("CBT", "PBO")) %>% round(2)
sim_tad(data = subset(tad, !is.na(tad$CDRS_12)), treatment_arm = c("FLX", "PBO")) %>% round(2)
sim_tad(data = subset(tad, !is.na(tad$CDRS_12)), treatment_arm = c("COMB", "PBO")) %>% round(2)
