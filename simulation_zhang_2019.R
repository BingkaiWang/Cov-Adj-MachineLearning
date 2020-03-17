set.seed(123)
library(SuperLearner) # for running machine learning algorithms
library(tidyverse)
library(xtable) # for creating latex code of the summary table

#' @param Y Continuous outcome
#' @param A Binary treatment indicator
#' @param W data frame of covariates
#' @param pi probability of being treated
#' @param n_fold number of folds for cross validation
#'
#' @return a matrix with 2 columns and 8 rows. 
#' Rows represent results from the unadjusted estimator, ANCOVA, lasso, ridge, 
#' GAM, rpart, random forest and an ensemble of these methods using SuperLearner.
#' Columns represent: sample variance of point estimates without CV (first column) 
#' and with 5-fold CV (second column) using the "indirect method" descibed in Zhang et al. (2019).
#' @export
#'
#' @examples
#' n <- 200
#' A <- rbinom(n, 1, 0.5)
#' W <- data.frame(matrix(rnorm(n*5), nrow = n))
#' Y <- A * (W[,1]-W[,2])^2 + (1-A) * (W[,3] + W[,4]) + rnorm(n)
#' cov_adj_variance(Y, A, W)
cov_adj_zhang_indirect <- function(Y, A, W, pi = 0.5, n_fold = 5){
  n <- length(Y)
  test_size <- floor(n/n_fold)
  SL.library <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ridge", "SL.gam", "SL.rpart", "SL.randomForest")
  var_table <- matrix(NA, nrow = 8, ncol = 2)
  rownames(var_table) <- c(SL.library, "SL.ensemble")
  colnames(var_table) <- c("non-CV", "CV")
  
  # non-CV variance
  SLfit.1 <- SuperLearner(Y[A==1], W[A==1,], newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
  SLfit.0 <- SuperLearner(Y[A==0], W[A==0,], newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
  
  for(i in 1:length(SL.library)){
    var_table[i, 1] <-  var(A*Y/pi - (1-A)*Y/pi - (A-pi)*(SLfit.1$library.predict[,i]/pi + SLfit.0$library.predict[,i]/(1-pi)))
  }
  var_table[8, 1] <- var(A * Y / pi  - (1-A) * Y/pi  - (A - pi) * (SLfit.1$SL.predict/pi + SLfit.0$SL.predict/(1-pi)))
  
  # CVed variance
  random_order <- sample(1:n, replace = F)
  random_partition <- map(1:n_fold, ~random_order[(.-1)*test_size + 1:test_size])
  random_partition[[n_fold]] <- random_order[((n_fold-1)*test_size + 1):n]
  cv_results <- map(1:n_fold, function(k){
    Y_train <- Y[-random_partition[[k]]]
    A_train <- A[-random_partition[[k]]]
    W_train <- W[-random_partition[[k]],]
    Y_test <- Y[random_partition[[k]]]
    A_test <- A[random_partition[[k]]]
    W_test <- W[random_partition[[k]], ]
    SLfit.1 <- SuperLearner(Y_train[A_train==1], W_train[A_train==1,], newX = W_test, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
    SLfit.0 <- SuperLearner(Y_train[A_train==0], W_train[A_train==0,], newX = W_test, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
    cbind(
      (A_test*Y_test/pi - (1- A_test)*Y_test/pi) %*% t(rep(1,length(SL.library))) - (A_test-pi)*(SLfit.1$library.predict/pi + SLfit.0$library.predict/(1-pi)),
      A_test * Y_test / pi  - (1- A_test) * Y_test/pi  - (A_test - pi) * (SLfit.1$SL.predict/pi + SLfit.0$SL.predict/(1-pi))
    )
  })
  combined_cv_resutls <- NULL
  for(k in 1:n_fold){
    combined_cv_resutls <- rbind(combined_cv_resutls, cv_results[[k]])
  }
  var_table[,2] <- apply(combined_cv_resutls, 2, var)
  return(var_table)
}

#' @param Y Continuous outcome
#' @param A Binary treatment indicator
#' @param W data frame of covariates
#' @param pi probability of being treated
#' @param n_fold number of folds for cross validation
#'
#' @return a matrix with 2 columns and 8 rows. 
#' Rows represent results from the unadjusted estimator, ANCOVA, lasso, ridge, 
#' GAM, rpart, random forest and an ensemble of these methods using SuperLearner.
#' Columns represent: sample variance of point estimates without CV (first column) 
#' and with 5-fold CV (second column) using the "direct method" descibed in Zhang et al. (2019).
#' @export
#'
#' @examples
#' n <- 200
#' A <- rbinom(n, 1, 0.5)
#' W <- data.frame(matrix(rnorm(n*5), nrow = n))
#' Y <- A * (W[,1]-W[,2])^2 + (1-A) * (W[,3] + W[,4]) + rnorm(n)
#' cov_adj_variance(Y, A, W)
cov_adj_zhang_direct <- function(Y, A, W, n_fold = 5){
  # psi(A,Y)/(A-pi) = 4*Y
  pi <- 0.5
  n <- length(Y)
  test_size <- floor(n/n_fold)
  SL.library <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ridge", "SL.gam", "SL.rpart", "SL.randomForest")
  var_table <- matrix(NA, nrow = 8, ncol = 2)
  rownames(var_table) <- c(SL.library, "SL.ensemble")
  colnames(var_table) <- c("non-CV", "CV")
  
  # non-CV variance
  SLfit <- SuperLearner (4 * Y, W, newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
  for(i in 1:length(SL.library)){
    var_table[i, 1] <-  var(A*Y/pi - (1-A)*Y/pi - (A-pi)* SLfit$library.predict[,i])
  }
  var_table[8, 1] <- var(A * Y / pi  - (1-A) * Y/pi  - (A - pi) * SLfit$SL.predict)
  
  # CVed variance
  random_order <- sample(1:n, replace = F)
  random_partition <- map(1:n_fold, ~random_order[(.-1)*test_size + 1:test_size])
  random_partition[[n_fold]] <- random_order[((n_fold-1)*test_size + 1):n]
  cv_results <- map(1:n_fold, function(k){
    Y_train <- Y[-random_partition[[k]]]
    A_train <- A[-random_partition[[k]]]
    W_train <- W[-random_partition[[k]],]
    Y_test <- Y[random_partition[[k]]]
    A_test <- A[random_partition[[k]]]
    W_test <- W[random_partition[[k]], ]
    SLfit <- SuperLearner (4 * Y_train, W_train, newX = W_test, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
    cbind(
      (A_test*Y_test/pi - (1- A_test)*Y_test/pi) %*% t(rep(1,length(SL.library))) - (A_test-pi)*SLfit$library.predict,
      A_test * Y_test / pi  - (1- A_test) * Y_test/pi  - (A_test - pi) * SLfit$SL.predict
    )
  })
  combined_cv_resutls <- NULL
  for(k in 1:n_fold){
    combined_cv_resutls <- rbind(combined_cv_resutls, cv_results[[k]])
  }
  var_table[,2] <- apply(combined_cv_resutls, 2, var)
  return(var_table)
}

# data analysis ----------------
# MCI (Donepezil vs Placebo)
load("ACDS.rdata")
d <- subset(d, (d$arm != "Vitamin E") & (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
A <- d$arm == "Donepezil"
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
MCI_Donepezil_indirect <- cov_adj_zhang_indirect(Y, A, W)
MCI_Donepezil_direct <- cov_adj_zhang_direct(Y, A, W)
MCI_Donepezil <- cbind(MCI_Donepezil_indirect, MCI_Donepezil_direct)

# MCI (Vitamin E vs Placebo)
load("ACDS.rdata")
d <- subset(d, (d$arm != "Donepezil") & (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
A <- d$arm == "Vitamin E"
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
MCI_VE_indirect <- cov_adj_zhang_indirect(Y, A, W)
MCI_VE_direct <- cov_adj_zhang_direct(Y, A, W)
MCI_VE <- cbind(MCI_VE_indirect, MCI_VE_direct)

# CLEAR III
clearIII <- readxl::read_xlsx(path = "CLEAR_III_master_file_9.26.2018.xlsx")
clearIII <- clearIII %>%
  select("study_arm","glasgow_rankin_365", "randomization_gcs_Total",
         "stabct_ich_volume_rc", "randomization_nihss_total",
         "er_present_nihss_total", "stabct_ivh_volume_rc")
clearIII$randomization_nihss_total[which(is.na(clearIII$randomization_nihss_total))] <- median(clearIII$randomization_nihss_total, na.rm = T)
clearIII$er_present_nihss_total[which(is.na(clearIII$er_present_nihss_total))] <- median(clearIII$er_present_nihss_total, na.rm = T)
clearIII$randomization_gcs_Total <- as.numeric(clearIII$randomization_gcs_Total)
clearIII$randomization_gcs_Total[which(is.na(clearIII$randomization_gcs_Total))] <- median(clearIII$randomization_gcs_Total, na.rm = T)
clearIII <- filter(clearIII, !is.na(glasgow_rankin_365))
Y <- clearIII$glasgow_rankin_365
W <- select(clearIII, randomization_gcs_Total, stabct_ich_volume_rc, randomization_nihss_total,
            er_present_nihss_total, stabct_ivh_volume_rc)
A <- clearIII$study_arm
clear3_indirect <- cov_adj_zhang_indirect(Y, A, W)
clear3_direct <- cov_adj_zhang_direct(Y, A, W)
clear3 <- cbind(clear3_indirect, clear3_direct)

saveRDS(object = rbind(MCI_Donepezil, MCI_VE, clear3), file = "zhang_data_analysis.rds")
# xtable(rbind(MCI_Donepezil, MCI_VE, clear3))

# simulation 1 (zhang et al. Sample with replacement) ------------
n <- 200
n_sim <- 500
# MCI simulation
load("ACDS.rdata")
d <- subset(d,  (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
MCI_sim <- map(1:n_sim, function(j){
  sample_indi <- sample(1:length(Y), size = n, replace = T)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  MCI_indirect <- cov_adj_zhang_indirect(Y_sample, A_sample, W_sample)
  MCI_direct <- cov_adj_zhang_direct(Y_sample, A_sample, W_sample)
  cbind(MCI_indirect, MCI_direct)
})
MCI_sim_result <- MCI_sim[[1]]
for(i in 1:length(MCI_sim_result)){
  MCI_sim_result[i] <- median(map_dbl(MCI_sim, ~.[i]))
}

# clearIII simulation
clearIII <- readxl::read_xlsx(path = "CLEAR_III_master_file_9.26.2018.xlsx")
clearIII <- clearIII %>%
  select("study_arm","glasgow_rankin_365", "randomization_gcs_Total",
         "stabct_ich_volume_rc", "randomization_nihss_total",
         "er_present_nihss_total", "stabct_ivh_volume_rc")
clearIII$randomization_nihss_total[which(is.na(clearIII$randomization_nihss_total))] <- median(clearIII$randomization_nihss_total, na.rm = T)
clearIII$er_present_nihss_total[which(is.na(clearIII$er_present_nihss_total))] <- median(clearIII$er_present_nihss_total, na.rm = T)
clearIII$randomization_gcs_Total <- as.numeric(clearIII$randomization_gcs_Total)
clearIII$randomization_gcs_Total[which(is.na(clearIII$randomization_gcs_Total))] <- median(clearIII$randomization_gcs_Total, na.rm = T)
clearIII <- filter(clearIII, !is.na(glasgow_rankin_365))
Y <- clearIII$glasgow_rankin_365
W <- select(clearIII, randomization_gcs_Total, stabct_ich_volume_rc, randomization_nihss_total,
            er_present_nihss_total, stabct_ivh_volume_rc)
clear3_sim <- map(1:n_sim, function(j){
  sample_indi <- sample(1:length(Y), size = n, replace = T)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  clear3_indirect <- cov_adj_zhang_indirect(Y_sample, A_sample, W_sample)
  clear3_direct <- cov_adj_zhang_direct(Y_sample, A_sample, W_sample)
  clear3 <- cbind(clear3_indirect, clear3_direct)
})
clear3_sim_result <- clear3_sim[[1]]
for(i in 1:length(clear3_sim_result)){
  clear3_sim_result[i] <- median(map_dbl(clear3_sim, ~.[i]))
}

# TADS
load("TADS.rdata")
tad <- subset(tad, !is.na(tad$CDRS_12))
Y <- tad$change_score
W <- select(tad, age, gender, CDRS_baseline, CGI, CGAS,RADS, depression_episode)
TADS_sim <- map(1:n_sim, function(j){
  sample_indi <- sample(1:length(Y), size = n, replace = T)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  TADS_indirect <- cov_adj_zhang_indirect(Y_sample, A_sample, W_sample)
  TADS_direct <- cov_adj_zhang_direct(Y_sample, A_sample, W_sample)
  TADS <- cbind(TADS_indirect, TADS_direct)
})
TADS_sim_result <- TADS_sim[[1]]
for(i in 1:length(TADS_sim_result)){
  TADS_sim_result[i] <- median(map_dbl(TADS_sim, ~.[i]))
}

saveRDS(object = rbind(MCI_sim_result, clear3_sim_result, TADS_sim_result), file = "zhang_sim_analysis.rds")


# simulation 2: resample without replacement -----------
n <- 200
n_sim <- 500
# MCI simulation
load("ACDS.rdata")
d <- subset(d,  (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
MCI_sim <- map(1:n_sim, function(j){
  sample_indi <- sample(1:length(Y), size = n, replace = F)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  MCI_indirect <- cov_adj_zhang_indirect(Y_sample, A_sample, W_sample)
  MCI_direct <- cov_adj_zhang_direct(Y_sample, A_sample, W_sample)
  MCI <- cbind(MCI_indirect, MCI_direct)
})
MCI_sim_result <- MCI_sim[[1]]
for(i in 1:length(MCI_sim_result)){
  MCI_sim_result[i] <- median(map_dbl(MCI_sim, ~.[i]))
}

# clearIII simulation 
clearIII <- readxl::read_xlsx(path = "CLEAR_III_master_file_9.26.2018.xlsx")
clearIII <- clearIII %>% 
  select("study_arm","glasgow_rankin_365", "randomization_gcs_Total", 
         "stabct_ich_volume_rc", "randomization_nihss_total",
         "er_present_nihss_total", "stabct_ivh_volume_rc")
clearIII$randomization_nihss_total[which(is.na(clearIII$randomization_nihss_total))] <- median(clearIII$randomization_nihss_total, na.rm = T)
clearIII$er_present_nihss_total[which(is.na(clearIII$er_present_nihss_total))] <- median(clearIII$er_present_nihss_total, na.rm = T)
clearIII$randomization_gcs_Total <- as.numeric(clearIII$randomization_gcs_Total)
clearIII$randomization_gcs_Total[which(is.na(clearIII$randomization_gcs_Total))] <- median(clearIII$randomization_gcs_Total, na.rm = T)
clearIII <- filter(clearIII, !is.na(glasgow_rankin_365))
Y <- clearIII$glasgow_rankin_365
W <- select(clearIII, randomization_gcs_Total, stabct_ich_volume_rc, randomization_nihss_total,
            er_present_nihss_total, stabct_ivh_volume_rc)
clear3_sim <- map(1:n_sim, function(j){
  sample_indi <- sample(1:length(Y), size = n, replace = F)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  clear3_indirect <- cov_adj_zhang_indirect(Y_sample, A_sample, W_sample)
  clear3_direct <- cov_adj_zhang_direct(Y_sample, A_sample, W_sample)
  clear3 <- cbind(clear3_indirect, clear3_direct)
})
clear3_sim_result <- clear3_sim[[1]]
for(i in 1:length(clear3_sim_result)){
  clear3_sim_result[i] <- median(map_dbl(clear3_sim, ~.[i]))
}

# TADS
load("TADS.rdata")
tad <- subset(tad, !is.na(tad$CDRS_12))
Y <- tad$change_score
W <- select(tad, age, gender, CDRS_baseline, CGI, CGAS,RADS, depression_episode)
TADS_sim <- map(1:n_sim, function(j){
  sample_indi <- sample(1:length(Y), size = n, replace = F)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  TADS_indirect <- cov_adj_zhang_indirect(Y_sample, A_sample, W_sample)
  TADS_direct <- cov_adj_zhang_direct(Y_sample, A_sample, W_sample)
  TADS <- cbind(TADS_indirect, TADS_direct)
})
TADS_sim_result <- TADS_sim[[1]]
for(i in 1:length(TADS_sim_result)){
  TADS_sim_result[i] <- median(map_dbl(TADS_sim, ~.[i]))
}

saveRDS(object = rbind(MCI_sim_result, clear3_sim_result, TADS_sim_result), file = "zhang_sim_noreplacement.rds")



##### presenting results ------------
sim_results <- readRDS("zhang_data_analysis.rds")
for(i in 2:8){sim_results[i,] <- sim_results[1,]/sim_results[i,]}
sim_results[1,] <- rep(1, 4)
for(i in 10:16){sim_results[i,] <- sim_results[9,]/sim_results[i,]}
sim_results[9,] <- rep(1, 4)
for(i in 18:24){sim_results[i,] <- sim_results[17,]/sim_results[i,]}
sim_results[17,] <- rep(1, 4)
round(sim_results,2)
xtable(sim_results)


sim_results <- readRDS("zhang_sim_analysis.rds")
for(i in 2:8){sim_results[i,] <- sim_results[1,]/sim_results[i,]}
sim_results[1,] <- rep(1, 4)
for(i in 10:16){sim_results[i,] <- sim_results[9,]/sim_results[i,]}
sim_results[9,] <- rep(1, 4)
for(i in 18:24){sim_results[i,] <- sim_results[17,]/sim_results[i,]}
sim_results[17,] <- rep(1, 4)
round(sim_results,2)
xtable(sim_results)


sim_results <- readRDS("zhang_sim_noreplacement.rds")
for(i in 2:8){sim_results[i,] <- sim_results[1,]/sim_results[i,]}
sim_results[1,] <- rep(1, 4)
for(i in 10:16){sim_results[i,] <- sim_results[9,]/sim_results[i,]}
sim_results[9,] <- rep(1, 4)
for(i in 18:24){sim_results[i,] <- sim_results[17,]/sim_results[i,]}
sim_results[17,] <- rep(1, 4)
round(sim_results,2)
xtable(sim_results)
