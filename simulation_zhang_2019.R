set.seed(123)
library(SuperLearner) # for running machine learning algorithms
library(tidyverse)
library(xtable) # for creating latex code of the summary table
library(readstata13)
center <- function(x){x-mean(x)}

#' @param Y Continuous outcome
#' @param A Binary treatment indicator
#' @param W data frame of covariates
#' @param pi probability of being treated
#' @param n_fold number of folds for cross validation
#'
#' @return a matrix with 2 columns and 8 rows. 
#' Rows represent results from the unadjusted estimator, ANCOVA, lasso, ridge, 
#' GAM, rpart, random forest and an ensemble of these methods using SuperLearner.
#' Columns represent: point estimate of ATE (first column) and its asymptotic variance (second column) without CV,
#' point estimate of ATE (third column) and its asymptotic variance (fourth column) with 20-fold CV
#' using the "indirect method" descibed in Zhang et al. (2019).
#' @export
#'
#' @examples
#' n <- 200
#' A <- rbinom(n, 1, 0.5)
#' W <- data.frame(matrix(rnorm(n*5), nrow = n))
#' Y <- A * (W[,1]-W[,2])^2 + (1-A) * (W[,3] + W[,4]) + rnorm(n)
#' cov_adj_zhang_indirect(Y, A, W)
cov_adj_zhang_indirect <- function(Y, A, W, pi = 0.5, n_fold = 20){
  n <- length(Y)
  test_size <- floor(n/n_fold)
  SL.library <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ridge", "SL.gam", "SL.rpart", "SL.randomForest")
  var_table <- matrix(NA, nrow = 8, ncol = 4)
  rownames(var_table) <- c(SL.library, "SL.ensemble")
  colnames(var_table) <- c("est-no-CV", "var-no-CV", "est-CV", "var-CV")
  bar_A <- mean(A)
  psi_hat <- A * (Y - mean(Y[A==1]))/bar_A - (1-A) * (Y - mean(Y[A==0]))/(1-bar_A)
  
  # non-CV
  SLfit.1 <- SuperLearner(2 * Y[A==1], W[A==1,], newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
  SLfit.0 <- SuperLearner(2 * Y[A==0], W[A==0,], newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))

  for(i in 1:length(SL.library)){
    var_table[i, "est-no-CV"] <- mean(Y[A==1]) - mean(Y[A==0]) - mean((A-bar_A)*(SLfit.1$library.predict[,i] + SLfit.0$library.predict[,i]))
    var_table[i, "var-no-CV"] <-  var(psi_hat - (A-bar_A)*(center(SLfit.1$library.predict[,i] + SLfit.0$library.predict[,i])))
  }
  var_table[8, "est-no-CV"] <- mean(Y[A==1]) - mean(Y[A==0]) - mean((A-bar_A)*(SLfit.1$SL.predict + SLfit.0$SL.predict))
  var_table[8, "var-no-CV"] <- var(psi_hat  - (A-bar_A)*(center(SLfit.1$SL.predict + SLfit.0$SL.predict)))
  
  # CVed
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
    SLfit.1 <- SuperLearner(2 * Y_train[A_train==1], W_train[A_train==1,], newX = W_test, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
    SLfit.0 <- SuperLearner(2 * Y_train[A_train==0], W_train[A_train==0,], newX = W_test, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
    cbind(
      SLfit.1$library.predict + SLfit.0$library.predict, SLfit.1$SL.predict + SLfit.0$SL.predict
    )
  })
  combined_cv_resutls <- NULL
  for(k in 1:n_fold){
    combined_cv_resutls <- rbind(combined_cv_resutls, cv_results[[k]])
  }
  for(i in 1:8){
    var_table[i,"est-CV"] <- mean(Y[A==1]) - mean(Y[A==0]) - mean((A[random_order] - bar_A) * center(combined_cv_resutls[,i]))
    var_table[i,"var-CV"] <- var(psi_hat[random_order] - (A[random_order] - bar_A) * center(combined_cv_resutls[,i]))
  }
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
#' cov_adj_zhang_direct(Y - mean(Y), A, W)
cov_adj_zhang_direct <- function(Y, A, W, pi = 0.5, n_fold = 20){
  n <- length(Y)
  test_size <- floor(n/n_fold)
  SL.library <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ridge", "SL.gam", "SL.rpart", "SL.randomForest")
  var_table <- matrix(NA, nrow = 8, ncol = 4)
  rownames(var_table) <- c(SL.library, "SL.ensemble")
  colnames(var_table) <- c("est-no-CV", "var-no-CV", "est-CV", "var-CV")
  bar_A <- mean(A)
  psi_hat <- A * (Y - mean(Y[A==1]))/bar_A - (1-A) * (Y - mean(Y[A==0]))/(1-bar_A)
  
  # non-CV
  SLfit <- SuperLearner(4 * Y, W, newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
  for(i in 1:length(SL.library)){
    var_table[i, "est-no-CV"] <- mean(Y[A==1]) - mean(Y[A==0]) - mean((A-bar_A) * center(SLfit$library.predict[,i]))
    var_table[i, "var-no-CV"] <-  var(psi_hat - (A-bar_A) * center(SLfit$library.predict[,i]))
  }
  var_table[8, "est-no-CV"] <- mean(Y[A==1]) - mean(Y[A==0]) - mean((A-bar_A) * center(SLfit$SL.predict))
  var_table[8, "var-no-CV"] <- var(psi_hat  - (A-bar_A) * center(SLfit$SL.predict))
  
  # CVed
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
    SLfit <- SuperLearner(4 * Y_train, W_train, newX = W_test, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
    cbind(SLfit$library.predict, SLfit$SL.predict)
  })
  combined_cv_resutls <- NULL
  for(k in 1:n_fold){
    combined_cv_resutls <- rbind(combined_cv_resutls, cv_results[[k]])
  }
  for(i in 1:8){
    var_table[i,"est-CV"] <- mean(Y[A==1]) - mean(Y[A==0]) - mean((A[random_order] - bar_A) * center(combined_cv_resutls[,i]))
    var_table[i,"var-CV"] <- var(psi_hat[random_order] - (A[random_order] - bar_A) * center(combined_cv_resutls[,i]))
  }
  return(var_table)
}

# data analysis ----------------
# MCI (Donepezil vs Placebo)
load("data/ACDS.rdata")
d <- subset(d, (d$arm != "Vitamin E") & (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
A <- d$arm == "Donepezil"
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
MCI_Donepezil_indirect <- cov_adj_zhang_indirect(Y, A, W)
MCI_Donepezil_direct <- cov_adj_zhang_direct(Y, A, W)
MCI_Donepezil <- cbind(MCI_Donepezil_indirect, MCI_Donepezil_direct)

# MCI (Vitamin E vs Placebo)
load("data/ACDS.rdata")
d <- subset(d, (d$arm != "Donepezil") & (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
A <- d$arm == "Vitamin E"
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
MCI_VE_indirect <- cov_adj_zhang_indirect(Y, A, W)
MCI_VE_direct <- cov_adj_zhang_direct(Y, A, W)
MCI_VE <- cbind(MCI_VE_indirect, MCI_VE_direct)

# CLEAR III
clearIII <- readxl::read_xlsx(path = "data/CLEAR_III_master_file_9.26.2018.xlsx")
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

# TADS
load("data/TADS.rdata")
tad <- subset(tad, !is.na(tad$CDRS_12))
Y <- tad$change_score
A <- tad$treatment %in% c("FLX", "COMB")
W <- select(tad, age, gender, CDRS_baseline, CGI, CGAS,RADS, depression_episode)
TADS_indirect <- cov_adj_zhang_indirect(Y, A, W)
TADS_direct <- cov_adj_zhang_direct(Y, A, W)
TADS <- cbind(TADS_indirect, TADS_direct)


# ASI
merged_ASI <- read.dta13("data/merged_ASI.dta")
merged_ASI <- merged_ASI %>% select(studyid:race, reduc_stim, ASIpsy_B: stimulant_B) %>% filter(!is.na(reduc_stim))
Y <- merged_ASI$reduc_stim
A <- merged_ASI$arm == "Treatment"
W <- merged_ASI %>% select(sex:age, ASIpsy_B: stimulant_B)
W$sex <- W$sex == "Female"
for(j in 2:10){
  W[is.na(W[,j]),j] <- median(W[,j], na.rm = T)
}
ASI_indirect <- cov_adj_zhang_indirect(Y, A, W)
ASI_direct <- cov_adj_zhang_direct(Y, A, W)
ASI <- cbind(ASI_indirect, ASI_direct)

saveRDS(object = rbind(MCI_Donepezil, MCI_VE, clear3, TADS, ASI), file = "zhang_data_analysis.rds")

# simulation 1 (zhang et al. Sample with replacement) ------------
library(foreach) # for parallel programming
library(doParallel) # for parallel programming
cl <- makeCluster(16)
registerDoParallel(cl)
simulation_wrapper <- function(Y, W, pi = 0.5, n_fold = 5, resample_size = 200, n_sim = 1000, replacement = T){
  simulation_result <- foreach(i=1:n_sim, .packages=c('tidyverse','SuperLearner'), .export = c('center', 'cov_adj_zhang_indirect', 'cov_adj_zhang_direct')) %dopar% {
    sample_indi <- sample(1:length(Y), size = resample_size, replace = replacement)
    Y_sample <- Y[sample_indi]
    W_sample <- W[sample_indi, ]
    A_sample <- rbinom(resample_size, size = 1, prob = pi)
    indirect <- cov_adj_zhang_indirect(Y_sample, A_sample, W_sample, pi = pi, n_fold = n_fold)
    direct <- cov_adj_zhang_direct(Y_sample, A_sample, W_sample, pi = pi, n_fold = n_fold)
    cbind(indirect, direct)
  }
  sim_result <- matrix(NA, nrow = 8, ncol = 10)
  colnames(sim_result) <- rep(c("Bias", "SD", "RE", "SE", "SE(CV)"), 2)
  rownames(sim_result) <- rownames(simulation_result[[1]])
  for(i in 1:8){
    sim_result[i,1] <- mean(map_dbl(simulation_result, ~.[i,3]))
    sim_result[i,2] <- sd(map_dbl(simulation_result, ~.[i,3]))
    sim_result[i,3] <- (sim_result[1,2]/sim_result[i,2])^2
    sim_result[i,4] <- sqrt(median(map_dbl(simulation_result, ~.[i,2]))/200)
    sim_result[i,5] <- sqrt(median(map_dbl(simulation_result, ~.[i,4]))/200)
    sim_result[i,6] <- mean(map_dbl(simulation_result, ~.[i,7]))
    sim_result[i,7] <- sd(map_dbl(simulation_result, ~.[i,7]))
    sim_result[i,8] <- (sim_result[1,7]/sim_result[i,7])^2
    sim_result[i,9] <- sqrt(median(map_dbl(simulation_result, ~.[i,6]))/200)
    sim_result[i,10] <- sqrt(median(map_dbl(simulation_result, ~.[i,8]))/200)
  }
  return(sim_result)
}
# MCI simulation
load("data/ACDS.rdata")
d <- subset(d,  (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)
MCI_sim_result_with_replace <- simulation_wrapper(Y,W, replacement = T)
saveRDS(MCI_sim_result_with_replace, file = "MCI_sim_result_with_replace.rds")
MCI_sim_result_without_replace <- simulation_wrapper(Y,W, replacement = F)
saveRDS(MCI_sim_result_without_replace, file = "MCI_sim_result_without_replace.rds")


# clearIII simulation
clearIII <- readxl::read_xlsx(path = "data/CLEAR_III_master_file_9.26.2018.xlsx")
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
clear3_sim_result_with_replace <- simulation_wrapper(Y,W, replacement = T)
saveRDS(clear3_sim_result_with_replace, file = "clear3_sim_result_with_replace.rds")
clear3_sim_result_without_replace <- simulation_wrapper(Y,W, replacement = F)
saveRDS(clear3_sim_result_without_replace, file = "clear3_sim_result_without_replace.rds")

# TADS
load("data/TADS.rdata")
tad <- subset(tad, !is.na(tad$CDRS_12))
Y <- tad$change_score
W <- select(tad, age, gender, CDRS_baseline, CGI, CGAS,RADS, depression_episode)
TADS_sim_result_with_replace <- simulation_wrapper(Y,W, replacement = T)
saveRDS(TADS_sim_result_with_replace, file = "TADS_sim_result_with_replace.rds")
TADS_sim_result_without_replace <- simulation_wrapper(Y,W, replacement = F)
saveRDS(TADS_sim_result_without_replace, file = "TADS_sim_result_without_replace.rds")

# ASI
merged_ASI <- read.dta13("data/merged_ASI.dta")
merged_ASI <- merged_ASI %>% select(studyid:race, reduc_stim, ASIpsy_B: stimulant_B) %>% filter(!is.na(reduc_stim))
Y <- merged_ASI$reduc_stim
A <- merged_ASI$arm == "Treatment"
W <- merged_ASI %>% select(sex:age, ASIpsy_B: stimulant_B)
W$sex <- W$sex == "Female"
for(j in 2:10){
  W[is.na(W[,j]),j] <- median(W[,j], na.rm = T)
}
ASI_sim_result_with_replace <- simulation_wrapper(Y,W, replacement = T)
saveRDS(ASI_sim_result_with_replace, file = "ASI_sim_result_with_replace.rds")
ASI_sim_result_without_replace <- simulation_wrapper(Y,W, replacement = F)
saveRDS(ASI_sim_result_without_replace, file = "ASI_sim_result_without_replace.rds")


# ##### presenting results ------------
sim_results <- readRDS("zhang_data_analysis.rds") %>% .[,c(2,4,6,8)]
for(i in 2:8){sim_results[i,] <- sim_results[1,]/sim_results[i,]}
sim_results[1,] <- rep(1, 4)
for(i in 10:16){sim_results[i,] <- sim_results[9,]/sim_results[i,]}
sim_results[9,] <- rep(1, 4)
for(i in 18:24){sim_results[i,] <- sim_results[17,]/sim_results[i,]}
sim_results[17,] <- rep(1, 4)
for(i in 26:32){sim_results[i,] <- sim_results[25,]/sim_results[i,]}
sim_results[25,] <- rep(1, 4)
for(i in 34:40){sim_results[i,] <- sim_results[33,]/sim_results[i,]}
sim_results[33,] <- rep(1, 4)
round(sim_results,2)
xtable(sim_results)

# Simulation with replacement
readRDS("MCI_sim_result_with_replace.rds") %>% round(2) %>% xtable
readRDS("clear3_sim_result_with_replace.rds") %>% round(2) %>% xtable
readRDS("TADS_sim_result_with_replace.rds") %>% round(2) %>% xtable
readRDS("ASI_sim_result_with_replace.rds") %>% round(2) %>% xtable

# Simulation without replacement
readRDS("MCI_sim_result_without_replace.rds") %>% round(2) %>% xtable
readRDS("clear3_sim_result_without_replace.rds") %>% round(2) %>% xtable
readRDS("TADS_sim_result_without_replace.rds") %>% round(2) %>% xtable
readRDS("ASI_sim_result_without_replace.rds") %>% round(2) %>% xtable
