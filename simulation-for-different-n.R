set.seed(123)
library(SuperLearner) # for running machine learning algorithms
library(tidyverse)
library(foreach) # for parallel programming
library(doParallel) # for parallel programming
cl <- makeCluster(16)
registerDoParallel(cl)

#' @param Y Continuous outcome
#' @param A Binary treatment indicator
#' @param W data frame of covariates
#' @param pi probability of being treated
#' @param n_fold number of folds for cross validation
#'
#' @return a matrix with 4 columns and 8 rows. 
#' Rows represent results from the unadjusted estimator, ANCOVA, lasso, ridge, 
#' GAM, rpart, random forest and an ensemble of these methods using SuperLearner.
#' Columns represent: point estimate of ATE without CV (first column) and with 5-fold CV (second column),
#' sample variance of point estimates without CV (third column) and with 5-fold CV (fourth column).
#' @export
#'
#' @examples
#' n <- 200
#' A <- rbinom(n, 1, 0.5)
#' W <- data.frame(matrix(rnorm(n*5), nrow = n))
#' Y <- A * (W[,1]-W[,2])^2 + (1-A) * (W[,3] + W[,4]) + rnorm(n)
#' cov_adj_variance(Y, A, W)
cov_adj_variance <- function(Y, A, W, pi = 0.5, n_fold = 5){
  n <- length(Y)
  test_size <- floor(n/n_fold)
  SL.library <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ridge", "SL.gam", "SL.rpart", "SL.randomForest")
  var_table <- matrix(NA, nrow = 8, ncol = 4)
  rownames(var_table) <- c(SL.library, "SL.ensemble")
  colnames(var_table) <- c("non-CV", "CV", "median-var-non-CV", "median-var-CV")
  
  # non-CV variance
  SLfit.1 <- SuperLearner(Y[A==1], W[A==1,], newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
  SLfit.0 <- SuperLearner(Y[A==0], W[A==0,], newX = W, family = "gaussian", SL.library = SL.library, cvControl = list(V=5))
  
  for(i in 1:length(SL.library)){
    var_table[i, 1] <- mean(A*Y/pi - (1-A)*Y/pi - (A-pi)*(SLfit.1$library.predict[,i]/pi + SLfit.0$library.predict[,i]/(1-pi)))
    var_table[i, 3] <- var(A*Y/pi - (1-A)*Y/pi - (A-pi)*(SLfit.1$library.predict[,i]/pi + SLfit.0$library.predict[,i]/(1-pi)))
    
  }
  var_table[8, 1] <- mean(A*Y/pi - (1-A)*Y/pi - (A-pi)*(SLfit.1$library.predict[,i]/pi + SLfit.0$library.predict[,i]/(1-pi)))
  var_table[8, 3] <- var(A * Y / pi  - (1-A) * Y/pi  - (A - pi) * (SLfit.1$SL.predict/pi + SLfit.0$SL.predict/(1-pi)))
  
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
  var_table[,2] <- apply(combined_cv_resutls, 2, mean)
  var_table[,4] <- apply(combined_cv_resutls, 2, var)
  return(var_table)
}

n_sim <- 500
SL.library <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ridge", "SL.gam", "SL.rpart", "SL.randomForest")
# MCI simulation
load("ACDS.rdata")
d <- subset(d,  (!is.na(d$Y18)))
Y <- d$Y18 - d$Y0
W <- select(d, age, female, cototscr, mmscore, adtotscr, gdstot, Y0)

MCI_sim_results <- matrix(NA, nrow = 32, ncol = 4)
colnames(MCI_sim_results) <- c("200", "500", "1000", "2000")
rownames(MCI_sim_results) <- rep(c(SL.library, "ensemble"), 4)

n <- 200
MCI_sim_200 <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','SuperLearner')) %dopar% {
  sample_indi <- sample(1:length(Y), size = n, replace = T)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  MCI_indirect <- cov_adj_variance(Y_sample, A_sample, W_sample)
  as.vector(MCI_indirect)
}
MCI_sim_results[1:16,1] <- apply(MCI_sim_200[1:16,], 1, var)
MCI_sim_results[17:32,1] <- apply(MCI_sim_200[17:32,], 1, median)

n <- 500
MCI_sim_500 <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','SuperLearner')) %dopar% {
  sample_indi <- sample(1:length(Y), size = n, replace = T)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  MCI_indirect <- cov_adj_variance(Y_sample, A_sample, W_sample)
  as.vector(MCI_indirect)
}
MCI_sim_results[1:16,2] <- apply(MCI_sim_500[1:16,], 1, var)
MCI_sim_results[17:32,2] <- apply(MCI_sim_500[17:32,], 1, median)

saveRDS(MCI_sim_results, file = "simulation-for-different-n.rds")
print("n = 500 finished.")

n <- 1000
MCI_sim_1000 <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','SuperLearner')) %dopar% {
  sample_indi <- sample(1:length(Y), size = n, replace = T)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  MCI_indirect <- cov_adj_variance(Y_sample, A_sample, W_sample)
  as.vector(MCI_indirect)
}
MCI_sim_results[1:16,3] <- apply(MCI_sim_1000[1:16,], 1, var)
MCI_sim_results[17:32,3] <- apply(MCI_sim_1000[17:32,], 1, median)

saveRDS(MCI_sim_results, file = "simulation-for-different-n.rds")
print("n = 1000 finished.")

n <- 2000
MCI_sim_2000 <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','SuperLearner')) %dopar% {
  sample_indi <- sample(1:length(Y), size = n, replace = T)
  Y_sample <- Y[sample_indi]
  W_sample <- W[sample_indi, ]
  A_sample <- rbinom(n, size = 1, prob = 0.5)
  MCI_indirect <- cov_adj_variance(Y_sample, A_sample, W_sample)
  as.vector(MCI_indirect)
}
MCI_sim_results[1:16,4] <- apply(MCI_sim_2000[1:16,], 1, var)
MCI_sim_results[17:32,4] <- apply(MCI_sim_2000[17:32,], 1, median)

saveRDS(MCI_sim_results, file = "simulation-for-different-n.rds")
print("n = 2000 finished.")


#### presenting results -----------
library(cowplot)
MCI_sim_results <- readRDS("simulation-for-different-n.rds")
mean_no_CV <- data.frame(MCI_sim_results[1:8,]) %>% mutate(method = rownames(MCI_sim_results)[1:8])
colnames(mean_no_CV) <- c("200", "500", "1000", "2000", "Methods")
for(i in 8:1){mean_no_CV[i,1:4] <- mean_no_CV[1,1:4]/mean_no_CV[i,1:4]}
mean_no_CV <- gather(mean_no_CV, "n", "RE", 1:4)
mean_no_CV$n <- factor(mean_no_CV$n, levels = c("200", "500", "1000", "2000"))
p_mean_no_CV <- ggplot(mean_no_CV) +
  geom_line(aes(x = n, y = RE, linetype = Methods, color = Methods, group = Methods), size = 1.5) +
  geom_hline(yintercept = 1) +
  theme_bw() + 
  theme(text = element_text(size=12),
        legend.text = element_text(size=11)) +
  xlab("Resampling size (with replacement)") + ylim(c(0,10.5)) +
  ggtitle("RE (sample variance) without CV")

mean_CV <- data.frame(MCI_sim_results[9:16,]) %>% mutate(method = rownames(MCI_sim_results)[1:8])
colnames(mean_CV) <- c("200", "500", "1000", "2000", "Methods")
for(i in 8:1){mean_CV[i,1:4] <- mean_CV[1,1:4]/mean_CV[i,1:4]}
mean_CV <- gather(mean_CV, "n", "RE", 1:4)
mean_CV$n <- factor(mean_CV$n, levels = c("200", "500", "1000", "2000"))
p_mean_CV <- ggplot(mean_CV) +
  geom_line(aes(x = n, y = RE, linetype = Methods, color = Methods, group = Methods), size = 1.5) +
  geom_hline(yintercept = 1) +
  theme_bw() + 
  theme(text = element_text(size=12),
        legend.text = element_text(size=11)) +
  xlab("Resampling size (with replacement)") + ylim(c(0,10.5)) +
  ggtitle("RE (sample variance) with 5-fold CV")

p_summary <- plot_grid(p_mean_no_CV, p_mean_CV, ncol = 2)
save_plot(plot = p_summary, filename = "RE-different-n-1.png", base_aspect_ratio = 3)

median_no_CV <- data.frame(MCI_sim_results[17:24,]) %>% mutate(method = rownames(MCI_sim_results)[1:8])
colnames(median_no_CV) <- c("200", "500", "1000", "2000", "Methods")
for(i in 8:1){median_no_CV[i,1:4] <- median_no_CV[1,1:4]/median_no_CV[i,1:4]}
median_no_CV <- gather(median_no_CV, "n", "RE", 1:4)
median_no_CV$n <- factor(median_no_CV$n, levels = c("200", "500", "1000", "2000"))
p_median_no_CV <- ggplot(median_no_CV) +
  geom_line(aes(x = n, y = RE, linetype = Methods, color = Methods, group = Methods), size = 1.5) +
  geom_hline(yintercept = 1) +
  theme_bw() + 
  theme(text = element_text(size=12),
        legend.text = element_text(size=11)) +
  xlab("Resampling size (with replacement)") + ylim(c(0,10.5)) +
  ggtitle("RE (median variance) without CV")

median_CV <- data.frame(MCI_sim_results[25:32,]) %>% mutate(method = rownames(MCI_sim_results)[1:8])
colnames(median_CV) <- c("200", "500", "1000", "2000", "Methods")
for(i in 8:1){median_CV[i,1:4] <- median_CV[1,1:4]/median_CV[i,1:4]}
median_CV <- gather(median_CV, "n", "RE", 1:4)
median_CV$n <- factor(median_CV$n, levels = c("200", "500", "1000", "2000"))
p_median_CV <- ggplot(median_CV) +
  geom_line(aes(x = n, y = RE, linetype = Methods, color = Methods, group = Methods), size = 1.5) +
  geom_hline(yintercept = 1) +
  theme_bw() + 
  theme(text = element_text(size=12),
        legend.text = element_text(size=11)) +
  xlab("Resampling size (with replacement)") + ylim(c(0,10.5)) +
  ggtitle("RE (median variance) with 5-fold CV")

p_summary <- plot_grid(p_median_no_CV, p_median_CV, ncol = 2)
save_plot(plot = p_summary, filename = "RE-different-n-2.png", base_aspect_ratio = 3)

