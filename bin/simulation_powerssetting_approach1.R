#' This code is to perform the power simulation for CRF project.
#' Three ML methods were utitlized, including glm, RF, XGBoost
#' 
#' @author Shi Yujia
#' @date 2024.11.07
#' 
#' 

suppressPackageStartupMessages({
    library(data.table)
    library(grf)
    library(xgboost)
    library(pwr)

    source("/mnt/md1/yujia/project/2024-11-06_CRF_assignment/power_simulation/bin/simulation_powers_function.R")
})

rsq <- function (x, y) cor(x, y) ^ 2

compute_stats <- function(tau.pred) {
  
  #' This function is to compute zval, pval and ajusted.p using the tau prediction
  #' 
  #' @param tau.pred: the prediction return by grf::predict
  #' 
  #' @return: a dataframe containing four columns, including original tau, tau z value, tau p-value, tau adjusted p-value
  
  tau.zval <- tau.pred$predictions/sqrt(tau.pred$variance.estimates) 
  tau.pval <- pnorm(abs(tau.zval), lower.tail = FALSE) * 2
  tau.p.adjust <- p.adjust(tau.pval, method = "BH")
  stats <- cbind(tau.pred$predictions, tau.zval, tau.pval, tau.p.adjust)
  colnames(stats) <- c("tau", "tau.zval", "tau.pval", "tau.p.adjust")
  return(as.data.frame(stats))

}

mu.mu <- c(function.8, function.5, function.4, function.7, function.3, function.1, function.2, function.6)
tau.tau <- c(function.1, function.2, function.3, function.4, function.5, function.6, function.7, function.8)

set.seed(1)
n_list <- c(60, 120, 180, 240, 300, 360, 420, 480) 
power_corr_mat.rf <- matrix(0, nrow = 8, ncol = 8)
power_corr_mat.xgb <- matrix(0, nrow = 8, ncol = 8)
power_corr_mat.glm <- matrix(0, nrow = 8, ncol = 8)

for (f in 1:8){
    
    message(paste0("Currently we are running scenario ", f, "."))
    power_corr_c.rf <- c()
    power_corr_c.xgb <- c()
    power_corr_c.glm <- c()

    for (n in n_list){

        message(paste0("Currently we are running n at ", n, "."))
        tau_corr_list.rf <- c()
        tau_corr_list.xgb <- c()
        tau_corr_list.glm <- c()

        for (k in 1:200){
            
            # Generate data.
            p <- 20 # We assume 10 binary and 10 continuous covariates

            # simulate the x
            dat.x <- as.data.frame(matrix(nrow = n, ncol = p))

            odd.cols <- colnames(dat.x[, c(TRUE, FALSE)])
            even.cols <- colnames(dat.x[, c(FALSE, TRUE)])

            # simulate the odd columns
            for (i in odd.cols){
                dat.x[, i] <- rnorm(n = n, sd = 1)
            }

            # simulate the even columns
            for (i in even.cols){
                dat.x[, i] <- rbinom(n = n, size = 1, prob = 0.5)
            }

            W <- rbinom(n, 1, 0.5)

            mu <- mu.mu[f][[1]]
            tau <- tau.tau[f][[1]]

            # simulate the outcome
            mean.y <- mu(dat.x) + (W - 0.5) * tau(dat.x)
            outcome.y  <- rnorm(n, mean = mean.y, sd = 0.25)

            # true tau
            true_tau <- tau(dat.x)
            true_meany <- mu(dat.x)

            # Approach 1: Follow the NHB paper. (Quite easy to get a very high power)
            X_regress <- cbind(dat.x, W)

            # Method 1: Random Forest
            r.forest <- regression_forest(X_regress[1:(n%/%2), ], outcome.y[1:(n%/%2)], num.trees = 100)
            pred_test <- predict(r.forest, X_regress[((n%/%2)+1):n, ], estimate.variance = TRUE)$predictions

            tt_rf <- cor.test(pred_test, outcome.y[((n%/%2)+1):n]) # Use same function as the NHB paper.
            tau_corr_list.rf <- c(tau_corr_list.rf, tt_rf$p.value)

            # Method 2: XGBoost
            r.xgb <- xgboost(data = as.matrix(X_regress[1:(n%/%2), ]), label = outcome.y[1:(n%/%2)], nrounds = 50, verbose = F)
            pred_test <- predict(r.xgb, as.matrix(X_regress[((n%/%2)+1):n, ]))

            tt.xgb <- cor.test(pred_test, outcome.y[((n%/%2)+1):n]) # Use same function as the NHB paper.
            tau_corr_list.xgb <- c(tau_corr_list.xgb, tt.xgb$p.value)

            # Method 3: glm
            X_glm <- data.frame(outcome.y, X_regress)
            colnames(X_glm)[1] <- "Y"
            r.glm <- glm(Y ~ ., family = gaussian, data = X_glm[1:(n%/%2), ])
            pred_test <- predict(r.glm, X_regress[((n%/%2)+1):n, ])

            tt.glm <- cor.test(pred_test, outcome.y[((n%/%2)+1):n]) # Use same function as the NHB paper.
            tau_corr_list.glm <- c(tau_corr_list.glm, tt.glm$p.value)

            # message("\n=============\n")
            # break
        }

        # break
        power_corr_c.rf <- c(power_corr_c.rf, length(tau_corr_list.rf[tau_corr_list.rf < 0.05])/200)
        power_corr_c.xgb <- c(power_corr_c.xgb, length(tau_corr_list.xgb[tau_corr_list.xgb < 0.05])/200)
        power_corr_c.glm <- c(power_corr_c.glm, length(tau_corr_list.glm[tau_corr_list.glm < 0.05])/200)
        
    }
    # break
    power_corr_mat.rf[f,] <- power_corr_c.rf
    power_corr_mat.xgb[f,] <- power_corr_c.xgb
    power_corr_mat.glm[f,] <- power_corr_c.glm
    
}

power_corr_mat.rf <- as.data.frame(power_corr_mat.rf)
power_corr_mat.xgb <- as.data.frame(power_corr_mat.xgb)
power_corr_mat.glm <- as.data.frame(power_corr_mat.glm)
fwrite(power_corr_mat.rf, "/mnt/md1/yujia/project/2024-11-06_CRF_assignment/power_simulation/res/test_Y_corr_powers_res_rf.csv")
fwrite(power_corr_mat.xgb, "/mnt/md1/yujia/project/2024-11-06_CRF_assignment/power_simulation/res/test_Y_corr_powers_res_xgb.csv")
fwrite(power_corr_mat.glm, "/mnt/md1/yujia/project/2024-11-06_CRF_assignment/power_simulation/res/test_Y_corr_powers_res_glm.csv")
