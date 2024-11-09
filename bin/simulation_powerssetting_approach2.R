#' This code is to perform the power simulation for CRF project.
#' Here, we test the most naive tau function, defined as 0.42 + 0.28 * X2 + 0.3 * X1
#' 
#' @author Shi Yujia
#' @date 2024.11.07
#' 
#' 

suppressPackageStartupMessages({
    library(data.table)
    library(grf)
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
power_corr_mat <- matrix(0, nrow = 8, ncol = 8)
power_testcalibra_mat <- matrix(0, nrow = 8, ncol = 8)

for (f in 1:8){
    
    message(paste0("Currently we are running scenario ", f, "."))
    power_testcalibra_c <- c()
    power_corr_c <- c()

    for (n in n_list){

        message(paste0("Currently we are running n at ", n, "."))
        test_calibration_res_list <- c()
        tau_corr_list <- c()

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
            
            # Approach 2: 
            causal.forest <- causal_forest(dat.x, outcome.y, W, W.hat = rep(0.5, length(W)))
            tau_pred <- predict(causal.forest, estimate.variance=TRUE)
            tau_pred <- compute_stats(tau_pred)
            
            cor_res <- cor.test(true_tau, tau_pred$tau)
            tau_corr_list <- c(tau_corr_list, cor_res$p.value)

            test_calibration_res <- test_calibration(causal.forest)
            test_calibration_res_list <- c(test_calibration_res_list, test_calibration_res[2,4])
        }

        power_testcalibra_c <- c(power_testcalibra_c, length(test_calibration_res_list[test_calibration_res_list < 0.05])/200)
        power_corr_c <- c(power_corr_c, length(tau_corr_list[tau_corr_list < 0.05])/200)
        
    }
    # break
    power_corr_mat[f,] <- power_corr_c
    power_testcalibra_mat[f, ] <- power_testcalibra_c
    
}

power_testcalibra_mat <- as.data.frame(power_testcalibra_mat)
power_corr_mat <- as.data.frame(power_corr_mat)

fwrite(power_testcalibra_mat, "/mnt/md1/yujia/project/2024-11-06_CRF_assignment/power_simulation/res/test_calibration_powers_res.csv")
fwrite(power_corr_mat, "/mnt/md1/yujia/project/2024-11-06_CRF_assignment/power_simulation/res/test_corr_powers_res.csv")
