#' This code is to generate the simulation data in instrumental setting.
#' The main idea comes from the power's paper "Some methods for heterogeneous treatment effect
#' estimation in high dimensions"
#' Also, Kai's thesis will be used as reference
#' 
#' @author: Yujia Shi
#' @date: 2022.01.07

#' ============================================================
#' 8 functions to generate the mean effect and treatment effect
#' ============================================================

# tau <- function.8
# p <- 20 # We assume 10 binary and 10 continuous covariates
# n <- 50000
# # simulate the x
# dat.x <- as.data.frame(matrix(nrow = n, ncol = p))

# odd.cols <- colnames(dat.x[, c(TRUE, FALSE)])
# even.cols <- colnames(dat.x[, c(FALSE, TRUE)])

# # simulate the odd columns
# for (i in odd.cols){
#     dat.x[, i] <- rnorm(n = n, sd = 1)
# }

# # simulate the even columns
# for (i in even.cols){
#     dat.x[, i] <- rbinom(n = n, size = 1, prob = 0.5)
# }

# mean(tau(dat.x))

function.1 <- function(x){
    res <- rep(0, dim(x)[1])
    res
}

function.2 <- function(x){
    res <- 0.35 * ifelse(x[, 1] > 1, 1, 0) 
    res
}

function.3 <- function(x){
    res <- 0.145 * (2 * x[, 1] - 4)
    res
}

function.4 <- function(x){
    res <- x[, 2] * x[, 4] * x[, 6] + 2 * x[, 2] * x[, 4] * (1 - x[, 6]) + 
            3 * x[, 2] * (1 - x[, 4]) * x[, 6] +
            4 * x[, 2] * (1 - x[, 4]) * (1 - x[, 6]) + 
            5 * (1 - x[, 2]) * x[, 4] * x[, 6] + 
            6 * (1 - x[, 2]) * x[, 4] * (1 - x[, 6]) +
            7 * (1 - x[, 2]) * (1 - x[, 4]) * x[, 6] + 
            8 * (1 - x[, 2]) * (1 - x[, 4]) * (1 - x[, 6])
    res*0.125
}

function.5 <- function(x){
    res <- x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2
    res*0.38
}

function.6 <- function(x){
    res <- 4 * ifelse(x[, 1] > 1, 1, 0) * ifelse(x[, 3] > 0, 1, 0) +
            4 * ifelse(x[, 5] > 1, 1, 0) * ifelse(x[, 7] > 0, 1, 0) +
            2 * x[, 8] * x[, 9]
    res * 0.9
}

function.7 <- function(x){
    res <- x[, 1] ** 2 + x[, 2] + x[, 3] ** 2 + 
            x[, 4] + x[, 5] ** 2 + x[, 6] + x[, 7] ** 2 + x[, 8] + x[, 9] ** 2 - 11
    res <- 0.5 * res
    res * 0.285
}

function.8 <- function(x){
    res <- 0.5 * (function.4(x) + function.6(x))
}
