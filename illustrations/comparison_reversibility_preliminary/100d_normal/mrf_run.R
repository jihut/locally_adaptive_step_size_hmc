rm(list = ls())

library(ggplot2)
library(doParallel)
source("implementation_scripts/standard_hmc/standard_hmc_functions.R")

grad_log_target_fun <- function(theta) {
  as.numeric(-theta)
}

hamiltonian_func <- function(theta, rho) {
  # -mvtnorm::dmvnorm(theta, mean = mu, sigma = sigma, log = T) + 0.5 * sum(rho ^ 2)
  0.5 * sum(theta ^ 2) + 0.5 * sum(rho ^ 2)
}

log_target_fun <- function(theta) {
  -0.5 * sum(theta ^ 2)
}

d <- 100 # 100-dimensional example
h <- pi / 2
L <- 1
delta <- 1.22
max_c <- 20
n_warmup_iterations <- 100000
n_sampling_iterations <- 100000
n_chains <- 10

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

final_run <- foreach::foreach(i = 1:n_chains) %dopar% {
  set.seed(i)
  theta <- rnorm(d)
  single_warmup_run <- adaptive_step_size_standard_HMC(
    micro_fun = micro_fun_mrf_error,
    n_samples = n_warmup_iterations, 
    h = h, 
    L = L, 
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = theta,
    return_c_for_L_equal_to_1 = TRUE
  )
  single_sampling_run <- adaptive_step_size_standard_HMC(
    micro_fun = micro_fun_mrf_error,
    n_samples = n_sampling_iterations, 
    h = h, 
    L = L,
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = single_warmup_run$samples_matrix[nrow(single_warmup_run$samples_matrix), ],
    norm = "L-infinity"
  ) 
}

final_result <- # concatenate same-named elements together across the chains
  lapply(names(final_run[[1]]), function(element_name) {
    elements <- lapply(final_run, '[[', element_name)
    if (is.matrix(elements[[1]])) {
      do.call(rbind, elements)
    } else {
      unlist(elements)
    }
  })
names(final_result) <- names(final_run[[1]])

cov(final_result$samples_matrix[, 1:3])
colMeans(final_result$samples_matrix)
table(final_result$reversibility_indicator)
table(final_result$max_c_vector)
table(final_result$min_c_vector)
table(as.numeric(final_result$c_matrix))
table(final_result$accept_indicator)
sum(final_result$n_evals_ode)
sum(final_result$n_evals_ode) / 1e6

final_samples_3d_array <- # concatenate same-named elements together across the chains
  lapply(names(final_run[[1]]), function(element_name) {
    elements <- lapply(final_run, '[[', element_name)
    if (is.matrix(elements[[1]])) {
      do.call(rbind, elements)
    } else {
      unlist(elements)
    }
  })

final_samples_3d_array <- abind::abind(lapply(final_run, '[[', "samples_matrix"), along = 3)
final_samples_3d_array <- aperm(final_samples_3d_array, perm = c(1, 3, 2))
dim(final_samples_3d_array)

rstan_monitor_summary <- rstan::monitor(final_samples_3d_array, warmup = 0)
rstan_monitor_summary$n_eff
rstan_monitor_summary$n_eff * 1e4 / sum(final_result$n_evals_ode)

min_ess_index <- which.min(rstan_monitor_summary$n_eff)
min_ess_index
round(mean(final_result$samples_matrix[, min_ess_index]), 4)
round(sd(final_result$samples_matrix[, min_ess_index]), 4)
min(rstan_monitor_summary$n_eff)
min(rstan_monitor_summary$n_eff) * 1e4 / sum(final_result$n_evals_ode) 

max_ess_index <- which.max(rstan_monitor_summary$n_eff)
max_ess_index
round(mean(final_result$samples_matrix[, max_ess_index]), 4)
round(sd(final_result$samples_matrix[, max_ess_index]), 4)
max(rstan_monitor_summary$n_eff)
max(rstan_monitor_summary$n_eff) * 1e4 / sum(final_result$n_evals_ode) 
