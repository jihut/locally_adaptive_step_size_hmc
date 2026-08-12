rm(list = ls())
library(doParallel)
library(ggplot2)
source("implementation_scripts/kinetic_langevin_ghmc/tuning_one_step_GHMC_functions.R")

ar1_model <- bridgestan::StanModel$new(
  "numerical_experiments/ar1_model/general_scripts/ar1.stan",
  "numerical_experiments/ar1_model/general_scripts/ar1_data.json",
  seed = 123
)

grad_log_target_fun <- function(theta) {
  ar1_model$log_density_gradient(theta)$gradient
}

hamiltonian_func <- function(theta, rho) {
  -ar1_model$log_density(theta, propto = FALSE) + 0.5 * sum(rho ^ 2)
} 

log_target_fun <- function(theta) ar1_model$log_density(theta, propto = FALSE)


cor_parameter <- 0.95
h <- 1 # macro step size of 5 for now
delta <- 0.06 # for m = 4
# m <- 5
# delta <- 0.089 # for m = 5
max_c <- 20
n_warmup_iterations <- 100000
n_sampling_iterations <- 100000
n_chains <- 10
nu <- 0.95

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

final_run <- foreach::foreach(i = 1:n_chains) %dopar% {
  ar1_model <- bridgestan::StanModel$new(
    "numerical_experiments/ar1_model/general_scripts/ar1.stan",
    "numerical_experiments/ar1_model/general_scripts/ar1_data.json",
    seed = 123
  )
  
  grad_log_target_fun <- function(theta) {
    ar1_model$log_density_gradient(theta)$gradient
  }
  
  hamiltonian_func <- function(theta, rho) {
    -ar1_model$log_density(theta, propto = FALSE) + 0.5 * sum(rho ^ 2)
  }
  
  log_target_fun <- function(theta) ar1_model$log_density(theta, propto = FALSE)
  
  set.seed(i)
  # theta <- rnorm(100)
  theta <- numeric(100)
  theta[1] <- rnorm(1)
  for (j in 2:100) {
    theta[j] <- cor_parameter * theta[j - 1] + rnorm(1, sd = sqrt(1 - cor_parameter ^ 2))
  }
  sink(paste0("numerical_experiments/ar1_model/kinetic_langevin_ghmc/rknf_error/log/log_nr_", i, ".txt"))
  print("Warmup")
  single_warmup_run <- adaptive_step_size_one_step_GHMC_sampling(
    micro_fun = micro_fun_rknf_error,
    n_samples = n_warmup_iterations, 
    h = h, 
    nu = nu,
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = theta,
    norm = "L-infinity"
  ) 
  print("Sampling")
  single_sampling_run <- adaptive_step_size_one_step_GHMC_sampling(
    micro_fun = micro_fun_rknf_error,
    n_samples = n_warmup_iterations, 
    h = h, 
    nu = nu,
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = single_warmup_run$samples_matrix[nrow(single_warmup_run$samples_matrix), ],
    norm = "L-infinity"
  ) 
  sink()
  single_sampling_run
}

parallel::stopCluster(init_cluster)

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
final_result

cov(final_result$samples_matrix)
cov(final_result$samples_matrix[, 1:3])
colMeans(final_result$samples_matrix)
table(final_result$reversibility_indicator)
table(final_result$accept_indicator)
table(final_result$c_vector)
sum(final_result$n_evals_ode)
final_c_table <- table(final_result$c_vector)
final_c_df <- data.frame(final_c_table)
colnames(final_c_df) <- c("final_c_vector", "Freq")

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
rstan_monitor_summary$n_eff * 1000000 / sum(final_result$n_evals_ode)

sink("numerical_experiments/ar1_model/kinetic_langevin_ghmc/rknf_error/ar1_model_test.txt")
print("reversibility")
print(table(final_result$reversibility_indicator))
print("metropolis acceptance")
print(table(final_result$accept_indicator))
print("ratio c")
print(final_c_df$Freq[final_c_df$final_c_vector == 0] / 
        sum(final_c_df$Freq[final_c_df$final_c_vector != -1]))
print("n evals ode")
print(sum(final_result$n_evals_ode))
print("min ess case")
min_ess_index <- which.min(rstan_monitor_summary$n_eff)
print(min_ess_index)
print(round(mean(final_result$samples_matrix[, min_ess_index]), 4))
print(round(sd(final_result$samples_matrix[, min_ess_index]), 4))
print(min(rstan_monitor_summary$n_eff))
print(min(rstan_monitor_summary$n_eff) * 1e6 / sum(final_result$n_evals_ode)) 
print("max ess case")
max_ess_index <- which.max(rstan_monitor_summary$n_eff)
print(max_ess_index)
print(round(mean(final_result$samples_matrix[, max_ess_index]), 4))
print(round(sd(final_result$samples_matrix[, max_ess_index]), 4))
print(max(rstan_monitor_summary$n_eff))
print(max(rstan_monitor_summary$n_eff) * 1e6 / sum(final_result$n_evals_ode)) 
sink()
