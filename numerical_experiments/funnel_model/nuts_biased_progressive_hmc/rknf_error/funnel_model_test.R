rm(list = ls())
library(doParallel)
library(ggplot2)
source("implementation_scripts/nuts_biased_progressive_hmc/stored_orbit_tuning_nuts_biased_progressive_hmc_functions.R")
source("numerical_experiments/funnel_model/general_scripts/funnel_10d_model.R")

# theta <- c(1, 1)
# rho <- c(0.25, 1)
h <- 0.25 # macro step size of 5 for now
m <- 7
delta <- 0.031 # energy tolerance
max_c <- 20
n_warmup_iterations <- 100000
n_sampling_iterations <- 100000
n_chains <- 10

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

final_run <- foreach::foreach(i = 1:n_chains) %dopar% {
  set.seed(i)
  theta <- numeric(d)
  theta[1] <- rnorm(1, mean = 0, sd = 3)
  theta[2:d] <- rnorm(d - 1, mean = 0, sd = exp(theta[1] / 2))
  sink(paste0("numerical_experiments/funnel_model/nuts_biased_progressive_hmc/rknf_error/log/log_nr_", i, ".txt"))
  print("Warmup")
  single_warmup_run <- adaptive_step_size_nuts_biased_progressive_HMC_sampling(
    micro_fun = micro_fun_rknf_error,
    n_samples = n_warmup_iterations, 
    h = h, 
    m = m, 
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = theta
  )
  print("Sampling")
  single_sampling_run <- adaptive_step_size_nuts_biased_progressive_HMC_sampling(
    micro_fun = micro_fun_rknf_error,
    n_samples = n_sampling_iterations, 
    h = h, 
    m = m, 
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = single_warmup_run$samples_matrix[n_warmup_iterations, ]
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
colMeans(final_result$samples_matrix)
table(final_result$dead_forward_vector)
table(final_result$dead_backward_vector)
table(pmax(final_result$dead_forward_vector, final_result$dead_backward_vector))
table(final_result$max_c_vector)
table(final_result$min_c_vector)
plot(final_result$max_orbit_energy_error_vector)
mean(final_result$max_orbit_energy_error_vector <= 0.2)
sum(final_result$n_evals_ode)
table(final_result$final_sample_index_vector)
sum(pmax(final_result$dead_forward_vector, final_result$dead_backward_vector) == 1) / sum(final_result$n_evals_ode)
final_c_vector <- as.numeric(final_result$c_matrix)
final_table_c <- table(final_c_vector)
final_c_df <- data.frame(final_table_c)
final_c_df
final_c_df$Freq[final_c_df$final_c_vector == 0] / 
  sum(final_c_df$Freq[final_c_df$final_c_vector != -1])

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

sink("numerical_experiments/funnel_model/nuts_biased_progressive_hmc/rknf_error/funnel_model_test.txt")
print("dead end")
print(table(pmax(final_result$dead_forward_vector, final_result$dead_backward_vector)))
print("ratio c")
print(final_c_df$Freq[final_c_df$final_c_vector == 0] / 
        sum(final_c_df$Freq[final_c_df$final_c_vector != -1]))
print("ratio energy")
print(mean(final_result$max_orbit_energy_error_vector <= 0.2))
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
print("theta1")
print(round(mean(final_result$samples_matrix[, 1]), 4))
print(round(sd(final_result$samples_matrix[, 1]), 4))
print((rstan_monitor_summary$n_eff[1]))
print((rstan_monitor_summary$n_eff)[1] * 1e6 / sum(final_result$n_evals_ode)) 
sink()
