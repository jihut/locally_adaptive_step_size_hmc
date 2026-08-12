rm(list = ls())
library(doParallel)
library(ggplot2)
source("implementation_scripts/multinomial_hmc/tuning_multinomial_HMC_functions.R")
source("numerical_experiments/funnel_model/general_scripts/funnel_10d_model.R")

# theta <- c(1, 1)
# rho <- c(0.25, 1)
set.seed(1)
h <- 0.25 # macro step size of 5 for now
L <- 16
delta <- 0.1 # energy tolerance
max_c <- 20
# moving_avg_parameter <- 10
prob_no_step_size_halving <- 0.8
prob_max_orbit_energy_error_below_threshold <- 0.8
max_orbit_energy_error_threshold <- 0.2
n_warmup_iterations <- 100000
n_warmup_iterations_between_updates <- 1000
n_sampling_iterations <- 100000
n_chains <- 10

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

final_run <- foreach::foreach(i = 1:n_chains) %dopar% {
  set.seed(i)
  theta <- numeric(d)
  theta[1] <- rnorm(1, mean = 0, sd = 3)
  theta[2:d] <- rnorm(d - 1, mean = 0, sd = exp(theta[1] / 2))
  sink(paste0("numerical_experiments/funnel_model/multinomial_hmc/mrf_error/log/log_nr_", i, ".txt"))
  print("Warmup")
  single_warmup_run <- adaptive_step_size_multinomial_HMC_warmup(
    micro_fun = micro_fun_mrf_error,
    n_samples = n_warmup_iterations, 
    # n_samples_between_updates = n_warmup_iterations_between_updates, 
    # moving_avg_indicator = FALSE,
    # moving_avg_parameter = moving_avg_parameter, 
    h = h, 
    L = L, 
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = theta,
    norm = "L-infinity",
    prob_no_step_size_halving = prob_no_step_size_halving,
    prob_max_orbit_energy_error_below_threshold = prob_max_orbit_energy_error_below_threshold,
    max_orbit_energy_error_threshold = max_orbit_energy_error_threshold
  )
  print("Sampling")
  single_sampling_run <- adaptive_step_size_multinomial_HMC_sampling(
    micro_fun = micro_fun_mrf_error,
    n_samples = n_sampling_iterations, 
    h = single_warmup_run$final_h, 
    L = single_warmup_run$final_L, 
    delta = single_warmup_run$final_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = single_warmup_run$samples_matrix[n_warmup_iterations, ],
    norm = "L-infinity"
  )
  sink()
  list(single_warmup_run = single_warmup_run, single_sampling_run = single_sampling_run)
}

parallel::stopCluster(init_cluster)

# Warmup results 

final_warmup_run <- lapply(final_run, '[[', 1)

final_warmup_result <- # concatenate same-named elements together across the chains
  lapply(names(final_warmup_run[[1]]), function(element_name) {
    elements <- lapply(final_warmup_run, '[[', element_name)
    if (is.matrix(elements[[1]])) {
      if (element_name != "c_matrix") {
        do.call(rbind, elements)
      } else {
        as.matrix(
          dplyr::bind_rows(lapply(elements, as.data.frame))
        )
      }
    } else {
      unlist(elements)
    }
  })

names(final_warmup_result) <- names(final_warmup_run[[1]])
final_warmup_result
range_final_warmup_h <- range(final_warmup_result$final_h)
range_final_warmup_delta <- range(final_warmup_result$final_delta)

# Sampling results

final_sampling_run <- lapply(final_run, '[[', 2)

final_sampling_result <- # concatenate same-named elements together across the chains
  lapply(names(final_sampling_run[[1]]), function(element_name) {
    elements <- lapply(final_sampling_run, '[[', element_name)
    if (is.matrix(elements[[1]])) {
      if (element_name != "c_matrix") {
        do.call(rbind, elements)
      } else {
        as.matrix(
          dplyr::bind_rows(lapply(elements, as.data.frame))
        )
      }
    } else {
      unlist(elements)
    }
  })

names(final_sampling_result) <- names(final_sampling_run[[1]])
final_sampling_result

cov(final_sampling_result$samples_matrix)
colMeans(final_sampling_result$samples_matrix)
table(final_sampling_result$dead_forward_vector)
table(final_sampling_result$dead_backward_vector)
table(pmax(final_sampling_result$dead_forward_vector, final_sampling_result$dead_backward_vector))
table(final_sampling_result$max_c_vector)
table(final_sampling_result$min_c_vector)
plot(final_sampling_result$max_orbit_energy_error_vector)
mean(final_sampling_result$max_orbit_energy_error_vector <= 0.2)
final_sampling_c_vector <- as.numeric(final_sampling_result$c_matrix)
final_table_sampling_c <- table(final_sampling_c_vector)
final_sampling_c_df <- data.frame(final_table_sampling_c)
final_sampling_c_df
final_sampling_c_df$Freq[final_sampling_c_df$final_sampling_c_vector == 0] / 
  sum(final_sampling_c_df$Freq[final_sampling_c_df$final_sampling_c_vector != -1])
sum(final_sampling_result$n_evals_ode)
table(final_sampling_result$final_sample_index_vector)
hist(final_sampling_result$final_sample_index_vector)
sum(pmax(final_sampling_result$dead_forward_vector, final_sampling_result$dead_backward_vector) == 1) / sum(final_sampling_result$n_evals_ode)

final_samples_3d_array <- abind::abind(lapply(final_sampling_run, '[[', "samples_matrix"), along = 3)
final_samples_3d_array <- aperm(final_samples_3d_array, perm = c(1, 3, 2))
dim(final_samples_3d_array)

rstan_monitor_summary <- rstan::monitor(final_samples_3d_array, warmup = 0)
rstan_monitor_summary$n_eff
rstan_monitor_summary$n_eff * 1000000 / sum(final_sampling_result$n_evals_ode)

sink("numerical_experiments/funnel_model/multinomial_hmc/mrf_error/funnel_model_tuning_test.txt")
print("range of tuned h")
print(range_final_warmup_h)
print("range of tuned delta")
print(range_final_warmup_delta)
print("dead end")
print(table(pmax(final_sampling_result$dead_forward_vector, final_sampling_result$dead_backward_vector)))
print("ratio c")
print(final_sampling_c_df$Freq[final_sampling_c_df$final_sampling_c_vector == 0] / 
        sum(final_sampling_c_df$Freq[final_sampling_c_df$final_sampling_c_vector != -1]))
print("ratio energy")
print(mean(final_sampling_result$max_orbit_energy_error_vector <= 0.2))
print("n evals ode")
print(sum(final_sampling_result$n_evals_ode))
print("min ess case")
min_ess_index <- which.min(rstan_monitor_summary$n_eff)
print(min_ess_index)
print(round(mean(final_sampling_result$samples_matrix[, min_ess_index]), 4))
print(round(sd(final_sampling_result$samples_matrix[, min_ess_index]), 4))
print(min(rstan_monitor_summary$n_eff))
print(min(rstan_monitor_summary$n_eff) * 1e6 / sum(final_sampling_result$n_evals_ode)) 
print("max ess case")
max_ess_index <- which.max(rstan_monitor_summary$n_eff)
print(max_ess_index)
print(round(mean(final_sampling_result$samples_matrix[, max_ess_index]), 4))
print(round(sd(final_sampling_result$samples_matrix[, max_ess_index]), 4))
print(max(rstan_monitor_summary$n_eff))
print(max(rstan_monitor_summary$n_eff) * 1e6 / sum(final_sampling_result$n_evals_ode)) 
print("theta1")
print(round(mean(final_sampling_result$samples_matrix[, 1]), 4))
print(round(sd(final_sampling_result$samples_matrix[, 1]), 4))
print((rstan_monitor_summary$n_eff[1]))
print((rstan_monitor_summary$n_eff)[1] * 1e6 / sum(final_sampling_result$n_evals_ode)) 
sink()
