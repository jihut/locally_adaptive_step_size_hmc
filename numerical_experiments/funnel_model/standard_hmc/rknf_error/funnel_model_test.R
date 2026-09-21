rm(list = ls())
library(doParallel)
library(ggplot2)
source("implementation_scripts/standard_hmc/standard_hmc_functions.R")
source("numerical_experiments/funnel_model/general_scripts/funnel_10d_model.R")

# h <- 0.25 
# delta <- 1.1 # this is better than 1.0
# L <- 16
h <- 0.5
L <- 4
delta <- 0.4
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
  dir.create("numerical_experiments/funnel_model/standard_hmc/rknf_error/log", showWarnings = F)
  sink(paste0("numerical_experiments/funnel_model/standard_hmc/rknf_error/log/log_nr_", i, ".txt"))
  print("Warmup")
  single_warmup_run <- adaptive_step_size_standard_HMC(
    micro_fun = micro_fun_rknf_error,
    n_samples = n_warmup_iterations, 
    h = h, 
    L = L,
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = theta,
    norm = "L-infinity"
  ) 
  print("Sampling")
  single_sampling_run <- adaptive_step_size_standard_HMC(
    micro_fun = micro_fun_rknf_error,
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
  sink()
  single_sampling_run
}

parallel::stopCluster(init_cluster)
saveRDS(final_run, "numerical_experiments/funnel_model/standard_hmc/rknf_error/funnel_model_test.RDS")


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
final_c_table <- table(as.numeric(final_result$c_matrix))
final_c_df <- data.frame(final_c_table)
colnames(final_c_df) <- c("final_c_vector", "Freq")
sum(final_result$n_evals_ode)

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

sink("numerical_experiments/funnel_model/standard_hmc/rknf_error/funnel_model_test.txt")
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
print("theta1")
print(round(mean(final_result$samples_matrix[, 1]), 4))
print(round(sd(final_result$samples_matrix[, 1]), 4))
print((rstan_monitor_summary$n_eff[1]))
print((rstan_monitor_summary$n_eff)[1] * 1e6 / sum(final_result$n_evals_ode)) 
sink()


min_ess_index <- which.min(rstan_monitor_summary$n_eff)
min_ess_index
round(mean(final_result$samples_matrix[, min_ess_index]), 4)
round(sd(final_result$samples_matrix[, min_ess_index]), 4)
min(rstan_monitor_summary$n_eff)
min(rstan_monitor_summary$n_eff) * 1e6 / sum(final_result$n_evals_ode) 

max_ess_index <- which.max(rstan_monitor_summary$n_eff)
max_ess_index
round(mean(final_result$samples_matrix[, max_ess_index]), 4)
round(sd(final_result$samples_matrix[, max_ess_index]), 4)
max(rstan_monitor_summary$n_eff)
max(rstan_monitor_summary$n_eff) * 1e6 / sum(final_result$n_evals_ode) 

max_c_table <- table(final_result$max_c_vector)
max_c_count_df <- data.frame(c = names(max_c_table), counts = as.numeric(max_c_table), 
                             ratio = as.numeric(max_c_table / sum(max_c_table)))
max_c_count_df$c <- factor(max_c_count_df$c, levels = max_c_count_df$c)
max_c_count_df
ggplot(max_c_count_df, aes(x = c, y = ratio, color = c)) + 
  geom_bar(stat = "identity", fill = "white") + 
  ylim(
    c(0, 1)
  )

min_c_table <- table(final_result$min_c_vector)
min_c_count_df <- data.frame(c = names(min_c_table), counts = as.numeric(min_c_table), 
                             ratio = as.numeric(min_c_table / sum(min_c_table)))
min_c_count_df$c <- factor(min_c_count_df$c, levels = min_c_count_df$c)
min_c_count_df
ggplot(min_c_count_df, aes(x = c, y = ratio, color = c)) + 
  geom_bar(stat = "identity", fill = "white") + 
  ylim(
    c(0, 1)
  )
