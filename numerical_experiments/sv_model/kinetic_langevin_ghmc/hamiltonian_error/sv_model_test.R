rm(list = ls())
library(doParallel)
library(ggplot2)
source("implementation_scripts/kinetic_langevin_ghmc/tuning_one_step_GHMC_functions.R")

sv_model <- bridgestan::StanModel$new(
  "numerical_experiments/sv_model/general_scripts/SVmodel.stan",
  "numerical_experiments/sv_model/general_scripts/SVmodel_data.json",
  seed = 123
)

d <- sv_model$param_num()

grad_log_target_fun <- function(theta) {
  # sv_model$log_density_gradient(theta)$gradient
  grad_eval <- tryCatch(
    sv_model$log_density_gradient(theta)$gradient,
    error = function(e) rep(1e12, d)
  )
  
  if (is.na(grad_eval[1])) {
    grad_eval <- rep(1e12, d)
  }
  
  grad_eval
}

hamiltonian_func <- function(theta, rho) {
  # -sv_model$log_density(theta, propto = FALSE) + 0.5 * sum(rho ^ 2)
  -tryCatch(
    sv_model$log_density(theta, propto = FALSE),
    error = function(e) -1e12
  ) + 
    0.5 * sum(rho ^ 2)
} 

log_target_fun <- function(theta) {
  tryCatch(
    sv_model$log_density(theta, propto = FALSE),
    error = function(e) -1e12
  )
}

cor_parameter <- 0.95
h <- 1 
delta <- 0.0002
max_c <- 20
n_warmup_iterations <- 10000
n_sampling_iterations <- 10000
n_chains <- 10
nu <- 0.95

initial_mean_vector <- readRDS("numerical_experiments/sv_model/general_scripts/initial_mean_vector.RDS")

relevant_indices = c(1, 2, 3, d)

init_cluster <- parallel::makeCluster(10)
doParallel::registerDoParallel(init_cluster)

final_run <- foreach::foreach(i = 1:n_chains) %dopar% {
  sv_model <- bridgestan::StanModel$new(
    "numerical_experiments/sv_model/general_scripts/SVmodel.stan",
    "numerical_experiments/sv_model/general_scripts/SVmodel_data.json",
    seed = 123
  )
  
  d <- sv_model$param_num()
  
  grad_log_target_fun <- function(theta) {
    # sv_model$log_density_gradient(theta)$gradient
    grad_eval <- tryCatch(
      sv_model$log_density_gradient(theta)$gradient,
      error = function(e) rep(1e12, d)
    )
    
    if (is.na(grad_eval[1])) {
      grad_eval <- rep(1e12, d)
    }
    
    grad_eval
  }
  
  hamiltonian_func <- function(theta, rho) {
    # -sv_model$log_density(theta, propto = FALSE) + 0.5 * sum(rho ^ 2)
    -tryCatch(
      sv_model$log_density(theta, propto = FALSE),
      error = function(e) -1e12
    ) + 
      0.5 * sum(rho ^ 2)
  } 
  
  log_target_fun <- function(theta) {
    tryCatch(
      sv_model$log_density(theta, propto = FALSE),
      error = function(e) -1e12
    )
  }
  
  set.seed(i)
  # theta <- rnorm(100)
  # theta <- rnorm(d) * 0.1
  theta <- rnorm(d, mean = initial_mean_vector, sd = 0.1)
  sink(paste0("numerical_experiments/sv_model/kinetic_langevin_ghmc/hamiltonian_error/log/log_nr_", i, ".txt"))
  print("Warmup")
  single_warmup_run <- adaptive_step_size_one_step_GHMC_sampling(
    micro_fun = micro_fun_hamiltonian_error,
    n_samples = n_warmup_iterations, 
    h = h, 
    nu = nu,
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = theta,
    norm = "L-infinity",
    relevant_indices = relevant_indices
  ) 
  print("Sampling")
  single_sampling_run <- adaptive_step_size_one_step_GHMC_sampling(
    micro_fun = micro_fun_hamiltonian_error,
    n_samples = n_sampling_iterations, 
    h = h, 
    nu = nu,
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    theta = single_warmup_run$final_theta,
    norm = "L-infinity",
    relevant_indices = relevant_indices
  ) 
  sink()
  single_sampling_run
}

saveRDS(final_run, "numerical_experiments/sv_model/kinetic_langevin_ghmc/hamiltonian_error/sv_model_test.RDS")

parallel::stopCluster(init_cluster)

final_run <- readRDS(
  "numerical_experiments/sv_model/kinetic_langevin_ghmc/hamiltonian_error/sv_model_test.RDS"
)

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

final_result$samples_matrix[, 1] <- (-1 + 2 * exp(final_result$samples_matrix[, 1])) /
  (1 + exp(final_result$samples_matrix[, 1])) # transform back to a parameter between -1 and 1
final_result$samples_matrix[, 2] <- exp(final_result$samples_matrix[, 2]) # transform back to sigma

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

final_samples_3d_array[, , 1] <- (-1 + 2 * exp(final_samples_3d_array[, , 1])) /
  (1 + exp(final_samples_3d_array[, , 1])) # transform back to a parameter between -1 and 1
final_samples_3d_array[, , 2] <- exp(final_samples_3d_array[, , 2]) # transform back to sigma

relevant_indices <- c(1, 2, 3, dim(final_samples_3d_array)[3])
rstan_monitor_summary <- rstan::monitor(final_samples_3d_array[, , relevant_indices], warmup = 0) # only consider the two parameters and the first and last latent 
rstan_monitor_summary$n_eff
rstan_monitor_summary$n_eff * 1000000 / sum(final_result$n_evals_ode)

sink("numerical_experiments/sv_model/kinetic_langevin_ghmc/hamiltonian_error/sv_model_test.txt")
print("reversibility")
print(table(final_result$reversibility_indicator))
print("metropolis acceptance")
print(table(final_result$accept_indicator))
print("ratio c")
print(final_c_df$Freq[final_c_df$final_c_vector == 0] / 
        sum(final_c_df$Freq[final_c_df$final_c_vector != -1]))
print("n evals ode")
print(sum(final_result$n_evals_ode))
print("overall")
print(round(colMeans(final_result$samples_matrix[, relevant_indices]), 4))
print(round(apply(final_result$samples_matrix[, relevant_indices], 2, sd), 4))
print((rstan_monitor_summary$n_eff))
print((rstan_monitor_summary$n_eff) * 1e6 / sum(final_result$n_evals_ode)) 
sink()
