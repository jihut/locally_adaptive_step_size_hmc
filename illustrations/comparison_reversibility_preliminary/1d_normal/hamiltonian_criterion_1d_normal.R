rm(list = ls())

library(ggplot2)
source("implementation_scripts/standard_hmc/standard_hmc_functions.R")

mu <- 0
inv_var <- 1

grad_log_target_fun <- function(theta) {
  as.numeric(-inv_var * (theta - mu))
}

log_target_fun <- function(theta) dnorm(theta, mean = mu, sd = 1 / sqrt(inv_var), log = T)

hamiltonian_func <- function(theta, rho) {
  -dnorm(theta, mean = mu, sd = 1 / sqrt(inv_var), log = T) + 0.5 * rho ^ 2
}

# Hamiltonian criterion

set.seed(1)
theta <- rnorm(1)
rho <- rnorm(1)
h <- pi / 2 
L <- 1
n_samples <- 100000
delta <- 0.2 # energy tolerance
# delta <- 0.005
max_c <- 20
grad_log_target_initial <- grad_log_target_fun(theta)

sample_run <- adaptive_step_size_standard_HMC(
  micro_fun = micro_fun_hamiltonian_error,
  n_samples = n_samples, 
  h = h, 
  L = L, 
  delta = delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  theta = theta,
  return_c_for_L_equal_to_1 = TRUE
)

table(sample_run$reversibility_indicator)
mean(sample_run$reversibility_indicator)
table(sample_run$accept_indicator)
mean(sample_run$accept_indicator)
sample_run_c_vector <- as.numeric(sample_run$c_matrix)
sample_run_c_table <- table(sample_run_c_vector)
sample_run_c_df <- data.frame(sample_run_c_table)
sample_run_c_df
sample_run_c_df$Freq[sample_run_c_df$sample_run_c_vector == 0] / 
  sum(sample_run_c_df$Freq[sample_run_c_df$sample_run_c_vector != -1])

sample_run_df_sim <- data.frame(initial_theta = sample_run$initial_theta_matrix[, 1], 
                     initial_rho = sample_run$initial_rho_matrix[, 1], 
                     num_forward_micro_steps = 2 ^ sample_run_c_vector,
                     reversibility = sample_run$reversibility_indicator,
                     accept = sample_run$accept_indicator) 
plot_sim_nr1 <- ggplot(sample_run_df_sim, aes(x = initial_theta, y = initial_rho, color = num_forward_micro_steps)) + 
  geom_point(alpha = 0.1) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5))
plot_sim_nr1
plot_sim_nr2 <- ggplot(sample_run_df_sim, aes(x = initial_theta, y = initial_rho, color = as.factor(reversibility))) + 
  geom_point(alpha = 0.1) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5))
plot_sim_nr2
plot_sim_nr3 <- ggplot(sample_run_df_sim, aes(x = initial_theta, y = initial_rho, color = as.factor(accept))) + 
  geom_point(alpha = 0.1) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5))
plot_sim_nr3

# Deterministic example

theta_grid <- seq(from = -10, to = 10, by = 0.05)
rho_grid <- seq(from = -5, to = 5, by = 0.05)
theta_rho_matrix <- expand.grid(theta = theta_grid, rho = rho_grid)

n_iterations <- nrow(theta_rho_matrix)

reversibility_vec <- numeric(n_iterations)
c_vec <- numeric(n_iterations)
accept_vec <- numeric(n_iterations)

for (i in 1:n_iterations) {
  if (i %% 10000 == 0) print(i)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_hamiltonian_error,
    h = h, 
    L = L, 
    delta = delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta_rho_matrix$theta[i], 
    rho = theta_rho_matrix$rho[i], 
    grad_log_target_initial = grad_log_target_fun(theta_rho_matrix$theta[i]),
    return_c_for_L_equal_to_1 = TRUE
  )
  reversibility_vec[i] <- single_run$weight_non_zero_indicator
  accept_vec[i] <- single_run$accept_indicator
  c_vec[i] <- single_run$c_vector[1]
}

df_deterministic <- data.frame(
  theta_rho_matrix,
  num_forward_micro_steps = 2 ^ c_vec,
  reversibility = reversibility_vec
)

table(reversibility_vec)
table(accept_vec)
range(c_vec)

plot_deterministic_nr1 <- ggplot(df_deterministic, aes(x = theta, y = rho, color = num_forward_micro_steps)) + 
  geom_point(alpha = 0.1) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5))
plot_deterministic_nr1
plot_deterministic_nr2 <- ggplot(df_deterministic, aes(x = theta, y = rho, color = as.factor(reversibility))) + 
  geom_point(alpha = 0.1) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5))
plot_deterministic_nr2
plot_deterministic_nr3 <- ggplot(df_deterministic, aes(x = theta, y = rho, color = as.factor(accept_indicator))) + 
  geom_point(alpha = 0.1) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5))
plot_deterministic_nr3

# Find the boundary for the case where the absolute difference in Hamiltonian is exactly equal to delta

# Number of forward steps = 1 (solutions found via Maple): 

boundary_nr1_fun_nr1 <- function(theta) {
  numerator <- 0.5 * h ^ 3 * theta - h * theta + 0.2 * sqrt(25 * h ^ 2 * theta ^ 2 - 40)
  denominator <- h ^ 2 
  numerator / denominator
}
boundary_nr1_fun_nr2 <- function(theta) {
  numerator <- 0.5 * h ^ 3 * theta - h * theta - 0.2 * sqrt(25 * h ^ 2 * theta ^ 2 - 40)
  denominator <- h ^ 2 
  numerator / denominator
}
boundary_nr1_fun_nr3 <- function(theta) {
  numerator <- 0.5 * h ^ 3 * theta - h * theta + 0.2 * sqrt(25 * h ^ 2 * theta ^ 2 + 40)
  denominator <- h ^ 2 
  numerator / denominator
}
boundary_nr1_fun_nr4 <- function(theta) {
  numerator <- 0.5 * h ^ 3 * theta - h * theta - 0.2 * sqrt(25 * h ^ 2 * theta ^ 2 + 40)
  denominator <- h ^ 2 
  numerator / denominator
}
plot(theta_rho_matrix[df_deterministic$num_forward_micro_steps == 1, ], cex = 0.1)
lines(theta_grid, boundary_nr1_fun_nr1(theta_grid), col = "red")
lines(theta_grid, boundary_nr1_fun_nr2(theta_grid), col = "red")
lines(theta_grid, boundary_nr1_fun_nr3(theta_grid), col = "red")
lines(theta_grid, boundary_nr1_fun_nr4(theta_grid), col = "red")

# Number of forward steps = 2

boundary_nr2_fun_nr1 <- function(theta) {
  numerator <- 0.025 * (10 * h ^ 5 * theta - 160 * h ^ 3 * theta + 320 * h * theta - 64 * sqrt(25 * h ^ 2 * theta ^ 2 - 160))
  denominator <- (h ^ 2 - 8) * h ^ 2 
  numerator / denominator
}
boundary_nr2_fun_nr2 <- function(theta) {
  numerator <- 0.025 * (10 * h ^ 5 * theta - 160 * h ^ 3 * theta + 320 * h * theta + 64 * sqrt(25 * h ^ 2 * theta ^ 2 - 160))
  denominator <- (h ^ 2 - 8) * h ^ 2 
  numerator / denominator
}
boundary_nr2_fun_nr3 <- function(theta) {
  numerator <- 0.025 * (10 * h ^ 5 * theta - 160 * h ^ 3 * theta + 320 * h * theta - 64 * sqrt(25 * h ^ 2 * theta ^ 2 + 160))
  denominator <- (h ^ 2 - 8) * h ^ 2 
  numerator / denominator
}
boundary_nr2_fun_nr4 <- function(theta) {
  numerator <- 0.025 * (10 * h ^ 5 * theta - 160 * h ^ 3 * theta + 320 * h * theta + 64 * sqrt(25 * h ^ 2 * theta ^ 2 + 160))
  denominator <- (h ^ 2 - 8) * h ^ 2 
  numerator / denominator
}
plot(theta_rho_matrix[df_deterministic$num_forward_micro_steps == 2, ], cex = 0.1)
lines(theta_grid, boundary_nr2_fun_nr1(theta_grid), col = "red")
lines(theta_grid, boundary_nr2_fun_nr2(theta_grid), col = "red")
lines(theta_grid, boundary_nr2_fun_nr3(theta_grid), col = "red")
lines(theta_grid, boundary_nr2_fun_nr4(theta_grid), col = "red")
lines(theta_grid, boundary_nr1_fun_nr1(theta_grid), col = "green")
lines(theta_grid, boundary_nr1_fun_nr2(theta_grid), col = "green")
lines(theta_grid, boundary_nr1_fun_nr3(theta_grid), col = "green")
lines(theta_grid, boundary_nr1_fun_nr4(theta_grid), col = "green")
