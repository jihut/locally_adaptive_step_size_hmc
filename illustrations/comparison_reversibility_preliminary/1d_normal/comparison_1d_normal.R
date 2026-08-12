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

theta_grid <- seq(from = -10, to = 10, by = 0.05)
rho_grid <- seq(from = -5, to = 5, by = 0.05)
theta_rho_matrix <- expand.grid(theta = theta_grid, rho = rho_grid)
n_iterations <- nrow(theta_rho_matrix)
L <- 1
max_c <- 20
h <- pi / 2

# Hamiltonian criterion

hamiltonian_reversibility_vec <- numeric(n_iterations)
hamiltonian_c_vec <- numeric(n_iterations)
hamiltonian_accept_vec <- numeric(n_iterations)
hamiltonian_delta <- 0.2
# hamiltonian_delta <- 0.005
for (i in 1:n_iterations) {
  if (i %% 10000 == 0) print(i)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_hamiltonian_error,
    h = h, 
    L = L, 
    delta = hamiltonian_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta_rho_matrix$theta[i], 
    rho = theta_rho_matrix$rho[i], 
    grad_log_target_initial = grad_log_target_fun(theta_rho_matrix$theta[i]),
    return_c_for_L_equal_to_1 = TRUE
  )
  hamiltonian_reversibility_vec[i] <- single_run$weight_non_zero_indicator
  hamiltonian_accept_vec[i] <- single_run$accept_indicator
  hamiltonian_c_vec[i] <- single_run$c_vector[1]
}

hamiltonian_df_deterministic <- data.frame(
  theta_rho_matrix,
  c = as.factor(hamiltonian_c_vec),
  reversibility = as.factor(hamiltonian_reversibility_vec)
)

table(hamiltonian_reversibility_vec)
table(hamiltonian_accept_vec)
range(hamiltonian_c_vec)

hamiltonian_plot_deterministic_nr1 <- ggplot(hamiltonian_df_deterministic, aes(x = theta, y = rho, color = c)) + 
  # geom_point(alpha = 0.1) +
  geom_raster(aes(x = theta, y = rho, fill = c)) +
  # labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", hamiltonian_delta)) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = pi / 2", ", delta = ", hamiltonian_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) +
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$")
  ) + 
  scale_fill_manual(
    values = c(
      "0" = "steelblue4", 
      "1" = "steelblue3",
      "2" = "steelblue2",
      "3" = "steelblue1",
      "4" = "steelblue"
    )
  )

hamiltonian_plot_deterministic_nr1
hamiltonian_plot_deterministic_nr2 <- ggplot(hamiltonian_df_deterministic, aes(x = theta, y = rho, color = as.factor(reversibility))) + 
  # geom_point(alpha = 0.1) + 
  geom_raster(aes(x = theta, y = rho, fill = reversibility)) +
  # labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", hamiltonian_delta)) + 
  labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = pi / 2", ", delta = ", hamiltonian_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) + 
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$"),
    color = "No reversibility issue",
    fill = "Reversible"
  ) + 
  scale_fill_discrete(
    labels = c("No", "Yes")
  )
hamiltonian_plot_deterministic_nr2

# ggplot(hamiltonian_df_deterministic, aes(x = theta, y = rho, z = c)) + 
#   geom_tile(
#     data = subset(hamiltonian_df_deterministic, reversibility == 0),
#     aes(fill = "red", alpha = 0.1)
#   ) + 
#   geom_contour(aes(color = after_stat(level))) + 
#   scale_color_viridis_c(name = "c") + 
#   xlim(c(-4, 4)) + 
#   ylim(c(-4, 4))

## Hamiltonian reversibility example

hamiltonian_c_matrix <- matrix(hamiltonian_c_vec, nrow = length(theta_grid), ncol = length(rho_grid))
hamiltonian_reversibility_matrix <- matrix(hamiltonian_reversibility_vec, nrow = length(theta_grid), ncol = length(rho_grid))
hamiltonian_reversibility_matrix <- ifelse(hamiltonian_reversibility_matrix == 0, 1, NA)
hamiltonian_c_levels <- sort(unique(hamiltonian_c_vec))
hamiltonian_c_colours   <- viridis::viridis(length(hamiltonian_c_levels))
contour(
  theta_grid,  
  rho_grid, 
  hamiltonian_c_matrix, 
  levels = 1 + 0:5,
  col = hamiltonian_c_colours,
  xlim = c(-10, 10),
  ylim = c(-5, 5),
  lwd = 2, 
  lty = 2,
  xlab = latex2exp::TeX("$\\theta$"),
  ylab = latex2exp::TeX("$\\rho$"),
  main = "Hamiltonian criterion"
)
image(
  theta_grid, rho_grid, hamiltonian_reversibility_matrix,
  col = adjustcolor("red", alpha.f = 0.2),
  add = TRUE
)
relevant_theta <- 1.45
relevant_rho <- -2
points(relevant_theta, relevant_rho, col = "blue", pch = 19)
relevant_forward_run <- micro_fun_hamiltonian_error(
  h = h, 
  delta = hamiltonian_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_theta, 
  rho = relevant_rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_theta)
)
relevant_forward_run$c
# lines(c(relevant_theta, relevant_forward_run$theta), c(relevant_rho, relevant_forward_run$rho), col = "red")
arrows(relevant_theta, relevant_rho, relevant_forward_run$theta, relevant_forward_run$rho, col = "darkblue")
# lines(c(relevant_forward_run$theta, relevant_forward_run$theta), c(relevant_forward_run$rho, -relevant_forward_run$rho), col = "yellow")
arrows(relevant_forward_run$theta, relevant_forward_run$rho, relevant_forward_run$theta, -relevant_forward_run$rho, col = "red")

relevant_backward_run <- micro_fun_hamiltonian_error(
  h = h, 
  delta = hamiltonian_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_forward_run$theta, 
  rho = -relevant_forward_run$rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_forward_run$theta)
)
# lines(c(relevant_forward_run$theta, relevant_backward_run$theta), c(-relevant_forward_run$rho, relevant_backward_run$rho), col = "red")
arrows(relevant_forward_run$theta, -relevant_forward_run$rho, relevant_backward_run$theta, relevant_backward_run$rho, col = "darkblue")
# lines(c(relevant_backward_run$theta, relevant_backward_run$theta), c(relevant_backward_run$rho, -relevant_backward_run$rho), col = "yellow")
arrows(relevant_backward_run$theta, relevant_backward_run$rho, relevant_backward_run$theta, -relevant_backward_run$rho, col = "red")
points(relevant_backward_run$theta, -relevant_backward_run$rho, col = "orange", pch = 19)

# Hamiltonian criterion - stochastic by sampling exactly from the joint target

hamiltonian_reversibility_vec_random <- numeric(n_iterations)
hamiltonian_c_vec_random <- numeric(n_iterations)
hamiltonian_accept_vec_random <- numeric(n_iterations)
hamiltonian_delta <- 0.2
# hamiltonian_delta <- 0.005
n_iterations_random <- 100000
set.seed(1)
for (i in 1:n_iterations_random) {
  if (i %% 10000 == 0) print(i)
  theta <- rnorm(1)
  rho <- rnorm(1)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_hamiltonian_error,
    h = h, 
    L = L, 
    delta = hamiltonian_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta, 
    rho = rho, 
    grad_log_target_initial = grad_log_target_fun(theta),
    return_c_for_L_equal_to_1 = TRUE
  )
  hamiltonian_reversibility_vec_random[i] <- single_run$weight_non_zero_indicator
  hamiltonian_accept_vec_random[i] <- single_run$accept_indicator
  hamiltonian_c_vec_random[i] <- single_run$c_vector[1]
}
table(hamiltonian_reversibility_vec_random)

# Original Hamiltonian criterion

og_hamiltonian_reversibility_vec <- numeric(n_iterations)
og_hamiltonian_c_vec <- numeric(n_iterations)
og_hamiltonian_accept_vec <- numeric(n_iterations)
og_hamiltonian_delta <- 0.2
for (i in 1:n_iterations) {
  if (i %% 10000 == 0) print(i)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_og_hamiltonian_error,
    h = h, 
    L = L, 
    delta = og_hamiltonian_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta_rho_matrix$theta[i], 
    rho = theta_rho_matrix$rho[i], 
    grad_log_target_initial = grad_log_target_fun(theta_rho_matrix$theta[i]),
    return_c_for_L_equal_to_1 = TRUE
  )
  og_hamiltonian_reversibility_vec[i] <- single_run$weight_non_zero_indicator
  og_hamiltonian_accept_vec[i] <- single_run$accept_indicator
  og_hamiltonian_c_vec[i] <- single_run$c_vector[1]
}

og_hamiltonian_df_deterministic <- data.frame(
  theta_rho_matrix,
  c = as.factor(og_hamiltonian_c_vec),
  reversibility = as.factor(og_hamiltonian_reversibility_vec)
)

table(og_hamiltonian_reversibility_vec)
table(og_hamiltonian_accept_vec)
range(og_hamiltonian_c_vec)

og_hamiltonian_plot_deterministic_nr1 <- ggplot(og_hamiltonian_df_deterministic, aes(x = theta, y = rho, color = c)) + 
  # geom_point(alpha = 0.1) +
  geom_raster(aes(x = theta, y = rho, fill = c)) +
  # labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", hamiltonian_delta)) + 
  labs(title = paste0("N(0, 1) - Original Hamiltonian criterion - L = ", L, ", h = pi / 2", ", delta = ", og_hamiltonian_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) +
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$")
  ) + 
  scale_fill_manual(
    values = c(
      "0" = "steelblue4", 
      "1" = "steelblue3",
      "2" = "steelblue2",
      "3" = "steelblue1",
      "4" = "steelblue"
    )
  )

og_hamiltonian_plot_deterministic_nr1
og_hamiltonian_plot_deterministic_nr2 <- ggplot(og_hamiltonian_df_deterministic, aes(x = theta, y = rho, color = as.factor(reversibility))) + 
  # geom_point(alpha = 0.1) + 
  geom_raster(aes(x = theta, y = rho, fill = reversibility)) +
  # labs(title = paste0("N(0, 1) - Hamiltonian criterion - L = ", L, ", h = ", h, ", delta = ", hamiltonian_delta)) + 
  labs(title = paste0("N(0, 1) - Original Hamiltonian criterion - L = ", L, ", h = pi / 2", ", delta = ", og_hamiltonian_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) + 
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$"),
    color = "No reversibility issue",
    fill = "Reversible"
  ) + 
  scale_fill_discrete(
    labels = c("No", "Yes")
  )
og_hamiltonian_plot_deterministic_nr2

# ggplot(hamiltonian_df_deterministic, aes(x = theta, y = rho, z = c)) + 
#   geom_tile(
#     data = subset(hamiltonian_df_deterministic, reversibility == 0),
#     aes(fill = "red", alpha = 0.1)
#   ) + 
#   geom_contour(aes(color = after_stat(level))) + 
#   scale_color_viridis_c(name = "c") + 
#   xlim(c(-4, 4)) + 
#   ylim(c(-4, 4))

## Original Hamiltonian reversibility example

og_hamiltonian_c_matrix <- matrix(og_hamiltonian_c_vec, nrow = length(theta_grid), ncol = length(rho_grid))
og_hamiltonian_reversibility_matrix <- matrix(og_hamiltonian_reversibility_vec, nrow = length(theta_grid), ncol = length(rho_grid))
og_hamiltonian_reversibility_matrix <- ifelse(og_hamiltonian_reversibility_matrix == 0, 1, NA)
og_hamiltonian_c_levels <- sort(unique(og_hamiltonian_c_vec))
og_hamiltonian_c_colours   <- viridis::viridis(length(og_hamiltonian_c_levels))
contour(
  theta_grid,  
  rho_grid, 
  og_hamiltonian_c_matrix, 
  levels = 1 + 0:5,
  col = og_hamiltonian_c_colours,
  xlim = c(-10, 10),
  ylim = c(-5, 5),
  lwd = 2, 
  lty = 2,
  xlab = latex2exp::TeX("$\\theta$"),
  ylab = latex2exp::TeX("$\\rho$"),
  main = "Original Hamiltonian criterion"
)
image(
  theta_grid, rho_grid, og_hamiltonian_reversibility_matrix,
  col = adjustcolor("red", alpha.f = 0.2),
  add = TRUE
)
relevant_theta <- 0.5
relevant_rho <- -0.75
points(relevant_theta, relevant_rho, col = "blue", pch = 19)
relevant_forward_run <- micro_fun_og_hamiltonian_error(
  h = h, 
  delta = og_hamiltonian_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_theta, 
  rho = relevant_rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_theta)
)
relevant_forward_run$c
# lines(c(relevant_theta, relevant_forward_run$theta), c(relevant_rho, relevant_forward_run$rho), col = "red")
arrows(relevant_theta, relevant_rho, relevant_forward_run$theta, relevant_forward_run$rho, col = "darkblue")
# lines(c(relevant_forward_run$theta, relevant_forward_run$theta), c(relevant_forward_run$rho, -relevant_forward_run$rho), col = "yellow")
arrows(relevant_forward_run$theta, relevant_forward_run$rho, relevant_forward_run$theta, -relevant_forward_run$rho, col = "red")

relevant_backward_run <- micro_fun_og_hamiltonian_error(
  h = h, 
  delta = og_hamiltonian_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_forward_run$theta, 
  rho = -relevant_forward_run$rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_forward_run$theta)
)
# lines(c(relevant_forward_run$theta, relevant_backward_run$theta), c(-relevant_forward_run$rho, relevant_backward_run$rho), col = "red")
arrows(relevant_forward_run$theta, -relevant_forward_run$rho, relevant_backward_run$theta, relevant_backward_run$rho, col = "darkblue")
# lines(c(relevant_backward_run$theta, relevant_backward_run$theta), c(relevant_backward_run$rho, -relevant_backward_run$rho), col = "yellow")
arrows(relevant_backward_run$theta, relevant_backward_run$rho, relevant_backward_run$theta, -relevant_backward_run$rho, col = "red")
points(relevant_backward_run$theta, -relevant_backward_run$rho, col = "orange", pch = 19)

# RKNF criterion

rknf_reversibility_vec <- numeric(n_iterations)
rknf_c_vec <- numeric(n_iterations)
rknf_accept_vec <- numeric(n_iterations)
rknf_delta <- 1
# rknf_delta <- 0.1
for (i in 1:n_iterations) {
  if (i %% 10000 == 0) print(i)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_rknf_error,
    h = h, 
    L = L, 
    delta = rknf_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta_rho_matrix$theta[i], 
    rho = theta_rho_matrix$rho[i], 
    grad_log_target_initial = grad_log_target_fun(theta_rho_matrix$theta[i]),
    return_c_for_L_equal_to_1 = TRUE
  )
  rknf_reversibility_vec[i] <- single_run$weight_non_zero_indicator
  rknf_accept_vec[i] <- single_run$accept_indicator
  rknf_c_vec[i] <- single_run$c_vector[1]
}

rknf_df_deterministic <- data.frame(
  theta_rho_matrix,
  c = as.factor(rknf_c_vec),
  reversibility = as.factor(rknf_reversibility_vec)
)

table(rknf_reversibility_vec)
table(rknf_accept_vec)
range(rknf_c_vec)

rknf_plot_deterministic_nr1 <- ggplot(rknf_df_deterministic, aes(x = theta, y = rho, color = c)) + 
  # geom_point(alpha = 0.1) + 
  geom_raster(aes(x = theta, y = rho, fill = c)) + 
  # labs(title = paste0("N(0, 1) - RKNF criterion - L = ", L, ", h = ", h, ", delta = ", rknf_delta)) + 
  labs(title = paste0("N(0, 1) - RKNF criterion - L = ", L, ", h = pi / 2", ", delta = ", rknf_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) +
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$")
  ) + 
  scale_fill_manual(
    values = c(
      "0" = "steelblue4", 
      "1" = "steelblue3"
    )
  )
rknf_plot_deterministic_nr1
rknf_plot_deterministic_nr2 <- ggplot(rknf_df_deterministic, aes(x = theta, y = rho, color = as.factor(reversibility))) + 
  # geom_point(alpha = 0.1) + 
  geom_raster(aes(x = theta, y = rho, fill = reversibility)) + 
  # labs(title = paste0("N(0, 1) - RKNF criterion - L = ", L, ", h = ", h, ", delta = ", rknf_delta)) + 
  labs(title = paste0("N(0, 1) - RKNF criterion - L = ", L, ", h = pi / 2", ", delta = ", rknf_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) + 
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$"),
    color = "No reversibility issue",
    fill = "Reversible"
  ) + 
  scale_fill_discrete(
    labels = c("No", "Yes")
  )
rknf_plot_deterministic_nr2

## RKNF reversibility example

rknf_c_matrix <- matrix(rknf_c_vec, nrow = length(theta_grid), ncol = length(rho_grid))
rknf_reversibility_matrix <- matrix(rknf_reversibility_vec, nrow = length(theta_grid), ncol = length(rho_grid))
rknf_reversibility_matrix <- ifelse(rknf_reversibility_matrix == 0, 1, NA)
rknf_c_levels <- sort(unique(rknf_c_vec))
rknf_c_colours   <- viridis::viridis(length(rknf_c_levels))
contour(
  theta_grid,  
  rho_grid, 
  rknf_c_matrix, 
  levels = 1 + 0:5,
  col = rknf_c_colours,
  xlim = c(-10, 10),
  ylim = c(-5, 5),
  lwd = 2, 
  lty = 2,
  xlab = latex2exp::TeX("$\\theta$"),
  ylab = latex2exp::TeX("$\\rho$"),
  main = "RKNF criterion"
)
image(
  theta_grid, rho_grid, rknf_reversibility_matrix,
  col = adjustcolor("red", alpha.f = 0.2),
  add = TRUE
)
relevant_theta <- -1
relevant_rho <- -2.7
points(relevant_theta, relevant_rho, col = "blue", pch = 19)
relevant_forward_run <- micro_fun_rknf_error(
  h = h, 
  delta = rknf_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_theta, 
  rho = relevant_rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_theta)
)
relevant_forward_run$c
# lines(c(relevant_theta, relevant_forward_run$theta), c(relevant_rho, relevant_forward_run$rho), col = "red")
arrows(relevant_theta, relevant_rho, relevant_forward_run$theta, relevant_forward_run$rho, col = "darkblue")
# lines(c(relevant_forward_run$theta, relevant_forward_run$theta), c(relevant_forward_run$rho, -relevant_forward_run$rho), col = "yellow")
arrows(relevant_forward_run$theta, relevant_forward_run$rho, relevant_forward_run$theta, -relevant_forward_run$rho, col = "red")

relevant_backward_run <- micro_fun_rknf_error(
  h = h, 
  delta = rknf_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_forward_run$theta, 
  rho = -relevant_forward_run$rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_forward_run$theta)
)
# lines(c(relevant_forward_run$theta, relevant_backward_run$theta), c(-relevant_forward_run$rho, relevant_backward_run$rho), col = "red")
arrows(relevant_forward_run$theta, -relevant_forward_run$rho, relevant_backward_run$theta, relevant_backward_run$rho, col = "darkblue")
# lines(c(relevant_backward_run$theta, relevant_backward_run$theta), c(relevant_backward_run$rho, -relevant_backward_run$rho), col = "yellow")
arrows(relevant_backward_run$theta, relevant_backward_run$rho, relevant_backward_run$theta, -relevant_backward_run$rho, col = "red")
points(relevant_backward_run$theta, -relevant_backward_run$rho, col = "orange", pch = 19)

# RKNF criterion - stochastic by sampling exactly from the joint target

rknf_reversibility_vec_random <- numeric(n_iterations)
rknf_c_vec_random <- numeric(n_iterations)
rknf_accept_vec_random <- numeric(n_iterations)
n_iterations_random <- 100000
set.seed(1)
for (i in 1:n_iterations_random) {
  if (i %% 10000 == 0) print(i)
  theta <- rnorm(1)
  rho <- rnorm(1)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_rknf_error,
    h = h, 
    L = L, 
    delta = rknf_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta, 
    rho = rho, 
    grad_log_target_initial = grad_log_target_fun(theta),
    return_c_for_L_equal_to_1 = TRUE
  )
  rknf_reversibility_vec_random[i] <- single_run$weight_non_zero_indicator
  rknf_accept_vec_random[i] <- single_run$accept_indicator
  rknf_c_vec_random[i] <- single_run$c_vector[1]
}
table(rknf_reversibility_vec_random)

# MRF criterion

mrf_reversibility_vec <- numeric(n_iterations)
mrf_c_vec <- numeric(n_iterations)
mrf_accept_vec <- numeric(n_iterations)
mrf_delta <- 1
for (i in 1:n_iterations) {
  if (i %% 10000 == 0) print(i)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_mrf_error,
    h = h, 
    L = L, 
    delta = mrf_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta_rho_matrix$theta[i], 
    rho = theta_rho_matrix$rho[i], 
    grad_log_target_initial = grad_log_target_fun(theta_rho_matrix$theta[i]),
    return_c_for_L_equal_to_1 = TRUE
  )
  mrf_reversibility_vec[i] <- single_run$weight_non_zero_indicator
  mrf_accept_vec[i] <- single_run$accept_indicator
  mrf_c_vec[i] <- single_run$c_vector[1]
}

mrf_df_deterministic <- data.frame(
  theta_rho_matrix,
  c = as.factor(mrf_c_vec + 1),
  reversibility = as.factor(mrf_reversibility_vec)
)

table(mrf_reversibility_vec)
table(mrf_accept_vec)
range(mrf_c_vec)

mrf_plot_deterministic_nr1 <- ggplot(mrf_df_deterministic, aes(x = theta, y = rho, color = c)) + 
  # geom_point(alpha = 0.1) + 
  geom_raster(aes(x = theta, y = rho, fill = c)) + 
  # labs(title = paste0("N(0, 1) - MRF criterion - L = ", L, ", h = ", h, ", delta = ", mrf_delta)) + 
  labs(title = paste0("N(0, 1) - MRF criterion - L = ", L, ", h = pi / 2", ", delta = ", mrf_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) +
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$")
  ) + 
  scale_fill_manual(
    values = c(
      "1" = "steelblue4", 
      "2" = "steelblue3"
    )
  )
mrf_plot_deterministic_nr1
mrf_plot_deterministic_nr2 <- ggplot(mrf_df_deterministic, aes(x = theta, y = rho, color = as.factor(reversibility))) + 
  # geom_point(alpha = 0.1) + 
  geom_raster(aes(x = theta, y = rho, fill = reversibility)) + 
  # labs(title = paste0("N(0, 1) - MRF criterion - L = ", L, ", h = ", h, ", delta = ", mrf_delta)) + 
  labs(title = paste0("N(0, 1) - MRF criterion - L = ", L, ", h = pi / 2", ", delta = ", mrf_delta)) + 
  xlim(c(-10, 10)) + 
  ylim(c(-5, 5)) + 
  labs(
    x = latex2exp::TeX("$\\theta$"),
    y = latex2exp::TeX("$\\rho$"),
    color = "No reversibility issue",
    fill = "Reversible"
  ) + 
  scale_fill_discrete(
    labels = c("No", "Yes")
  )
mrf_plot_deterministic_nr2

## MRF reversibility example

mrf_c_matrix <- matrix(mrf_c_vec, nrow = length(theta_grid), ncol = length(rho_grid))
mrf_reversibility_matrix <- matrix(mrf_reversibility_vec, nrow = length(theta_grid), ncol = length(rho_grid))
mrf_reversibility_matrix <- ifelse(mrf_reversibility_matrix == 0, 1, NA)
mrf_c_levels <- sort(unique(mrf_c_vec))
mrf_c_colours   <- viridis::viridis(length(mrf_c_levels))
contour(
  theta_grid,  
  rho_grid, 
  mrf_c_matrix, 
  levels = 1 + 0:5,
  col = mrf_c_colours,
  xlim = c(-10, 10),
  ylim = c(-5, 5),
  lwd = 2, 
  lty = 2,
  xlab = latex2exp::TeX("$\\theta$"),
  ylab = latex2exp::TeX("$\\rho$"),
  main = "MRF criterion"
)
image(
  theta_grid, rho_grid, mrf_reversibility_matrix,
  col = adjustcolor("red", alpha.f = 0.2),
  add = TRUE
)
relevant_theta <- -1.00
relevant_rho <- 1.75
points(relevant_theta, relevant_rho, col = "blue", pch = 19)
relevant_forward_run <- micro_fun_mrf_error(
  h = h, 
  delta = mrf_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_theta, 
  rho = relevant_rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_theta)
)
relevant_forward_run$c
# lines(c(relevant_theta, relevant_forward_run$theta), c(relevant_rho, relevant_forward_run$rho), col = "red")
arrows(relevant_theta, relevant_rho, relevant_forward_run$theta, relevant_forward_run$rho, col = "darkblue")
# lines(c(relevant_forward_run$theta, relevant_forward_run$theta), c(relevant_forward_run$rho, -relevant_forward_run$rho), col = "yellow")
arrows(relevant_forward_run$theta, relevant_forward_run$rho, relevant_forward_run$theta, -relevant_forward_run$rho, col = "red")

relevant_backward_run <- micro_fun_mrf_error(
  h = h, 
  delta = mrf_delta, 
  max_c = max_c, 
  log_target = log_target_fun, 
  grad_log_target = grad_log_target_fun, 
  hamiltonian = hamiltonian_func,
  theta = relevant_forward_run$theta, 
  rho = -relevant_forward_run$rho, 
  grad_log_target_initial = grad_log_target_fun(relevant_forward_run$theta)
)
# lines(c(relevant_forward_run$theta, relevant_backward_run$theta), c(-relevant_forward_run$rho, relevant_backward_run$rho), col = "red")
arrows(relevant_forward_run$theta, -relevant_forward_run$rho, relevant_backward_run$theta, relevant_backward_run$rho, col = "darkblue")
# lines(c(relevant_backward_run$theta, relevant_backward_run$theta), c(relevant_backward_run$rho, -relevant_backward_run$rho), col = "yellow")
arrows(relevant_backward_run$theta, relevant_backward_run$rho, relevant_backward_run$theta, -relevant_backward_run$rho, col = "red")
points(relevant_backward_run$theta, -relevant_backward_run$rho, col = "orange", pch = 19)

# MRF criterion - stochastic by sampling exactly from the joint target

mrf_reversibility_vec_random <- numeric(n_iterations)
mrf_c_vec_random <- numeric(n_iterations)
mrf_accept_vec_random <- numeric(n_iterations)
n_iterations_random <- 100000
set.seed(1)
for (i in 1:n_iterations_random) {
  if (i %% 10000 == 0) print(i)
  theta <- rnorm(1)
  rho <- rnorm(1)
  single_run <- single_iteration_standard_HMC(
    micro_fun = micro_fun_mrf_error,
    h = h, 
    L = L, 
    delta = mrf_delta, 
    max_c = max_c, 
    log_target = log_target_fun, 
    grad_log_target = grad_log_target_fun, 
    hamiltonian = hamiltonian_func,
    theta = theta, 
    rho = rho, 
    grad_log_target_initial = grad_log_target_fun(theta),
    return_c_for_L_equal_to_1 = TRUE
  )
  mrf_reversibility_vec_random[i] <- single_run$weight_non_zero_indicator
  mrf_accept_vec_random[i] <- single_run$accept_indicator
  mrf_c_vec_random[i] <- single_run$c_vector[1]
}
table(mrf_reversibility_vec_random)
###

png("illustrations/comparison_reversibility_preliminary/1d_normal/comparison_1d_normal.png", width = 3500, height = 3500/1.618, res = 300)
gridExtra::grid.arrange(
  hamiltonian_plot_deterministic_nr1, hamiltonian_plot_deterministic_nr2,
  rknf_plot_deterministic_nr1, rknf_plot_deterministic_nr2,
  mrf_plot_deterministic_nr1, mrf_plot_deterministic_nr2,
  nrow = 3, ncol = 2
)
dev.off()
