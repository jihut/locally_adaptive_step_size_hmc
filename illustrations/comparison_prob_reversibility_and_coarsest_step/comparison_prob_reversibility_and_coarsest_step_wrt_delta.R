rm(list = ls())
source("implementation_scripts/general_scripts/tuning_criterion_functions.R")
library(parallel)
library(doParallel)
library(ggplot2)
library(doRNG)
library(gridExtra)
library(grid)

mu <- c(0, 0)
sigma <- matrix(c(1, 0.95, 0.95, 1), nrow = 2, byrow = T)
inv_sigma <- solve(sigma)

grad_log_target_fun <- function(theta) {
  as.numeric(-inv_sigma %*% (theta - mu))
}

hamiltonian_func <- function(theta, rho) {
  # -mvtnorm::dmvnorm(theta, mean = mu, sigma = sigma, log = T) + 0.5 * sum(rho ^ 2)
  0.5 * ((theta - mu) %*% inv_sigma) %*% (theta - mu) + 0.5 * sum(rho ^ 2)
}

log_target_fun <- function(theta) -0.5 * ((theta - mu) %*% inv_sigma) %*% (theta - mu)

# Hamiltonian - fixed h and m, varying delta
h <- 1
m <- 5
n_iterations_per_value <- 100000

# delta_vec <- c(0.001, 0.005, seq(0.01, 0.09, 0.01), seq(0.1, 2.5, 0.1))
delta_vec <- exp(seq(from = log(1e-3), to = 1, length.out = 24))
init_cluster <- parallel::makeCluster(8)
doParallel::registerDoParallel(init_cluster)
hamiltonian_results <- foreach::foreach(k = 1:length(delta_vec), .combine = "rbind", .options.RNG = 123) %dorng% {
  reversibility_issue <- numeric(n_iterations_per_value)
  c_vec <- numeric(n_iterations_per_value)
  delta <- delta_vec[k]
  for (i in 1:n_iterations_per_value) {
    cor_parameter <- 0.95
    theta <- numeric(2)
    theta[1] <- rnorm(1)
    for (j in 2:2) {
      theta[j] <- cor_parameter * theta[j - 1] + rnorm(1, sd = sqrt(1 - cor_parameter ^ 2))
    }
    rho <- rnorm(2)
    
    forward_run <- micro_fun_hamiltonian_error(
      h = h,
      delta = delta,
      max_c = 20,
      log_target = log_target_fun, 
      grad_log_target = grad_log_target_fun,
      hamiltonian = hamiltonian_func,
      theta = theta,
      rho = rho,
      grad_log_target_initial = grad_log_target_fun(theta)
    )
    
    backward_run <- micro_fun_hamiltonian_error(
      h = h,
      delta = delta,
      max_c = 20,
      log_target = log_target_fun, 
      grad_log_target = grad_log_target_fun,
      hamiltonian = hamiltonian_func,
      theta = forward_run$theta,
      rho = -forward_run$rho,
      grad_log_target_initial = forward_run$grad_log_target
    )
    
    reversibility_issue[i] <- forward_run$c != backward_run$c
    c_vec[i] <- forward_run$c
  }
  
  prob_reversiblity_issue <- mean(reversibility_issue)
  prob_reversiblity_issue_conf_int <- binom.test(x = sum(reversibility_issue), n = n_iterations_per_value)
  prob_c_0 <- mean(c_vec == 0)
  prob_c_1 <- mean(c_vec == 1)
  prob_c_2 <- mean(c_vec == 2)
  prob_c_larger_2 <- mean(c_vec > 2)
  data.frame(
    prob_reversiblity_issue = prob_reversiblity_issue, 
    lower_prob_reversibility_issue = prob_reversiblity_issue_conf_int$conf.int[1],
    upper_prob_reversibility_issue = prob_reversiblity_issue_conf_int$conf.int[2],
    prob_c_0 = prob_c_0,
    prob_c_1 = prob_c_1,
    prob_c_2 = prob_c_2,
    prob_c_larger_2 = prob_c_larger_2
  )
}
hamiltonian_results 

# Runge-Kutta-Nyström flow

rknf_results <- foreach::foreach(k = 1:length(delta_vec), .combine = "rbind", .options.RNG = 123) %dorng% {
  reversibility_issue <- numeric(n_iterations_per_value)
  c_vec <- numeric(n_iterations_per_value)
  delta <- delta_vec[k]
  for (i in 1:n_iterations_per_value) {
    cor_parameter <- 0.95
    theta <- numeric(2)
    theta[1] <- rnorm(1)
    for (j in 2:2) {
      theta[j] <- cor_parameter * theta[j - 1] + rnorm(1, sd = sqrt(1 - cor_parameter ^ 2))
    }
    rho <- rnorm(2)
    
    forward_run <- micro_fun_rknf_error(
      h = h,
      delta = delta,
      max_c = 20,
      log_target = log_target_fun, 
      grad_log_target = grad_log_target_fun,
      hamiltonian = hamiltonian_func,
      theta = theta,
      rho = rho,
      grad_log_target_initial = grad_log_target_fun(theta)
    )
    
    backward_run <- micro_fun_rknf_error(
      h = h,
      delta = delta,
      max_c = 20,
      log_target = log_target_fun, 
      grad_log_target = grad_log_target_fun,
      hamiltonian = hamiltonian_func,
      theta = forward_run$theta,
      rho = -forward_run$rho,
      grad_log_target_initial = forward_run$grad_log_target
    )
    
    reversibility_issue[i] <- forward_run$c != backward_run$c
    c_vec[i] <- forward_run$c
  }
  
  prob_reversiblity_issue <- mean(reversibility_issue)
  prob_reversiblity_issue_conf_int <- binom.test(x = sum(reversibility_issue), n = n_iterations_per_value)
  prob_c_0 <- mean(c_vec == 0)
  prob_c_1 <- mean(c_vec == 1)
  prob_c_2 <- mean(c_vec == 2)
  prob_c_larger_2 <- mean(c_vec > 2)
  data.frame(
    prob_reversiblity_issue = prob_reversiblity_issue, 
    lower_prob_reversibility_issue = prob_reversiblity_issue_conf_int$conf.int[1],
    upper_prob_reversibility_issue = prob_reversiblity_issue_conf_int$conf.int[2],
    prob_c_0 = prob_c_0,
    prob_c_1 = prob_c_1,
    prob_c_2 = prob_c_2,
    prob_c_larger_2 = prob_c_larger_2
  )
}
rknf_results

# Multiresolution flow

mrf_results <- foreach::foreach(k = 1:length(delta_vec), .combine = "rbind", .options.RNG = 123) %dorng% {
  reversibility_issue <- numeric(n_iterations_per_value)
  c_vec <- numeric(n_iterations_per_value)
  delta <- delta_vec[k]
  for (i in 1:n_iterations_per_value) {
    cor_parameter <- 0.95
    theta <- numeric(2)
    theta[1] <- rnorm(1)
    for (j in 2:2) {
      theta[j] <- cor_parameter * theta[j - 1] + rnorm(1, sd = sqrt(1 - cor_parameter ^ 2))
    }
    rho <- rnorm(2)
    
    forward_run <- micro_fun_mrf_error(
      h = h,
      delta = delta,
      max_c = 20,
      log_target = log_target_fun, 
      grad_log_target = grad_log_target_fun,
      hamiltonian = hamiltonian_func,
      theta = theta,
      rho = rho,
      grad_log_target_initial = grad_log_target_fun(theta)
    )
    
    backward_run <- micro_fun_mrf_error(
      h = h,
      delta = delta,
      max_c = 20,
      log_target = log_target_fun, 
      grad_log_target = grad_log_target_fun,
      hamiltonian = hamiltonian_func,
      theta = forward_run$theta,
      rho = -forward_run$rho,
      grad_log_target_initial = forward_run$grad_log_target
    )
    
    reversibility_issue[i] <- forward_run$c != backward_run$c
    c_vec[i] <- forward_run$c
  }
  
  prob_reversiblity_issue <- mean(reversibility_issue)
  prob_reversiblity_issue_conf_int <- binom.test(x = sum(reversibility_issue), n = n_iterations_per_value)
  prob_c_0 <- mean(c_vec == 0)
  prob_c_1 <- mean(c_vec == 1)
  prob_c_2 <- mean(c_vec == 2)
  prob_c_larger_2 <- mean(c_vec > 2)
  data.frame(
    prob_reversiblity_issue = prob_reversiblity_issue, 
    lower_prob_reversibility_issue = prob_reversiblity_issue_conf_int$conf.int[1],
    upper_prob_reversibility_issue = prob_reversiblity_issue_conf_int$conf.int[2],
    prob_c_0 = prob_c_0,
    prob_c_1 = prob_c_1,
    prob_c_2 = prob_c_2,
    prob_c_larger_2 = prob_c_larger_2
  )
}
mrf_results

list_results <- list(
  hamiltonian_results = hamiltonian_results,
  rknf_results = rknf_results, 
  mrf_results = mrf_results
)
saveRDS(
  list_results, 
  "illustrations/comparison_prob_reversibility_and_coarsest_step/varying_delta_results.RDS"
)

parallel::stopCluster(init_cluster)

list_results <- readRDS("illustrations/comparison_prob_reversibility_and_coarsest_step/varying_delta_results.RDS")
hamiltonian_results <- list_results$hamiltonian_results
rknf_results <- list_results$rknf_results
mrf_results <- list_results$mrf_results

par(mfrow = c(3, 2))
plot(delta_vec, hamiltonian_results$prob_reversiblity_issue, type = "l", ylim = c(0, 0.18))
lines(delta_vec, hamiltonian_results$lower_prob_reversibility_issue, lty = 2, col = "red")
lines(delta_vec, hamiltonian_results$upper_prob_reversibility_issue, lty = 2, col = "red")
plot(delta_vec, hamiltonian_results$prob_c_0, type = "l", ylim = c(0, 1))
lines(delta_vec, hamiltonian_results$prob_c_1, col = "green")
lines(delta_vec, hamiltonian_results$prob_c_2, col = "blue")
lines(delta_vec, hamiltonian_results$prob_c_larger_2, col = "red")
plot(delta_vec, rknf_results$prob_reversiblity_issue, type = "l", ylim = c(0, 0.18))
lines(delta_vec, rknf_results$lower_prob_reversibility_issue, lty = 2, col = "red")
lines(delta_vec, rknf_results$upper_prob_reversibility_issue, lty = 2, col = "red")
plot(delta_vec, rknf_results$prob_c_0, type = "l", ylim = c(0, 1))
lines(delta_vec, rknf_results$prob_c_1, col = "green")
lines(delta_vec, rknf_results$prob_c_2, col = "blue")
lines(delta_vec, rknf_results$prob_c_larger_2, col = "red")
plot(delta_vec, mrf_results$prob_reversiblity_issue, type = "l", ylim = c(0, 0.18))
lines(delta_vec, mrf_results$lower_prob_reversibility_issue, lty = 2, col = "red")
lines(delta_vec, mrf_results$upper_prob_reversibility_issue, lty = 2, col = "red")
plot(delta_vec, mrf_results$prob_c_0, type = "l", ylim = c(0, 1))
lines(delta_vec, mrf_results$prob_c_1, col = "green")
lines(delta_vec, mrf_results$prob_c_2, col = "blue")
lines(delta_vec, mrf_results$prob_c_larger_2, col = "red")
par(mfrow = c(1, 1))

## GGPLOT style 

# Hamiltonian

hamiltonian_results_df <- cbind(delta_vec, hamiltonian_results)
plot_hamiltonian_prob_reversibility_issue <- ggplot(data = hamiltonian_results_df) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_reversiblity_issue
    ),
    color = "black"
  ) +
  geom_ribbon(
    aes(
      x = delta_vec,
      ymin = lower_prob_reversibility_issue,
      ymax = upper_prob_reversibility_issue
    ),
    alpha = 0.15,
    fill = "red",
    linetype = 0
  ) +
  ylim(
    c(0, 0.18)
  ) +
  labs(
    y = "Probability",
    x = "delta"
  ) + 
  scale_x_log10()
plot_hamiltonian_prob_reversibility_issue

plot_hamiltonian_prob_c <- ggplot(data = hamiltonian_results_df) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_0,
      color = "P(c = 0)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_1,
      color = "P(c = 1)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_2,
      color = "P(c = 2)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_larger_2,
      color = "P(c > 2)"
    )
  ) +
  labs(
    y = "Probability",
    x = "delta"
  ) +
  scale_color_manual(
    name = "c",
    values = c(
      "P(c = 0)" = "black",
      "P(c = 1)" = "green",
      "P(c = 2)" = "blue",
      "P(c > 2)" = "red"
    )
  ) + 
  scale_x_log10()
plot_hamiltonian_prob_c

# RKNF
rknf_results_df <- cbind(delta_vec, rknf_results)
plot_rknf_prob_reversibility_issue <- ggplot(data = rknf_results_df) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_reversiblity_issue
    ),
    color = "black"
  ) +
  geom_ribbon(
    aes(
      x = delta_vec,
      ymin = lower_prob_reversibility_issue,
      ymax = upper_prob_reversibility_issue
    ),
    alpha = 0.15,
    fill = "red",
    linetype = 0
  ) +
  ylim(
    c(0, 0.18)
  ) +
  labs(
    y = "Probability",
    x = "delta"
  ) + 
  scale_x_log10()
plot_rknf_prob_reversibility_issue

plot_rknf_prob_c <- ggplot(data = rknf_results_df) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_0,
      color = "P(c = 0)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_1,
      color = "P(c = 1)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_2,
      color = "P(c = 2)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_larger_2,
      color = "P(c > 2)"
    )
  ) +
  labs(
    y = "Probability",
    x = "delta"
  ) +
  scale_color_manual(
    name = "c",
    values = c(
      "P(c = 0)" = "black",
      "P(c = 1)" = "green",
      "P(c = 2)" = "blue",
      "P(c > 2)" = "red"
    )
  ) + 
  scale_x_log10()
plot_rknf_prob_c

# Multiresolution flow
mrf_results_df <- cbind(delta_vec, mrf_results)
plot_mrf_prob_reversibility_issue <- ggplot(data = mrf_results_df) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_reversiblity_issue
    ),
    color = "black"
  ) +
  geom_ribbon(
    aes(
      x = delta_vec,
      ymin = lower_prob_reversibility_issue,
      ymax = upper_prob_reversibility_issue
    ),
    alpha = 0.15,
    fill = "red",
    linetype = 0
  ) +
  ylim(
    c(0, 0.18)
  ) +
  labs(
    y = "Probability",
    x = "delta"
  ) + 
  scale_x_log10()
plot_mrf_prob_reversibility_issue

plot_mrf_prob_c <- ggplot(data = mrf_results_df) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_0,
      color = "P(c = 1)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_1,
      color = "P(c = 2)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_2,
      color = "P(c = 3)"
    )
  ) +
  geom_line(
    aes(
      x = delta_vec,
      y = prob_c_larger_2,
      color = "P(c > 3)"
    )
  ) +
  labs(
    y = "Probability",
    x = "delta"
  ) +
  scale_color_manual(
    name = "c",
    values = c(
      "P(c = 1)" = "green",
      "P(c = 2)" = "blue",
      "P(c = 3)" = "red",
      "P(c > 3)" = "orange"
    )
  ) + 
  scale_x_log10()
plot_mrf_prob_c

gridExtra::grid.arrange(
  plot_hamiltonian_prob_reversibility_issue, plot_hamiltonian_prob_c,
  plot_rknf_prob_reversibility_issue, plot_rknf_prob_c,
  plot_mrf_prob_reversibility_issue, plot_mrf_prob_c,
  nrow = 3,
  ncol = 2
)

# Alternative:

hamiltonian_row <- gridExtra::arrangeGrob(
  grid::textGrob("Hamiltonian criterion"),
  gridExtra::arrangeGrob(plot_hamiltonian_prob_reversibility_issue, plot_hamiltonian_prob_c, ncol = 2),
  ncol = 1,
  heights = c(0.15, 1)
)

rknf_row <- gridExtra::arrangeGrob(
  grid::textGrob("Runge-Kutta-Nyström flow criterion"),
  gridExtra::arrangeGrob(plot_rknf_prob_reversibility_issue, plot_rknf_prob_c, ncol = 2),
  ncol = 1,
  heights = c(0.15, 1)
)

mrf_row <- gridExtra::arrangeGrob(
  grid::textGrob("Multiresolution flow criterion"),
  gridExtra::arrangeGrob(plot_mrf_prob_reversibility_issue, plot_mrf_prob_c, ncol = 2),
  ncol = 1,
  heights = c(0.15, 1)
)

grid.arrange(hamiltonian_row, rknf_row, mrf_row, ncol = 1)

# # Hamiltonian
# hamiltonian_results_df <- cbind(delta_vec, hamiltonian_results)
# plot_hamiltonian_prob_reversibility_issue <- ggplot(data = hamiltonian_results_df) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_reversiblity_issue
#     ),
#     color = "black"
#   ) + 
#   geom_ribbon(
#     aes(
#       x = log(delta_vec), 
#       ymin = lower_prob_reversibility_issue,
#       ymax = upper_prob_reversibility_issue
#     ),
#     alpha = 0.15,
#     fill = "red",
#     linetype = 0
#   ) + 
#   ylim(
#     c(0, 0.18)
#   ) + 
#   labs(
#     y = "Probability",
#     x = "log(delta)"
#   )
# plot_hamiltonian_prob_reversibility_issue
# 
# plot_hamiltonian_prob_c <- ggplot(data = hamiltonian_results_df) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_0,
#       color = "P(c = 0)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_1,
#       color = "P(c = 1)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_2,
#       color = "P(c = 2)"
#     )
#   ) +
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_larger_2,
#       color = "P(c > 2)"
#     )
#   ) + 
#   labs(
#     y = "Probability",
#     x = "log(delta)"
#   ) + 
#   scale_color_manual(
#     name = "c",
#     values = c(
#       "P(c = 0)" = "black",
#       "P(c = 1)" = "green",
#       "P(c = 2)" = "blue",
#       "P(c > 2)" = "red"
#     )
#   )
# plot_hamiltonian_prob_c
# 
# # RKNF
# rknf_results_df <- cbind(delta_vec, rknf_results)
# plot_rknf_prob_reversibility_issue <- ggplot(data = rknf_results_df) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_reversiblity_issue
#     ),
#     color = "black"
#   ) + 
#   geom_ribbon(
#     aes(
#       x = log(delta_vec), 
#       ymin = lower_prob_reversibility_issue,
#       ymax = upper_prob_reversibility_issue
#     ),
#     alpha = 0.15,
#     fill = "red",
#     linetype = 0
#   ) + 
#   ylim(
#     c(0, 0.18)
#   ) + 
#   labs(
#     y = "Probability",
#     x = "log(delta)"
#   )
# plot_rknf_prob_reversibility_issue
# 
# plot_rknf_prob_c <- ggplot(data = rknf_results_df) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_0,
#       color = "P(c = 0)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_1,
#       color = "P(c = 1)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_2,
#       color = "P(c = 2)"
#     )
#   ) +
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_larger_2,
#       color = "P(c > 2)"
#     )
#   ) + 
#   labs(
#     y = "Probability",
#     x = "log(delta)"
#   ) + 
#   scale_color_manual(
#     name = "c",
#     values = c(
#       "P(c = 0)" = "black",
#       "P(c = 1)" = "green",
#       "P(c = 2)" = "blue",
#       "P(c > 2)" = "red"
#     )
#   )
# plot_rknf_prob_c
# 
# # Multiresolution flow
# mrf_results_df <- cbind(delta_vec, mrf_results)
# plot_mrf_prob_reversibility_issue <- ggplot(data = mrf_results_df) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_reversiblity_issue
#     ),
#     color = "black"
#   ) + 
#   geom_ribbon(
#     aes(
#       x = log(delta_vec), 
#       ymin = lower_prob_reversibility_issue,
#       ymax = upper_prob_reversibility_issue
#     ),
#     alpha = 0.15,
#     fill = "red",
#     linetype = 0
#   ) + 
#   ylim(
#     c(0, 0.18)
#   ) + 
#   labs(
#     y = "Probability",
#     x = "log(delta)"
#   )
# plot_mrf_prob_reversibility_issue
# 
# plot_mrf_prob_c <- ggplot(data = mrf_results_df) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_0,
#       color = "P(c = 1)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_1,
#       color = "P(c = 2)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_2,
#       color = "P(c = 3)"
#     )
#   ) +
#   geom_line(
#     aes(
#       x = log(delta_vec), 
#       y = prob_c_larger_2,
#       color = "P(c > 3)"
#     )
#   ) + 
#   labs(
#     y = "Probability",
#     x = "log(delta)"
#   ) + 
#   scale_color_manual(
#     name = "c",
#     values = c(
#       "P(c = 1)" = "black",
#       "P(c = 2)" = "green",
#       "P(c = 3)" = "blue",
#       "P(c > 3)" = "red"
#     )
#   )
# plot_mrf_prob_c
# 
# gridExtra::grid.arrange(
#   plot_hamiltonian_prob_reversibility_issue, plot_hamiltonian_prob_c,
#   plot_rknf_prob_reversibility_issue, plot_rknf_prob_c,
#   plot_mrf_prob_reversibility_issue, plot_mrf_prob_c,
#   nrow = 3, 
#   ncol = 2
# )
# 
# # Alternative: 
# 
# hamiltonian_row <- gridExtra::arrangeGrob(
#   grid::textGrob("Hamiltonian criterion"),
#   gridExtra::arrangeGrob(plot_hamiltonian_prob_reversibility_issue, plot_hamiltonian_prob_c, ncol = 2),
#   ncol = 1,
#   heights = c(0.15, 1)
# )
# 
# rknf_row <- gridExtra::arrangeGrob(
#   grid::textGrob("Runge-Kutta-Nyström flow criterion"),
#   gridExtra::arrangeGrob(plot_rknf_prob_reversibility_issue, plot_rknf_prob_c, ncol = 2),
#   ncol = 1,
#   heights = c(0.15, 1)
# )
# 
# mrf_row <- gridExtra::arrangeGrob(
#   grid::textGrob("Multiresolution flow criterion"),
#   gridExtra::arrangeGrob(plot_mrf_prob_reversibility_issue, plot_mrf_prob_c, ncol = 2),
#   ncol = 1,
#   heights = c(0.15, 1)
# )
# 
# grid.arrange(hamiltonian_row, rknf_row, mrf_row, ncol = 1)

# # Hamiltonian
# 
# hamiltonian_results_df <- cbind(delta_vec, hamiltonian_results)
# plot_hamiltonian_prob_reversibility_issue <- ggplot(data = hamiltonian_results_df) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_reversiblity_issue
#     ),
#     color = "black"
#   ) + 
#   geom_ribbon(
#     aes(
#       x = delta_vec, 
#       ymin = lower_prob_reversibility_issue,
#       ymax = upper_prob_reversibility_issue
#     ),
#     alpha = 0.15,
#     fill = "red",
#     linetype = 0
#   ) + 
#   ylim(
#     c(0, 0.18)
#   ) + 
#   labs(
#     y = "Probability",
#     x = "delta"
#   )
# plot_hamiltonian_prob_reversibility_issue
# 
# plot_hamiltonian_prob_c <- ggplot(data = hamiltonian_results_df) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_0,
#       color = "P(c = 0)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_1,
#       color = "P(c = 1)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_2,
#       color = "P(c = 2)"
#     )
#   ) +
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_larger_2,
#       color = "P(c > 2)"
#     )
#   ) + 
#   labs(
#     y = "Probability",
#     x = "delta"
#   ) + 
#   scale_color_manual(
#     name = "c",
#     values = c(
#       "P(c = 0)" = "black",
#       "P(c = 1)" = "green",
#       "P(c = 2)" = "blue",
#       "P(c > 2)" = "red"
#     )
#   )
# plot_hamiltonian_prob_c
# 
# # RKNF
# rknf_results_df <- cbind(delta_vec, rknf_results)
# plot_rknf_prob_reversibility_issue <- ggplot(data = rknf_results_df) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_reversiblity_issue
#     ),
#     color = "black"
#   ) + 
#   geom_ribbon(
#     aes(
#       x = delta_vec, 
#       ymin = lower_prob_reversibility_issue,
#       ymax = upper_prob_reversibility_issue
#     ),
#     alpha = 0.15,
#     fill = "red",
#     linetype = 0
#   ) + 
#   ylim(
#     c(0, 0.18)
#   ) + 
#   labs(
#     y = "Probability",
#     x = "delta"
#   )
# plot_rknf_prob_reversibility_issue
# 
# plot_rknf_prob_c <- ggplot(data = rknf_results_df) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_0,
#       color = "P(c = 0)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_1,
#       color = "P(c = 1)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_2,
#       color = "P(c = 2)"
#     )
#   ) +
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_larger_2,
#       color = "P(c > 2)"
#     )
#   ) + 
#   labs(
#     y = "Probability",
#     x = "delta"
#   ) + 
#   scale_color_manual(
#     name = "c",
#     values = c(
#       "P(c = 0)" = "black",
#       "P(c = 1)" = "green",
#       "P(c = 2)" = "blue",
#       "P(c > 2)" = "red"
#     )
#   )
# plot_rknf_prob_c
# 
# # Multiresolution flow
# mrf_results_df <- cbind(delta_vec, mrf_results)
# plot_mrf_prob_reversibility_issue <- ggplot(data = mrf_results_df) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_reversiblity_issue
#     ),
#     color = "black"
#   ) + 
#   geom_ribbon(
#     aes(
#       x = delta_vec, 
#       ymin = lower_prob_reversibility_issue,
#       ymax = upper_prob_reversibility_issue
#     ),
#     alpha = 0.15,
#     fill = "red",
#     linetype = 0
#   ) + 
#   ylim(
#     c(0, 0.18)
#   ) + 
#   labs(
#     y = "Probability",
#     x = "delta"
#   )
# plot_mrf_prob_reversibility_issue
# 
# plot_mrf_prob_c <- ggplot(data = mrf_results_df) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_0,
#       color = "P(c = 1)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_1,
#       color = "P(c = 2)"
#     )
#   ) + 
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_2,
#       color = "P(c = 3)"
#     )
#   ) +
#   geom_line(
#     aes(
#       x = delta_vec, 
#       y = prob_c_larger_2,
#       color = "P(c > 3)"
#     )
#   ) + 
#   labs(
#     y = "Probability",
#     x = "delta"
#   ) + 
#   scale_color_manual(
#     name = "c",
#     values = c(
#       "P(c = 1)" = "black",
#       "P(c = 2)" = "green",
#       "P(c = 3)" = "blue",
#       "P(c > 3)" = "red"
#     )
#   )
# plot_mrf_prob_c
# 
# gridExtra::grid.arrange(
#   plot_hamiltonian_prob_reversibility_issue, plot_hamiltonian_prob_c,
#   plot_rknf_prob_reversibility_issue, plot_rknf_prob_c,
#   plot_mrf_prob_reversibility_issue, plot_mrf_prob_c,
#   nrow = 3, 
#   ncol = 2
# )
# 
# # Alternative: 
# 
# hamiltonian_row <- gridExtra::arrangeGrob(
#   grid::textGrob("Hamiltonian criterion"),
#   gridExtra::arrangeGrob(plot_hamiltonian_prob_reversibility_issue, plot_hamiltonian_prob_c, ncol = 2),
#   ncol = 1,
#   heights = c(0.15, 1)
# )
# 
# rknf_row <- gridExtra::arrangeGrob(
#   grid::textGrob("Runge-Kutta-Nyström flow criterion"),
#   gridExtra::arrangeGrob(plot_rknf_prob_reversibility_issue, plot_rknf_prob_c, ncol = 2),
#   ncol = 1,
#   heights = c(0.15, 1)
# )
# 
# mrf_row <- gridExtra::arrangeGrob(
#   grid::textGrob("Multiresolution flow criterion"),
#   gridExtra::arrangeGrob(plot_mrf_prob_reversibility_issue, plot_mrf_prob_c, ncol = 2),
#   ncol = 1,
#   heights = c(0.15, 1)
# )
# 
# grid.arrange(hamiltonian_row, rknf_row, mrf_row, ncol = 1)
