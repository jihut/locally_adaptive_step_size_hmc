d <- 10

grad_log_target_fun <- function(theta) {
  c(
    -1 / 9 * theta[1] - (d - 1) / 2 + 1 / (2 * exp(theta[1])) * sum(theta[2:d] ^ 2),
    -theta[2:d] / exp(theta[1])
  )
}

hamiltonian_func <- function(theta, rho) {
  # -mvtnorm::dmvnorm(theta, mean = mu, sigma = sigma, log = T) + 0.5 * sum(rho ^ 2)
  -(-1 / (2 * 9) * (theta[1] ^ 2) - (d - 1) / 2 * theta[1] - 1 / (2 * exp(theta[1])) * sum(theta[2:d] ^ 2)) + 0.5 * sum(rho ^ 2)
} 

log_target_fun <- function(theta)   -1 / (2 * 9) * (theta[1] ^ 2) - (d - 1) / 2 * theta[1] - 1 / (2 * exp(theta[1])) * sum(theta[2:d] ^ 2)
