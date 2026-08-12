# Leapfrog functions 

# One single leapfrog step with step size defined as eps

single_leapfrog_integrator_step_identity_matrix <- function(grad_log_target, eps, theta, rho, grad_log_target_initial = NULL) {
  
  if (is.null(grad_log_target_initial)) {
    # Evaluate the gradient of log target at initial state if this is not given before hand 
    # This part is included in order to evaluate the gradient only once per leapfrog step when composing several steps by reusing the former evaluated gradient
    grad_log_target_initial <- grad_log_target(theta)
  }
  
  theta_prime <- theta + (
    eps * rho + (eps ^ 2) / 2 * grad_log_target_initial # update theta according to leapfrog 
  )
  
  grad_log_target_prime <- grad_log_target(theta_prime)
  
  rho_prime <- rho + eps / 2 * (grad_log_target_initial + grad_log_target_prime) # update rho according to leapfrog - need this step now since the proposed state is composed by several leapfrog steps
  
  list(theta_prime = theta_prime, rho_prime = rho_prime, grad_log_target_prime = grad_log_target_prime)
  
}

# Several leapfrog step with step size defined as eps and number of leapfrog steps taken with this step size as k

several_leapfrog_integrator_steps_identity_matrix <- function(grad_log_target, eps, k, theta, rho, grad_log_target_initial = NULL) {
  current_rho <- rho
  current_theta <- theta
  if (!is.null(grad_log_target_initial)) { # similar as above
    current_grad_log_target <- grad_log_target_initial
  } else {
    current_grad_log_target <- grad_log_target(theta)
  }
  for (j in 1:k) {
    one_leapfrog_step <- single_leapfrog_integrator_step_identity_matrix(grad_log_target, eps, current_theta, current_rho, current_grad_log_target) # run one single leapfrog step
    current_theta <- one_leapfrog_step$theta_prime # update initial state
    current_rho <- one_leapfrog_step$rho_prime
    current_grad_log_target <- one_leapfrog_step$grad_log_target_prime
  }
  list(theta_prime = current_theta, rho_prime = current_rho, grad_log_target_prime = current_grad_log_target)
}
