source("implementation_scripts/general_scripts/tuning_criterion_functions.R")

single_iteration_one_step_GHMC <- function(
    micro_fun, 
    h,
    nu, 
    delta,
    max_c,
    log_target,
    grad_log_target,
    hamiltonian,
    theta, 
    rho, 
    d,
    grad_log_target_initial = NULL, 
    norm = "L-infinity"
)

{
  
  current_rho <- rho
  current_theta <- theta
  if (!is.null(grad_log_target_initial)) {
    current_grad_log_target <- grad_log_target_initial
    n_evals_ode <- 0
  } else {
    grad_log_target_initial <- current_grad_log_target <- grad_log_target(theta)
    n_evals_ode <- 1
  }
  
  weight_non_zero_indicator <- 1
  c <- -1
  gamma_obs <- -1
  # Partial momentum refresh
  refreshed_rho <- nu * rho + sqrt(1 - nu ^ 2) * rnorm(d)
  
  forward_run_micro <- micro_fun(
    h = h, 
    delta = delta,
    max_c = max_c,
    log_target = log_target,
    grad_log_target = grad_log_target,
    hamiltonian = hamiltonian,
    theta = current_theta,
    rho = refreshed_rho, 
    grad_log_target_initial = current_grad_log_target,
    norm = norm,
    forward_micro_run_indicator = TRUE
  )
  n_evals_ode <- n_evals_ode + forward_run_micro$n_evals_ode
  
  if (forward_run_micro$l == 1) {
    backward_run_micro_l <- 1
  } else {
    # Check if running backwards based from the resulting state from above, if the same number of micro steps is the smallest to achieve same criterion
    
    backward_run_micro <- micro_fun(
      h = h, 
      delta = delta,
      max_c = forward_run_micro$c - 1, # Lemma 1 in WALNUTS article - for the same value of delta, l_backward will always be less than or equal to l_forward --> only need to check if the criterion is satisfied backward for l smaller than l_forward
      log_target = log_target_fun,
      grad_log_target = grad_log_target_fun,
      hamiltonian = hamiltonian, 
      theta = forward_run_micro$theta,
      rho = -forward_run_micro$rho, 
      grad_log_target_initial = forward_run_micro$grad_log_target,
      norm = norm,
      forward_micro_run_indicator = FALSE,
      theta_final_backward = forward_run_micro$theta_final_backward,
      rho_final_backward = forward_run_micro$rho_final_backward,
      grad_log_target_final_backward = forward_run_micro$grad_log_target_final_backward
    )
    n_evals_ode <- n_evals_ode + backward_run_micro$n_evals_ode
    
    backward_run_micro_l <- ifelse(backward_run_micro$max_c_indicator == 1, forward_run_micro$l, backward_run_micro$l)
  }
  
  if (forward_run_micro$l != backward_run_micro_l) { 
    # if this is not the case, due to deterministic variant of micro step distribution
    # --> Final acceptance probability of the final state after L macro steps will be zero
    # --> Need to throw away the whole trajectory at the end anyways, therefore stop the integration now to save computations
    weight_non_zero_indicator <- 0
  } else {
    current_theta <- forward_run_micro$theta
    current_rho <- forward_run_micro$rho
    current_grad_log_target <- forward_run_micro$grad_log_target
    c <- forward_run_micro$c
  }
  
  if (weight_non_zero_indicator == 0) {
    # If reject, one must apply a momentum flip to the initial state prior to the macro step due to the involution of GHMC
    final_theta <- theta
    final_rho <- -refreshed_rho
    final_grad_log_target <- grad_log_target_initial
    accept_indicator <- 0
    diff_hamiltonian <- 0L
    
  } else { # otherwise perform standard Metropolis to check if one should accept the final state after L macro steps. 
    hamiltonian_initial <- hamiltonian(theta, refreshed_rho)
    hamiltonian_final <- hamiltonian(current_theta, current_rho)
    diff_hamiltonian <- hamiltonian_initial - hamiltonian_final
    u <- runif(1)
    acc_prob <- min(1, exp(diff_hamiltonian))
    
    # The proposal is actually supposed to be 
    # c(new_state$theta, -new_state$rho) to get an involution. 
    # After the Metropolis step, one applies a momentum flip again. 
    # Therefore, if accept, this corresponds to keeping the same sign on the momentum from the leapfrog. 
    if (u < acc_prob) {
      final_theta <- forward_run_micro$theta
      final_rho <- forward_run_micro$rho
      final_grad_log_target <- forward_run_micro$grad_log_target
      accept_indicator <- 1
    } else {
      # If reject, one must apply a momentum flip to the initial state prior to the leapfrog step due to the above
      final_theta <- theta
      final_rho <- -refreshed_rho
      final_grad_log_target <- grad_log_target_initial
      accept_indicator <- 0
    }
    
  }
  
  final_state <- list(
    theta = final_theta,
    rho = final_rho, 
    updated_rho = refreshed_rho, 
    grad_log_target = final_grad_log_target,
    weight_non_zero_indicator = weight_non_zero_indicator,
    c = c,
    accept_indicator = accept_indicator,
    n_evals_ode = n_evals_ode,
    diff_hamiltonian = diff_hamiltonian,
    gamma_obs = forward_run_micro$gamma_obs
  )
  
}

adaptive_step_size_one_step_GHMC_warmup <- function(
    micro_fun, 
    n_samples, 
    n_samples_before_adaptive = 100,
    h,
    nu, 
    delta,
    max_c,
    log_target,
    grad_log_target,
    hamiltonian,
    theta, 
    rho, 
    grad_log_target_initial = NULL, 
    norm = "L-infinity",
    target_acc_prob = 0.8,
    prob_no_step_size_halving = 0.8
) {
  
  d <- length(theta) # dimension of the parameter vector
  samples_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the actual samples
  samples_matrix[1, ] <- theta
  
  initial_theta_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial theta at each iteration
  initial_rho_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial rho/momentum at each iteration
  updated_rho_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the partially refreshed momentum before taking a micro step
  
  accept_indicator <- numeric(n_samples)
  reversibility_indicator <- numeric(n_samples)
  num_forward_steps_vector <- numeric(n_samples)
  diff_hamiltonian_vector <- numeric(n_samples)
  
  h_vector <- numeric(n_samples)
  delta_vector <- numeric(n_samples)
  
  current_theta <- theta
  current_grad_log_target <- grad_log_target(theta)
  n_evals_ode <- 1
  current_delta <- delta
  current_h <- h
  
  hamiltonian <- function(theta, rho) {
    -log_target(theta) + 0.5 * sum(rho ^ 2)
  }
  diff_hamiltonian_threshold <- log(target_acc_prob)
  gamma_obs_vector <- numeric(n_samples)
  # kappa_vector <- numeric(n_samples_between_updates)
  kappa_vector <- numeric(n_samples)
  n_updates <- 1
  count_samples_between_updates <- 0
  current_rho <- rnorm(d)
  
  for (i in 1:n_samples) {
    # if (i %% 1000 == 0) print(i)
    print(i)
    initial_rho_matrix[i, ] <- current_rho
    initial_theta_matrix[i, ] <- current_theta
    new_state <- single_iteration_one_step_GHMC(
      micro_fun = micro_fun, 
      h = current_h,
      nu = nu, 
      delta = current_delta,
      max_c = max_c,
      log_target = log_target,
      grad_log_target = grad_log_target,
      hamiltonian = hamiltonian,
      theta = current_theta, 
      rho = current_rho, 
      d = d,
      grad_log_target_initial = current_grad_log_target, 
      norm = "L-infinity"
    )
    n_evals_ode <- n_evals_ode + new_state$n_evals_ode
    samples_matrix[i, ] <- current_theta <- new_state$theta
    accept_indicator[i] <- new_state$accept_indicator
    current_grad_log_target <- new_state$grad_log_target
    current_rho <- new_state$rho
    updated_rho <- new_state$updated_rho
    updated_rho_matrix[i, ] <- updated_rho
    reversibility_indicator[i] <- new_state$weight_non_zero_indicator
    accept_indicator[i] <- new_state$accept_indicator
    diff_hamiltonian_vector[i] <- new_state$diff_hamiltonian
    gamma_obs_vector[i] <- new_state$gamma_obs
    kappa_vector[i] <- new_state$diff_hamiltonian / current_delta
    
    if (i >= n_samples_before_adaptive) {
      current_kappa_vector <- kappa_vector[1:i]
      current_kappa_vector <- current_kappa_vector[current_kappa_vector != 0L]
      # current_delta <- (diff_hamiltonian_threshold / quantile(current_kappa_vector, probs = target_acc_prob))
      current_delta <- (diff_hamiltonian_threshold / quantile(current_kappa_vector, probs = 1 - target_acc_prob))
      # current_log_h <- 1 / 3 * log(current_delta / quantile(gamma_obs_vector[1:i], probs = prob_no_step_size_halving))
      # current_h <- exp(current_log_h)
      current_h <- (current_delta / quantile(gamma_obs_vector[1:i], probs = prob_no_step_size_halving)) ^ (1 / 3)
      if (i %% 1000 == 0) {
        print(paste0("i: ", i))
        print(paste0("h: ", current_h))
        print(paste0("delta: ", current_delta))
      }
    }
  }
  
  list(
    samples_matrix = samples_matrix, 
    initial_theta_matrix = initial_theta_matrix, 
    initial_rho_matrix = initial_rho_matrix,
    updated_rho_matrix = updated_rho_matrix, 
    diff_hamiltonian_vector = diff_hamiltonian_vector,
    reversibility_indicator = reversibility_indicator,
    accept_indicator = accept_indicator,
    final_h = current_h,
    final_delta = current_delta, 
    kappa_vector = kappa_vector,
    gamma_obs_vector = gamma_obs_vector,
    n_evals_ode = n_evals_ode
  )
  
}

adaptive_step_size_one_step_GHMC_sampling <- function(
    micro_fun, 
    n_samples, 
    h,
    nu, 
    delta,
    max_c,
    log_target,
    grad_log_target,
    theta, 
    norm = "L-infinity",
    relevant_indices = NULL # indices of the parameters to store in case one does not want to store all dimensions
) {
  
  d <- length(theta) # dimension of the parameter vector
  
  if (is.null(relevant_indices)) {
    relevant_indices <- 1:d
    samples_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the actual samples
    samples_matrix[1, ] <- theta
    
    initial_theta_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial theta at each iteration
    initial_rho_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial rho/momentum at each iteration
    updated_rho_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the partially refreshed momentum before taking a micro step 
  } else {
    d_relevant <- length(relevant_indices)
    samples_matrix <- matrix(nrow = n_samples, ncol = d_relevant) # matrix to store the actual samples
    samples_matrix[1, ] <- theta[relevant_indices]
    
    initial_theta_matrix <- matrix(nrow = n_samples, ncol = d_relevant) # matrix to store the initial theta at each iteration
    initial_rho_matrix <- matrix(nrow = n_samples, ncol = d_relevant) # matrix to store the initial rho/momentum at each iteration
    updated_rho_matrix <- matrix(nrow = n_samples, ncol = d_relevant) # matrix to store the partially refreshed momentum before taking a micro step 
  }
  
  accept_indicator <- numeric(n_samples)
  reversibility_indicator <- numeric(n_samples)
  c_vector <- numeric(n_samples)
  diff_hamiltonian_vector <- numeric(n_samples)
  
  h_vector <- numeric(n_samples)
  delta_vector <- numeric(n_samples)
  
  current_theta <- theta
  current_grad_log_target <- grad_log_target(theta)
  n_evals_ode <- 1
  
  hamiltonian <- function(theta, rho) {
    -log_target(theta) + 0.5 * sum(rho ^ 2)
  }
  
  current_rho <- rnorm(d)
  
  for (i in 1:n_samples) {
    # if (i %% 1000 == 0) print(i)
    print(i)
    initial_rho_matrix[i, ] <- current_rho[relevant_indices]
    initial_theta_matrix[i, ] <- current_theta[relevant_indices]
    new_state <- single_iteration_one_step_GHMC(
      micro_fun = micro_fun, 
      h = h,
      nu = nu, 
      delta = delta,
      max_c = max_c,
      log_target = log_target,
      grad_log_target = grad_log_target,
      hamiltonian = hamiltonian,
      theta = current_theta, 
      rho = current_rho, 
      d = d,
      grad_log_target_initial = current_grad_log_target, 
      norm = "L-infinity"
    )
    n_evals_ode <- n_evals_ode + new_state$n_evals_ode
    samples_matrix[i, ] <- new_state$theta[relevant_indices]
    current_theta <- new_state$theta
    accept_indicator[i] <- new_state$accept_indicator
    current_grad_log_target <- new_state$grad_log_target
    current_rho <- new_state$rho
    updated_rho <- new_state$updated_rho
    updated_rho_matrix[i, ] <- updated_rho[relevant_indices]
    reversibility_indicator[i] <- new_state$weight_non_zero_indicator
    accept_indicator[i] <- new_state$accept_indicator
    diff_hamiltonian_vector[i] <- new_state$diff_hamiltonian
    c_vector[i] <- new_state$c
    
  }
  
  list(
    samples_matrix = samples_matrix, 
    initial_theta_matrix = initial_theta_matrix, 
    initial_rho_matrix = initial_rho_matrix,
    updated_rho_matrix = updated_rho_matrix, 
    diff_hamiltonian_vector = diff_hamiltonian_vector,
    reversibility_indicator = reversibility_indicator,
    accept_indicator = accept_indicator,
    c_vector = c_vector, 
    n_evals_ode = n_evals_ode,
    final_theta = current_theta
  )
}
