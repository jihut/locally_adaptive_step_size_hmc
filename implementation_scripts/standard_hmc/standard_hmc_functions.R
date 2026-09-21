source("implementation_scripts/general_scripts/tuning_criterion_functions.R")

single_iteration_standard_HMC <- function(
    micro_fun, h, L, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial = NULL, norm = "L-infinity",
    return_c_for_L_equal_to_1 = FALSE
) {
  
  current_rho <- rho
  current_theta <- theta
  if (!is.null(grad_log_target_initial)) {
    current_grad_log_target <- grad_log_target_initial
    n_evals_ode <- 0
  } else {
    grad_log_target_initial <- current_grad_log_target <- grad_log_target(theta)
    n_evals_ode <- 1
  }
  c_vector <- rep(-1, L)
  c_vector_index <- 1
  weight_non_zero_indicator <- 1
  for(i in 1:L) { # run L macro steps of size h to get the final proposed state
    # For each macro step of size h, need to find the number of micro steps to achieve the desired error criterion
    forward_run_micro <- micro_fun(
      h = h, 
      delta = delta,
      max_c = max_c,
      log_target = log_target,
      grad_log_target = grad_log_target,
      hamiltonian = hamiltonian,
      theta = current_theta,
      rho = current_rho, 
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
      if (L == 1 & return_c_for_L_equal_to_1 == TRUE) {
        c_vector[1] <- forward_run_micro$c
      }
      weight_non_zero_indicator <- 0
      break
    } else {
      current_theta <- forward_run_micro$theta
      current_rho <- forward_run_micro$rho
      current_grad_log_target <- forward_run_micro$grad_log_target
      c_vector[c_vector_index] <- forward_run_micro$c
      c_vector_index <- c_vector_index + 1
    }
  }
  
  if (weight_non_zero_indicator == 0) { # return initial state if one of the macro steps leads to a situation where the minimum number of micro steps is not the same forward and backward
    
    final_theta <- theta
    final_rho <- rho
    final_grad_log_target <- grad_log_target_initial
    accept_indicator <- 0
    
    
  } else { # otherwise perform standard Metropolis to check if one should accept the final state after L macro steps. 
    hamiltonian_initial <- hamiltonian(theta, rho)
    hamiltonian_final <- hamiltonian(current_theta, current_rho)
    
    u <- runif(1)
    acc_prob <- min(1, exp(hamiltonian_initial - hamiltonian_final))
    
    if (u < acc_prob) {
      
      final_theta <- current_theta
      final_rho <- current_rho
      final_grad_log_target <- current_grad_log_target
      accept_indicator <- 1
      
      
    } else {
      
      final_theta <- theta
      final_rho <- rho
      final_grad_log_target <- grad_log_target_initial
      accept_indicator <- 0
      
    }
    
  }
  
  final_state <- list(theta = final_theta,
                      rho = final_rho,
                      grad_log_target = final_grad_log_target,
                      weight_non_zero_indicator = weight_non_zero_indicator,
                      c_vector = c_vector, 
                      max_c = max(c_vector[1:(c_vector_index - 1)]),
                      min_c = min(c_vector[1:(c_vector_index - 1)]),
                      accept_indicator = accept_indicator,
                      n_evals_ode = n_evals_ode)
  
  final_state
  
}

# N iterations of standard HMC where each iteration is based on L macro steps of step size h with adaptive micro step size

adaptive_step_size_standard_HMC <- function(
    micro_fun,
    n_samples, 
    h, 
    L, 
    delta, 
    max_c, 
    log_target, 
    grad_log_target, 
    theta,
    norm = "L-infinity",
    return_c_for_L_equal_to_1 = FALSE,
    relevant_indices = NULL, # indices of the parameters to store in case one does not want to store all dimensions
    randomized_L = FALSE # indicates if the number of leapfrog steps is randomized according to a Poisson distribution with expectation equal to L or not
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
  
  max_c_vector <- numeric(n_samples)
  min_c_vector <- numeric(n_samples)
  reversibility_indicator <- numeric(n_samples) 
  accept_indicator <- numeric(n_samples)
  if (randomized_L) {
    c_matrix <- matrix(nrow = n_samples, ncol = 5 * L) # safeguard in case one gets L that is much larger than the mean
    L_vector <- rpois(n = n_samples, lambda = L)
    L_vector <- ifelse(L_vector == 0, 1, L_vector) # at least one leapfrog step needs tbe taken
  } else {
    c_matrix <- matrix(nrow = n_samples, ncol = L)
  }

  current_theta <- theta
  current_grad_log_target <- grad_log_target(theta)
  n_evals_ode <- 1
  
  hamiltonian <- function(theta, rho) {
    -log_target(theta) + 0.5 * sum(rho ^ 2)
  }
  
  for (i in 1:n_samples) {
    # if (i %% 1000 == 0) print(i)
    print(i)
    current_rho <- rnorm(d)
    initial_rho_matrix[i, ] <- current_rho[relevant_indices]
    initial_theta_matrix[i, ] <- current_theta[relevant_indices]
    if (randomized_L) {
      new_state <- single_iteration_standard_HMC( # obtain one single iteration using fixed length Hamiltonian with adaptive step size
        micro_fun = micro_fun, 
        h = h, 
        L = L_vector[i], 
        delta = delta, 
        max_c = max_c, 
        log_target = log_target_fun, 
        grad_log_target = grad_log_target_fun, 
        hamiltonian = hamiltonian,
        theta = current_theta, 
        rho = current_rho, 
        grad_log_target_initial = current_grad_log_target,
        norm = norm,
        return_c_for_L_equal_to_1 = return_c_for_L_equal_to_1
      )
      c_matrix[i, 1:L_vector[i]] <- new_state$c_vector
    } else {
      new_state <- single_iteration_standard_HMC( # obtain one single iteration using fixed length Hamiltonian with adaptive step size
        micro_fun = micro_fun, 
        h = h, 
        L = L, 
        delta = delta, 
        max_c = max_c, 
        log_target = log_target_fun, 
        grad_log_target = grad_log_target_fun, 
        hamiltonian = hamiltonian,
        theta = current_theta, 
        rho = current_rho, 
        grad_log_target_initial = current_grad_log_target,
        norm = norm,
        return_c_for_L_equal_to_1 = return_c_for_L_equal_to_1
      )
      c_matrix[i, ] <- new_state$c_vector
    }
    n_evals_ode <- n_evals_ode + new_state$n_evals_ode
    samples_matrix[i, ] <- new_state$theta[relevant_indices]
    current_theta <- new_state$theta
    max_c_vector[i] <- new_state$max_c
    min_c_vector[i] <- new_state$min_c
    accept_indicator[i] <- new_state$accept_indicator
    reversibility_indicator[i] <- new_state$weight_non_zero_indicator
    current_grad_log_target <- new_state$grad_log_target
  }
  
  if (randomized_L) {
    list(samples_matrix = samples_matrix, initial_theta_matrix = initial_theta_matrix, initial_rho_matrix = initial_rho_matrix,
         max_c_vector = max_c_vector, min_c_vector = min_c_vector, reversibility_indicator = reversibility_indicator,
         accept_indicator = accept_indicator,
         c_matrix = c_matrix, 
         L_vector = L_vector,
         n_evals_ode = n_evals_ode,
         final_theta = current_theta)  
  } else {
    list(samples_matrix = samples_matrix, initial_theta_matrix = initial_theta_matrix, initial_rho_matrix = initial_rho_matrix,
         max_c_vector = max_c_vector, min_c_vector = min_c_vector, reversibility_indicator = reversibility_indicator,
         accept_indicator = accept_indicator,
         c_matrix = c_matrix, 
         n_evals_ode = n_evals_ode,
         final_theta = current_theta)
  }
  
}
