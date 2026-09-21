source("implementation_scripts/general_scripts/tuning_criterion_functions.R")

single_iteration_standard_biased_progressive_HMC <- function(
    micro_fun, h, m, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial = NULL, d, norm = "L-infinity"
) {
  if (!is.null(grad_log_target_initial)) {
    current_grad_log_target <- grad_log_target_initial
    n_evals_ode <- 0
  } else {
    current_grad_log_target <- grad_log_target(theta)
    n_evals_ode <- 1
  }
  
  n_states_in_orbit <- 2 ^ m
  B_vector <- sample(0:1, m, replace = T)
  b <- sum(B_vector * 2 ^ (1:m - 1)) # rightmost index
  orbit_matrix <- matrix(c(theta, rho), nrow = 1, ncol = 2 * d)
  grad_log_target_matrix <- matrix(current_grad_log_target, nrow = 1, ncol = d)
  hamiltonian_vector <- hamiltonian_initial <- hamiltonian(theta, rho)
  log_weight_vector <- log_weight_initial <- -hamiltonian_vector
  c_vector <- c(-1)
  gamma_obs_vector <- c(-1)
  final_state <- orbit_matrix[1, ]
  final_grad_log_target <- grad_log_target_matrix[1, ]
  final_index <- 0
  current_a <- 0
  current_b <- 0
  backwards_yet_indicator <- 0
  
  dead_forward_indicator <- 0 # if one of the weights when doing forward steps becomes zero because the minimum number of micro steps forward and backward to achieve same criterion does not coincide
  dead_backward_indicator <- 0 # similar as above, but for backward
  
  for (k in (0:(m - 1))) {
    # if (dead_forward_indicator == 1 & dead_forward_indicator == 1) break # dead ends, weights beyond will be zero
    num_of_steps <- 2 ^ k
    current_orbit_matrix <- matrix(0, nrow = num_of_steps, ncol = 2 * d)
    current_grad_log_target_matrix <- matrix(0, nrow = num_of_steps, ncol = d)
    # current_hamiltonian_vector <- numeric(num_of_steps)
    current_hamiltonian_vector <- rep(NA, num_of_steps)
    current_log_weight_vector <- rep(-1e12, num_of_steps)
    current_c_vector <- rep(-1, num_of_steps)
    current_gamma_obs_vector <- rep(-1, num_of_steps)
    
    if (B_vector[k + 1] == 1) { # forward
      
      if (dead_forward_indicator == 0) {
        
        current_state <- tail(orbit_matrix, 1)[1, ]
        current_theta <- current_state[1:d]
        current_rho <- current_state[(d + 1):(2 * d)]
        current_grad_log_target <- tail(grad_log_target_matrix, 1)[1, ]
        
        for (j in 1:num_of_steps) {
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
          current_gamma_obs_vector[j] <- forward_run_micro$gamma_obs
          n_evals_ode <- n_evals_ode + forward_run_micro$n_evals_ode
          
          if (forward_run_micro$l == 1) {
            backward_run_micro_l <- 1
          } else {
            backward_run_micro <- micro_fun(
              h = h, 
              delta = delta,
              max_c = forward_run_micro$c - 1, # Lemma 1 in WALNUTS article - for the same value of delta, l_backward will always be less than or equal to l_forward --> only need to check if the criterion is satisfied backward for l smaller than l_forward
              log_target = log_target_fun,
              grad_log_target = grad_log_target,
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
            dead_forward_indicator <- 1
            break
          } else {
            current_orbit_matrix[j, ] <- c(forward_run_micro$theta, forward_run_micro$rho)
            current_grad_log_target_matrix[j, ] <- forward_run_micro$grad_log_target
            current_hamiltonian_vector[j] <- hamiltonian(forward_run_micro$theta, forward_run_micro$rho)
            current_log_weight_vector[j] <- -current_hamiltonian_vector[j]
            current_c_vector[j] <- forward_run_micro$c
            
            current_theta <- forward_run_micro$theta
            current_rho <- forward_run_micro$rho
            current_grad_log_target <- forward_run_micro$grad_log_target
          }
          
        }
        
        current_max_log_weight <- max(current_log_weight_vector)
        scaled_weight_vector_extension <- exp(current_log_weight_vector - current_max_log_weight)
        sample_index_extension <- sample(1:num_of_steps, 1, replace = T, scaled_weight_vector_extension / sum(scaled_weight_vector_extension))
        metropolis_acceptance <- min(1, sum(exp(current_log_weight_vector - current_max_log_weight)) / 
                                       sum(exp(log_weight_vector - current_max_log_weight)))
        if (runif(1) <= metropolis_acceptance) {
          final_state <- current_orbit_matrix[sample_index_extension, ]
          final_grad_log_target <- current_grad_log_target_matrix[sample_index_extension, ]
          final_index <- current_b + sample_index_extension
        }
        
      }
      
      orbit_matrix <- rbind(orbit_matrix, current_orbit_matrix)
      grad_log_target_matrix <- rbind(grad_log_target_matrix, current_grad_log_target_matrix)
      hamiltonian_vector <- c(hamiltonian_vector, current_hamiltonian_vector)
      log_weight_vector <- c(log_weight_vector, current_log_weight_vector)
      c_vector <- c(c_vector, current_c_vector)
      gamma_obs_vector <- c(gamma_obs_vector, current_gamma_obs_vector)
      current_b <- current_b + num_of_steps
      
    } else if (B_vector[k + 1] == 0) { # backward
      
      if (dead_backward_indicator == 0) {
        
        current_state <- head(orbit_matrix, 1)[1, ]
        current_theta <- current_state[1:d]
        if (backwards_yet_indicator == 0) {
          current_rho <- -current_state[(d + 1):(2 * d)] # flip momentum to go backwards
          backwards_yet_indicator <- 1
        } else {
          current_rho <- current_state[(d + 1):(2 * d)]
        }
        current_grad_log_target <- head(grad_log_target_matrix, 1)[1, ]
        
        for (j in 1:num_of_steps) {
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
          current_gamma_obs_vector[num_of_steps + 1 - j] <- forward_run_micro$gamma_obs
          n_evals_ode <- n_evals_ode + forward_run_micro$n_evals_ode
          
          if (forward_run_micro$l == 1) {
            backward_run_micro_l <- 1
          } else {
            backward_run_micro <- micro_fun(
              h = h, 
              delta = delta,
              max_c = forward_run_micro$c - 1, # Lemma 1 in WALNUTS article - for the same value of delta, l_backward will always be less than or equal to l_forward --> only need to check if the criterion is satisfied backward for l smaller than l_forward
              log_target = log_target_fun,
              grad_log_target = grad_log_target,
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
            dead_backward_indicator <- 1
            break
          } else {
            current_orbit_matrix[num_of_steps + 1 - j, ] <- c(forward_run_micro$theta, forward_run_micro$rho)
            current_grad_log_target_matrix[num_of_steps + 1 - j, ] <- forward_run_micro$grad_log_target
            current_hamiltonian_vector[num_of_steps + 1 - j] <- hamiltonian(forward_run_micro$theta, forward_run_micro$rho)
            current_log_weight_vector[num_of_steps + 1 - j] <- -current_hamiltonian_vector[num_of_steps + 1 - j]
            current_c_vector[num_of_steps + 1 - j] <- forward_run_micro$c
            
            current_theta <- forward_run_micro$theta
            current_rho <- forward_run_micro$rho
            current_grad_log_target <- forward_run_micro$grad_log_target
          }
          
        }
        
        current_max_log_weight <- max(current_log_weight_vector)
        scaled_weight_vector_extension <- exp(current_log_weight_vector - current_max_log_weight)
        sample_index_extension <- sample(num_of_steps:1, 1, replace = T, scaled_weight_vector_extension / sum(scaled_weight_vector_extension))
        metropolis_acceptance <- min(1, sum(exp(current_log_weight_vector - current_max_log_weight)) / 
                                       sum(exp(log_weight_vector - current_max_log_weight)))
        if (runif(1) <= metropolis_acceptance) {
          final_state <- current_orbit_matrix[num_of_steps + 1 - sample_index_extension, ]
          final_grad_log_target <- current_grad_log_target_matrix[num_of_steps + 1 - sample_index_extension, ]
          final_index <- current_a - sample_index_extension
        } 
        
      }
      
      orbit_matrix <- rbind(current_orbit_matrix, orbit_matrix)
      grad_log_target_matrix <- rbind(current_grad_log_target_matrix, grad_log_target_matrix)
      hamiltonian_vector <- c(current_hamiltonian_vector, hamiltonian_vector)
      log_weight_vector <- c(current_log_weight_vector, log_weight_vector)
      c_vector <- c(current_c_vector, c_vector)
      gamma_obs_vector <- c(current_gamma_obs_vector, gamma_obs_vector)
      current_a <- current_a - num_of_steps
      
    }
    
  }
  
  initial_index_orbit_matrix <- 2 ^ m - b
  c_vector_max <- max(c_vector)
  c_vector_min <- min(c_vector[c_vector != -1])
  
  list(
    B_vector = B_vector, 
    b = b,
    a = b - n_states_in_orbit + 1,
    final_index = final_index, 
    theta = final_state[1:d],
    rho = final_state[(d + 1):(2 * d)],
    grad_log_target = final_grad_log_target,
    orbit_matrix = orbit_matrix,
    log_weight_vector = log_weight_vector,
    hamiltonian_vector = hamiltonian_vector,
    c_vector = c_vector,
    dead_forward_indicator = dead_forward_indicator,
    dead_backward_indicator = dead_backward_indicator,
    max_c = ifelse(is.finite(c_vector_max), c_vector_max, -1),
    min_c = ifelse(is.finite(c_vector_min), c_vector_min, -1),
    # max_orbit_energy_error = max(hamiltonian_vector[hamiltonian_vector != 0]) - min(hamiltonian_vector[hamiltonian_vector != 0]),
    max_orbit_energy_error = max(hamiltonian_vector, na.rm = T) - min(hamiltonian_vector, na.rm = T),
    n_evals_ode = n_evals_ode,
    gamma_obs_vector = gamma_obs_vector[gamma_obs_vector != -1]
  )
  
}

adaptive_step_size_biased_progressive_HMC_warmup <- function(
    micro_fun,
    n_samples, 
    # n_samples_between_updates, 
    # moving_avg_indicator = FALSE, 
    # moving_avg_parameter = 5, 
    n_samples_before_adaptive = 100, 
    h, 
    m, 
    delta, 
    max_c, 
    log_target, 
    grad_log_target, 
    theta,
    norm = "L-infinity",
    prob_no_step_size_halving = 0.8,
    prob_max_orbit_energy_error_below_threshold = 0.8,
    max_orbit_energy_error_threshold = 0.2
) {
  
  # if (n_samples %% n_samples_between_updates != 0) {
  #   stop("Please provide n_samples and n_samples_between_updates so that the ratio is an integer!")
  # }
  # 
  # if (n_samples / n_samples_between_updates < moving_avg_parameter & moving_avg_indicator) {
  #   stop("The number of updates is less than the moving average parameter!")
  # }
  # 
  # if (n_samples_before_adaptive %% n_samples_between_updates != 0) {
  #   stop("Please provide n_samples_before_adaptive and n_samples_between_updates so that the ratio is an integer!")
  # }
  # 
  # n_updates <- n_samples / n_samples_between_updates
  L <- (2 ^ m) - 1 
  
  path_length <- h * L
  d <- length(theta) # dimension of the parameter vector
  samples_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the actual samples
  samples_matrix[1, ] <- theta
  
  initial_theta_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial theta at each iteration
  initial_rho_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial rho/momentum at each iteration
  
  dead_forward_vector <- numeric(n_samples)
  dead_backward_vector <- numeric(n_samples)
  num_forward_steps_vector <- numeric(n_samples)
  final_sample_index_vector <- numeric(n_samples)
  max_orbit_energy_error_vector <- numeric(n_samples)
  max_c_vector <- numeric(n_samples)
  min_c_vector <- numeric(n_samples)
  # h_vector <- numeric(n_updates)
  # delta_vector <- numeric(n_updates)
  h_vector <- numeric(n_samples)
  delta_vector <- numeric(n_samples)
  
  current_theta <- theta
  current_grad_log_target <- grad_log_target(theta)
  n_evals_ode <- 1
  current_delta <- delta
  current_h <- h
  current_m <- m
  
  hamiltonian <- function(theta, rho) {
    -log_target(theta) + 0.5 * sum(rho ^ 2)
  }
  
  gamma_obs_vector <- c(NULL)
  # kappa_vector <- numeric(n_samples_between_updates)
  kappa_vector <- numeric(n_samples)
  n_updates <- 1
  count_samples_between_updates <- 0
  
  for (i in 1:n_samples) {
    if (i %% 1000 == 0) print(i)
    current_rho <- rnorm(d)
    initial_rho_matrix[i, ] <- current_rho
    initial_theta_matrix[i, ] <- current_theta
    new_state <- single_iteration_standard_biased_progressive_HMC( # obtain one single iteration using fixed length Hamiltonian with adaptive step size
      micro_fun = micro_fun, 
      h = current_h, 
      m = current_m, 
      delta = current_delta, 
      max_c = max_c, 
      log_target = log_target_fun, 
      grad_log_target = grad_log_target, 
      hamiltonian = hamiltonian,
      theta = current_theta, 
      rho = current_rho, 
      grad_log_target_initial = current_grad_log_target,
      d = d,
      norm = norm
    )
    n_evals_ode <- n_evals_ode + new_state$n_evals_ode
    samples_matrix[i, ] <- current_theta <- new_state$theta
    dead_forward_vector[i] <- new_state$dead_forward_indicator
    dead_backward_vector[i] <- new_state$dead_backward_indicator
    num_forward_steps_vector[i] <- new_state$b
    final_sample_index_vector[i] <- new_state$final_index
    max_orbit_energy_error_vector[i] <- new_state$max_orbit_energy_error
    max_c_vector[i] <- new_state$max_c
    min_c_vector[i] <- new_state$min_c
    current_grad_log_target <- new_state$grad_log_target
    
    # count_before_update <- i %% n_samples_between_updates
    gamma_obs_vector <- c(gamma_obs_vector, new_state$gamma_obs_vector)
    # kappa_vector[ifelse(count_before_update == 0, n_samples_between_updates, count_before_update)] <-
    #   new_state$max_orbit_energy_error / current_delta
    # if (count_before_update == 0) {
    #   delta_vector[n_updates] <- (max_orbit_energy_error_threshold / quantile(kappa_vector, probs = prob_max_orbit_energy_error_below_threshold))
    #   if (n_updates >= moving_avg_parameter & moving_avg_indicator) {
    #     current_delta <- mean(tail(delta_vector[delta_vector != 0], moving_avg_parameter))
    #     h_vector[n_updates] <- (current_delta / quantile(gamma_obs_vector, probs = prob_no_step_size_halving)) ^ (1 / 3)
    #     current_h <- mean(tail(h_vector[h_vector != 0], moving_avg_parameter))
    #   } else {
    #     current_delta <- delta_vector[n_updates]
    #     h_vector[n_updates] <- (current_delta / quantile(gamma_obs_vector, probs = prob_no_step_size_halving)) ^ (1 / 3)
    #     current_h <- h_vector[n_updates]
    #   }
    #   current_m <- ceiling(log(path_length / current_h, base = 2))
    #   print(paste0("i: ", i))
    #   print(paste0("h: ", current_h))
    #   print(paste0("delta: ", current_delta))
    #   print(paste0("m: ", current_m))
    #   gamma_obs_vector <- c(NULL)
    #   kappa_vector <- numeric(n_samples_between_updates)
    #   n_updates <- n_updates + 1
    # }
    kappa_vector[i] <-
      new_state$max_orbit_energy_error / current_delta
    
    if (i >= n_samples_before_adaptive) {
      current_delta <- (max_orbit_energy_error_threshold / quantile(kappa_vector[1:i], probs = prob_max_orbit_energy_error_below_threshold))
      # current_log_h <- 1 / 3 * log(current_delta / quantile(gamma_obs_vector[1:i], probs = prob_no_step_size_halving))
      # current_h <- exp(current_log_h)
      current_h <- (current_delta / quantile(gamma_obs_vector, probs = prob_no_step_size_halving)) ^ (1 / 3)
      current_m <- ceiling(log(path_length / current_h, base = 2))
      if (i %% 1000 == 0) {
        print(paste0("i: ", i))
        print(paste0("h: ", current_h))
        print(paste0("delta: ", current_delta))
        print(paste0("m: ", current_m))
      }
    }
  }
  
  list(samples_matrix = samples_matrix, initial_theta_matrix = initial_theta_matrix, initial_rho_matrix = initial_rho_matrix, dead_forward_vector = dead_forward_vector, dead_backward_vector = dead_backward_vector, 
       num_forward_steps_vector = num_forward_steps_vector, final_sample_index_vector = final_sample_index_vector,
       max_orbit_energy_error_vector = max_orbit_energy_error_vector, max_c_vector = max_c_vector, min_c_vector = min_c_vector, 
       final_h = current_h, final_delta = current_delta, final_m = current_m, kappa_vector = kappa_vector, gamma_obs_vector = gamma_obs_vector, 
       h_vector = h_vector, delta_vector = delta_vector,
       n_evals_ode = n_evals_ode)
}

adaptive_step_size_biased_progressive_HMC_sampling <- function(
    micro_fun,
    n_samples, 
    h, 
    m, 
    delta, 
    max_c, 
    log_target, 
    grad_log_target, 
    theta,
    norm = "L-infinity"
) {
  L <- 2 ^ m
  d <- length(theta) # dimension of the parameter vector
  samples_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the actual samples
  samples_matrix[1, ] <- theta
  
  initial_theta_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial theta at each iteration
  initial_rho_matrix <- matrix(nrow = n_samples, ncol = d) # matrix to store the initial rho/momentum at each iteration
  
  dead_forward_vector <- numeric(n_samples)
  dead_backward_vector <- numeric(n_samples)
  num_forward_steps_vector <- numeric(n_samples)
  final_sample_index_vector <- numeric(n_samples)
  max_orbit_energy_error_vector <- numeric(n_samples)
  max_c_vector <- numeric(n_samples)
  min_c_vector <- numeric(n_samples)
  c_matrix <- matrix(nrow = n_samples, ncol = L)
  B_matrix <- matrix(nrow = n_samples, ncol = m)
  
  current_theta <- theta
  current_grad_log_target <- grad_log_target(theta)
  n_evals_ode <- 1
  
  hamiltonian <- function(theta, rho) {
    (-log_target(theta) + 0.5 * sum(rho ^ 2))[1]
  }
  
  for (i in 1:n_samples) {
    if (i %% 1000 == 0) print(i)
    current_rho <- rnorm(d)
    initial_rho_matrix[i, ] <- current_rho
    initial_theta_matrix[i, ] <- current_theta
    new_state <- single_iteration_standard_biased_progressive_HMC( # obtain one single iteration using fixed length Hamiltonian with adaptive step size
      micro_fun = micro_fun, 
      h = h, 
      m = m, 
      delta = delta, 
      max_c = max_c, 
      log_target = log_target_fun, 
      grad_log_target = grad_log_target, 
      hamiltonian = hamiltonian,
      theta = current_theta, 
      rho = current_rho, 
      grad_log_target_initial = current_grad_log_target,
      d = d,
      norm = norm
    )
    n_evals_ode <- n_evals_ode + new_state$n_evals_ode
    samples_matrix[i, ] <- current_theta <- new_state$theta
    dead_forward_vector[i] <- new_state$dead_forward_indicator
    dead_backward_vector[i] <- new_state$dead_backward_indicator
    num_forward_steps_vector[i] <- new_state$b
    final_sample_index_vector[i] <- new_state$final_index
    max_orbit_energy_error_vector[i] <- new_state$max_orbit_energy_error
    max_c_vector[i] <- new_state$max_c
    min_c_vector[i] <- new_state$min_c
    c_matrix[i, ] <- new_state$c_vector
    B_matrix[i, ] <- new_state$B_vector
    
    current_grad_log_target <- new_state$grad_log_target
  }
  
  list(samples_matrix = samples_matrix, initial_theta_matrix = initial_theta_matrix, initial_rho_matrix = initial_rho_matrix, dead_forward_vector = dead_forward_vector, dead_backward_vector = dead_backward_vector, 
       num_forward_steps_vector = num_forward_steps_vector, final_sample_index_vector = final_sample_index_vector,
       max_orbit_energy_error_vector = max_orbit_energy_error_vector, max_c_vector = max_c_vector, min_c_vector = min_c_vector, c_matrix = c_matrix, B_matrix = B_matrix, 
       n_evals_ode = n_evals_ode)
}

