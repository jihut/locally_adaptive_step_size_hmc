source("implementation_scripts/general_scripts/tuning_criterion_functions.R")

# One iteration of multinomial HMC based on L number of macro steps with a macro step size of h - concurrent format
single_iteration_multinomial_HMC <- function(
    micro_fun, h, L, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial = NULL, d, norm = "L-infinity"
) {
  
  # d = dimension of theta vector
  
  b <- sample(0:L, 1, replace = T) # sample uniformly the rightmost index of the orbit
  a <- b - L # given L and b --> get the leftmost index of the orbit
  initial_index <- L + 1 - b # this will be the row index of the initial sample in the orbit
  # print(paste0("b:", b))
  # print(paste0("a:", a))
  # print(paste0("initial_index:", initial_index))
  
  # i.e. need to run b macro steps of step size h forward from (theta, rho) and abs(a) macro steps of size h backwards from (theta, rho)
  
  final_index <- 0
  final_theta <- theta
  final_rho <- rho
  # hamiltonian_vector <- numeric(L + 1)
  hamiltonian_vector <- rep(NA, L + 1)
  c_vector <- rep(-1, L)
  gamma_obs_vector <- rep(-1, L)
  
  dead_forward_indicator <- 0 # if one of the weights when doing forward steps becomes zero because the minimum number of micro steps forward and backward to achieve same criterion does not coincide
  dead_backward_indicator <- 0 # similar as above, but for backward
  
  current_rho <- rho
  current_theta <- theta
  if (!is.null(grad_log_target_initial)) {
    current_grad_log_target <- grad_log_target_initial
    n_evals_ode <- 0
  } else {
    current_grad_log_target <- grad_log_target(theta)
    n_evals_ode <- 1
  }
  final_grad_log_target <- current_grad_log_target
  current_hamiltonian <- hamiltonian(current_theta, current_rho)
  # current_weight <- exp(-current_hamiltonian)
  current_weight <- 1 # scale using the initial value to use "log sum exp" trick
  sum_weights <- current_weight
  log_current_weight <- 0
  hamiltonian_vector[1] <- current_hamiltonian
  hamiltonian_vector_index <- 2
  # First: Forward macro steps
  if (b > 0) {
    for (i in 1:b) {
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
      gamma_obs_vector[hamiltonian_vector_index - 1] <- forward_run_micro$gamma_obs
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
        new_theta <- forward_run_micro$theta
        new_rho <- forward_run_micro$rho
        new_hamiltonian <- hamiltonian(new_theta, new_rho)
        hamiltonian_vector[hamiltonian_vector_index] <- new_hamiltonian
        c_vector[hamiltonian_vector_index - 1] <- forward_run_micro$c
        hamiltonian_vector_index <- hamiltonian_vector_index + 1
        log_new_weight <- (current_hamiltonian - new_hamiltonian + log_current_weight) # see Section 3.3 and B.2 in WALNUTS
        new_weight <- exp(log_new_weight)
        sum_weights <- sum_weights + new_weight
        if (runif(1) < new_weight / sum_weights) {
          final_index <- i
          final_theta <- new_theta
          final_rho <- new_rho
          final_grad_log_target <- forward_run_micro$grad_log_target
        }
        current_theta <- new_theta
        current_rho <- new_rho
        current_hamiltonian <- new_hamiltonian
        log_current_weight <- log_new_weight
        current_grad_log_target <- forward_run_micro$grad_log_target
      }
      
    }
  }
  
  
  # Next: Prepare for backward macro steps
  
  current_rho <- -rho # backwards now, need to flip the initial momentum
  current_theta <- theta 
  if (!is.null(grad_log_target_initial)) {
    current_grad_log_target <- grad_log_target_initial
  } else {
    current_grad_log_target <- grad_log_target(theta)
    n_evals_ode <- n_evals_ode + 1
  }
  current_hamiltonian <- hamiltonian(current_theta, current_rho)
  # current_weight <- exp(-current_hamiltonian)
  current_weight <- 1 # scale using the initial value to use log sum exp trick
  log_current_weight <- 0
  
  if (a < 0) {
    for (i in 1:abs(a)) {
      
      # For each macro step of size h, need to find the number of micro steps to achieve the desired error criterion
      forward_run_micro <- micro_fun( # NOTE: Forward with initial momentum flipped, so technically backward in the global sense!
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
      gamma_obs_vector[hamiltonian_vector_index - 1] <- forward_run_micro$gamma_obs
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
        new_theta <- forward_run_micro$theta
        new_rho <- forward_run_micro$rho
        new_hamiltonian <- hamiltonian(new_theta, new_rho)
        hamiltonian_vector[hamiltonian_vector_index] <- new_hamiltonian
        c_vector[hamiltonian_vector_index - 1] <- forward_run_micro$c
        hamiltonian_vector_index <- hamiltonian_vector_index + 1
        log_new_weight <- (current_hamiltonian - new_hamiltonian + log_current_weight) # see Section 3.3 and B.2 in WALNUTS
        new_weight <- exp(log_new_weight)
        sum_weights <- sum_weights + new_weight
        if (runif(1) < new_weight / sum_weights) {
          final_index <- -i
          final_theta <- new_theta
          final_rho <- new_rho
          final_grad_log_target <- forward_run_micro$grad_log_target
        }
        current_theta <- new_theta
        current_rho <- new_rho
        current_hamiltonian <- new_hamiltonian
        log_current_weight <- log_new_weight
        current_grad_log_target <- forward_run_micro$grad_log_target
      }
      
    }
  }
  
  list(
    b = b, 
    final_index = final_index, 
    theta = final_theta,
    rho = final_rho,
    grad_log_target = final_grad_log_target,
    dead_forward_indicator = dead_forward_indicator,
    dead_backward_indicator = dead_backward_indicator,
    hamiltonian_vector = hamiltonian_vector,
    # max_orbit_energy_error = max(hamiltonian_vector[hamiltonian_vector != 0]) - min(hamiltonian_vector[hamiltonian_vector != 0]),
    max_orbit_energy_error = max(hamiltonian_vector, na.rm = T) - min(hamiltonian_vector, na.rm = T),
    c_vector = c_vector, 
    max_c = max(c_vector[1:(hamiltonian_vector_index - 2)]),
    min_c = min(c_vector[1:(hamiltonian_vector_index - 2)]),
    n_evals_ode = n_evals_ode,
    # gamma_obs = max(gamma_obs_vector[gamma_obs_vector != -1])
    # gamma_obs = mean(gamma_obs_vector[gamma_obs_vector != -1])
    # gamma_obs = median(gamma_obs_vector[gamma_obs_vector != -1])
    gamma_obs_vector = gamma_obs_vector[gamma_obs_vector != -1]
  )
  
}

# N iterations of multinomial HMC where each iteration is based on L macro steps of step size h with adaptive micro step size

adaptive_step_size_multinomial_HMC_warmup <- function(
    micro_fun,
    n_samples, 
    # n_samples_between_updates, 
    # moving_avg_indicator = FALSE, 
    # moving_avg_parameter = 5, 
    n_samples_before_adaptive = 100, 
    h, 
    L, 
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
  current_L <- L
  
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
    new_state <- single_iteration_multinomial_HMC( # obtain one single iteration using fixed length Hamiltonian with adaptive step size
      micro_fun = micro_fun, 
      h = current_h, 
      L = current_L, 
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
    #   current_L <- ceiling(path_length / current_h)
    #   print(paste0("i: ", i))
    #   print(paste0("h: ", current_h))
    #   print(paste0("delta: ", current_delta))
    #   print(paste0("L: ", current_L))
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
      current_L <- ceiling(path_length / current_h)
      if (i %% 1000 == 0) {
        print(paste0("i: ", i))
        print(paste0("h: ", current_h))
        print(paste0("delta: ", current_delta))
        print(paste0("L: ", current_L))
      }
    }
  }
  
  list(samples_matrix = samples_matrix, initial_theta_matrix = initial_theta_matrix, initial_rho_matrix = initial_rho_matrix, dead_forward_vector = dead_forward_vector, dead_backward_vector = dead_backward_vector, 
       num_forward_steps_vector = num_forward_steps_vector, final_sample_index_vector = final_sample_index_vector,
       max_orbit_energy_error_vector = max_orbit_energy_error_vector, max_c_vector = max_c_vector, min_c_vector = min_c_vector, 
       final_h = current_h, final_delta = current_delta, final_L = current_L, kappa_vector = kappa_vector, gamma_obs_vector = gamma_obs_vector, 
       h_vector = h_vector, delta_vector = delta_vector,
       n_evals_ode = n_evals_ode)
}

adaptive_step_size_multinomial_HMC_sampling <- function(
    micro_fun,
    n_samples, 
    h, 
    L, 
    delta, 
    max_c, 
    log_target, 
    grad_log_target, 
    theta,
    norm = "L-infinity"
) {
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
  
  current_theta <- theta
  current_grad_log_target <- grad_log_target(theta)
  n_evals_ode <- 1
  
  hamiltonian <- function(theta, rho) {
    -log_target(theta) + 0.5 * sum(rho ^ 2)
  }
  
  for (i in 1:n_samples) {
    if (i %% 1000 == 0) print(i)
    current_rho <- rnorm(d)
    initial_rho_matrix[i, ] <- current_rho
    initial_theta_matrix[i, ] <- current_theta
    new_state <- single_iteration_multinomial_HMC( # obtain one single iteration using fixed length Hamiltonian with adaptive step size
      micro_fun = micro_fun, 
      h = h, 
      L = L, 
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
    current_grad_log_target <- new_state$grad_log_target
  }
  
  list(samples_matrix = samples_matrix, initial_theta_matrix = initial_theta_matrix, initial_rho_matrix = initial_rho_matrix, dead_forward_vector = dead_forward_vector, dead_backward_vector = dead_backward_vector, 
       num_forward_steps_vector = num_forward_steps_vector, final_sample_index_vector = final_sample_index_vector,
       max_orbit_energy_error_vector = max_orbit_energy_error_vector, max_c_vector = max_c_vector, min_c_vector = min_c_vector, c_matrix = c_matrix, 
       n_evals_ode = n_evals_ode)
}
