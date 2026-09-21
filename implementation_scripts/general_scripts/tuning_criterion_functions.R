source("implementation_scripts/general_scripts/leapfrog_functions.R")

# Hamiltonian criterion

micro_fun_hamiltonian_error <- function(h, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial, ...) {
  theta_initial <- theta
  rho_initial <- rho
  hamiltonian_initial <- hamiltonian(theta_initial, rho_initial)
  max_c_indicator <- 1 # indicator showing if max_c is reached without achieving the desired criterion
  n_evals_ode <- 0
  for (c in 0:max_c) {
    # print("##########")
    # print(paste0("c: ", c))
    l <- 2 ^ c
    # print(paste0("l: ", l))
    eps <- h / l # micro step size
    current_state <- several_leapfrog_integrator_steps_identity_matrix( # based on the micro step size above, run now l leapfrog steps with this to get one full macro step
      grad_log_target, 
      eps,
      l,
      theta_initial,
      rho_initial,
      grad_log_target_initial
    )
    n_evals_ode <- n_evals_ode + l
    current_hamiltonian <- hamiltonian(current_state$theta_prime, current_state$rho_prime) # check if the number of micro steps covering one macro step gives the desired Hamiltonian error
    # print(paste0("Absolute difference in Hamiltonian between initial and final state: ", abs(hamiltonian_initial - current_hamiltonian)))
    obs_criterion <- abs(hamiltonian_initial - current_hamiltonian)
    obs_criterion <- ifelse(is.na(obs_criterion), Inf, obs_criterion)
    # print(c)
    # print(obs_criterion)
    if (c == 0) gamma_obs <- obs_criterion / h ^ 3
    if (obs_criterion < delta) {
      max_c_indicator <- 0
      break
    }
  }
  list(l = l, c = c, max_c_indicator = max_c_indicator, theta = current_state$theta_prime, rho = current_state$rho_prime, grad_log_target = current_state$grad_log_target_prime,
       n_evals_ode = n_evals_ode,
       gamma_obs = gamma_obs
       # gamma_obs = abs(hamiltonian_initial - current_hamiltonian) * (l ^ 2) / (h ^ 3)
       # gamma_obs = abs(hamiltonian_initial - current_hamiltonian) / (h ^ 3)
  )
}

# Hamiltonian criterion - original

micro_fun_og_hamiltonian_error <- function(h, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial, ...) {
  
  theta_initial <- theta
  rho_initial <- rho
  hamiltonian_initial <- hamiltonian(theta_initial, rho_initial)
  max_c_indicator <- 1 # indicator showing if max_c is reached without achieving the desired criterion
  n_evals_ode <- 0
  
  for (c in 0:max_c) {
    
    # print("##########")
    # print(paste0("c: ", c))
    l <- 2 ^ c
    # print(paste0("l: ", l))
    eps <- h / l # micro step size
    current_rho <- rho 
    current_theta <- theta
    current_grad_log_target <- grad_log_target_initial
    
    hamiltonian_vec <- numeric(l + 1) 
    hamiltonian_vec[1] <- hamiltonian_initial
    
    for (j in 1:l) {
      new_state <- single_leapfrog_integrator_step_identity_matrix(grad_log_target, eps, current_theta, current_rho, current_grad_log_target) # one single standard leapfrog micro step
      current_theta <- new_state$theta_prime
      current_rho <- new_state$rho_prime
      current_grad_log_target <- new_state$grad_log_target_prime
      current_hamiltonian <- hamiltonian(current_theta, current_rho)
      hamiltonian_vec[j + 1] <- current_hamiltonian
    }
    n_evals_ode <- n_evals_ode + l
    max_hamiltonian <- max(hamiltonian_vec)
    min_hamiltonian <- min(hamiltonian_vec)
    obs_criterion <- max_hamiltonian - min_hamiltonian
    obs_criterion <- ifelse(is.na(obs_criterion), Inf, obs_criterion)
    # print(c)
    # print(obs_criterion)
    # print(hamiltonian_vec)
    if (c == 0) gamma_obs <- obs_criterion / h ^ 3
    if (obs_criterion < delta) {
      max_c_indicator <- 0
      break
    }
  }
  list(l = l, c = c, max_c_indicator = max_c_indicator, theta = current_theta, rho = current_rho, grad_log_target = current_grad_log_target,
       n_evals_ode = n_evals_ode,
       gamma_obs = gamma_obs
       # gamma_obs = abs(hamiltonian_initial - current_hamiltonian) * (l ^ 2) / (h ^ 3)
       # gamma_obs = abs(hamiltonian_initial - current_hamiltonian) / (h ^ 3)
  )
}

# Runge-Kutta-Nyström flow (RKNF) criterion

micro_fun_rknf_error <- function(h, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial, norm = "L-infinity", ...) {
  
  if (!(norm %in% c("L1", "L2", "L-infinity"))) {
    stop("Norm must be either: L1, L2 or L-infinity!")
  } else {
    
    if (norm == "L1") {
      norm_fun <- function(x) {
        sum(abs(x))
      }
    } else if (norm == "L2") {
      norm_fun <- function(x) {
        sqrt(sum(x ^ 2))
      }
    } else {
      norm_fun <- function(x) {
        max(abs(x))
      }
    }
    
  }
  
  theta_initial <- theta
  rho_initial <- rho
  if (is.null(grad_log_target_initial)) {
    grad_log_target_initial <- grad_log_target(theta)
  }
  max_c_indicator <- 1 # indicator showing if max_c is reached without achieving the desired criterion
  n_evals_ode <- 0
  for (c in 0:max_c) {
    # print("##########")
    # print(paste0("c: ", c))
    l <- 2 ^ c
    # print(paste0("l: ", l))
    eps <- h / l # micro step size
    current_rho <- rho 
    current_theta <- theta
    current_grad_log_target <- grad_log_target_initial
    criterion_satisfied_indicator <- 1 # indicator that shows if max(flow error forward, flow error backward) < delta for all micro steps or not
    
    for (j in 1:l) { # run l micro steps
      new_state <- single_leapfrog_integrator_step_identity_matrix(grad_log_target, eps, current_theta, current_rho, current_grad_log_target) # one single standard leapfrog micro step
      theta_intermediate <- 1 / 2 * (current_theta + new_state$theta_prime) + # theta_m
        eps / 8 * (current_rho - new_state$rho_prime)
      grad_log_target_theta_intermediate <- grad_log_target(theta_intermediate)
      n_evals_ode <- n_evals_ode + 2
      
      forward_theta_double_prime <- current_theta + eps * current_rho + # third order embedded based on already computed leapfrog from above
        (eps ^ 2) * (1 / 3 * current_grad_log_target +  1 / 6 * new_state$grad_log_target_prime)
      forward_rho_double_prime <- current_rho + 
        eps * (1 / 6 * current_grad_log_target + 2 / 3 * grad_log_target_theta_intermediate + 1 / 6 * new_state$grad_log_target_prime)
      forward_error <- norm_fun(c(new_state$theta_prime, new_state$rho_prime) - c(forward_theta_double_prime, forward_rho_double_prime))
      
      backward_theta_double_prime <- new_state$theta_prime + eps * (-new_state$rho_prime) + # similar as above, but now based on the new state obtained from the single micro leapfrog step
        (eps ^ 2) * (1 / 3 * new_state$grad_log_target_prime + 1 / 6 * current_grad_log_target)
      backward_rho_double_prime <- -new_state$rho_prime + 
        eps * (1 / 6 * new_state$grad_log_target_prime + 2 / 3 * grad_log_target_theta_intermediate + 1 / 6 * current_grad_log_target)
      backward_error <- norm_fun(c(current_theta, -current_rho) - c(backward_theta_double_prime, backward_rho_double_prime))
      obs_criterion <- max(forward_error, backward_error)
      obs_criterion <- ifelse(is.na(obs_criterion), Inf, obs_criterion)
      if (c == 0) gamma_obs <- obs_criterion / h ^ 3
      if (obs_criterion > delta) { # flow error condition
        criterion_satisfied_indicator <- 0 # set the indicator to zero if the condition above does not hold for a single micro leapfrog step --> need to try a different value of l right away
        break
      } else { # prepare the next initial state for the next micro leapfrog step
        current_theta <- new_state$theta_prime
        current_rho <- new_state$rho_prime
        current_grad_log_target <- new_state$grad_log_target_prime
      }
      
    }
    
    if (criterion_satisfied_indicator == 1) {
      max_c_indicator <- 0
      break # if for a given value of l, all the l micro leapfrog steps satisfy the condition --> the given l (or c where l = 2 ^ c) is ok, no need to run further 
    }
    
  }
  list(l = l, c = c, max_c_indicator = max_c_indicator, theta = current_theta, rho = current_rho, grad_log_target = current_grad_log_target,
       n_evals_ode = n_evals_ode,
       gamma_obs = gamma_obs)
}

# Multiresolution flow (MRF) criterion

micro_fun_mrf_error <- function(h, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial, norm = "L-infinity",
                                theta_final_backward = NULL, rho_final_backward = NULL, grad_log_target_final_backward = NULL, forward_micro_run_indicator = TRUE, ...) {
  
  # NOTE: c = 0 here corresponds to c = 1 in the text!
  # The general sampler implementation does not require any additional modifications
  # if l is returned instead of l_finer (due to the condition forward_run_micro$l == 1 in the samplers), 
  # i.e. this corresponds to starting from c = 0 and returning c and l.
  # But this implies that 2^(c+1) is the correct number of micro steps taken for the returning state. 
  # In other words, to get the correct value of c defined in the text, one needs to add by 1 here, i.e. c + 1. 
  
  if (!(norm %in% c("L1", "L2", "L-infinity"))) {
    stop("Norm must be either: L1, L2 or L-infinity!")
  } else {
    
    if (norm == "L1") {
      norm_fun <- function(x) {
        sum(abs(x))
      }
    } else if (norm == "L2") {
      norm_fun <- function(x) {
        sqrt(sum(x ^ 2))
      }
    } else {
      norm_fun <- function(x) {
        max(abs(x))
      }
    }
    
  }
  
  theta_initial <- theta
  rho_initial <- rho
  if (is.null(grad_log_target_initial)) {
    grad_log_target_initial <- grad_log_target(theta)
  }
  max_c_indicator <- 1 # indicator showing if max_c is reached without achieving the desired criterion
  n_evals_ode <- 0
  for (c in 0:max_c) {
    # print("##########")
    # print(paste0("c: ", c))
    l <- 2 ^ c
    # print(paste0("l: ", l))
    eps <- h / l # micro step size
    l_finer <- 2 ^ (c + 1)
    eps_finer <- h / l_finer
    
    if (c == 0) {
      new_coarse_state <- several_leapfrog_integrator_steps_identity_matrix( # based on the micro step size above, run now l leapfrog steps with this to get one full macro step
        grad_log_target, 
        eps,
        l,
        theta_initial,
        rho_initial,
        grad_log_target_initial
      )
      n_evals_ode <- n_evals_ode + l
      theta_coarse <- new_coarse_state$theta_prime
      rho_coarse <- new_coarse_state$rho_prime
      grad_log_target_coarse <- new_coarse_state$grad_log_target_prime
    }
    
    if (!forward_micro_run_indicator & c == max_c) {
      theta_finer <- theta_final_backward
      rho_finer <- rho_final_backward
      grad_log_target_finer <- grad_log_target_final_backward
    } else {
      new_finer_state <- several_leapfrog_integrator_steps_identity_matrix( # based on the finer micro step size above, run now l_finer leapfrog steps with this to get one full macro step
        grad_log_target, 
        eps_finer,
        l_finer,
        theta_initial,
        rho_initial,
        grad_log_target_initial
      )
      n_evals_ode <- n_evals_ode + l_finer
      theta_finer <- new_finer_state$theta_prime
      rho_finer <- new_finer_state$rho_prime
      grad_log_target_finer <- new_finer_state$grad_log_target_prime  
    }
    
    forward_error <- norm_fun(c(theta_coarse, rho_coarse) - c(theta_finer, rho_finer)) # check the difference between the two obtained states
    
    backward_state_finer <- several_leapfrog_integrator_steps_identity_matrix( # based on the micro step size above, run now l leapfrog steps with this to get one full macro step backward from the forward state obtained from the finer micro step size
      grad_log_target, 
      eps,
      l,
      theta_finer,
      -rho_finer, # backward
      grad_log_target_finer
    )
    n_evals_ode <- n_evals_ode + l
    theta_finer_backward <- backward_state_finer$theta_prime
    rho_finer_backward <- backward_state_finer$rho_prime
    grad_log_target_finer_backward <- backward_state_finer$grad_log_target_prime
    
    backward_error <- norm_fun(c(theta_finer_backward, rho_finer_backward) - c(theta_initial, -rho_initial)) # negative here because backward
    obs_criterion <- max(backward_error, forward_error)
    obs_criterion <- ifelse(is.na(obs_criterion), Inf, obs_criterion)
    if (c == 0) gamma_obs <- obs_criterion / h ^ 3
    if (obs_criterion < delta) {
      max_c_indicator <- 0
      break
    } else {
      theta_coarse <- theta_finer
      rho_coarse <- rho_finer
      grad_log_target_coarse <- grad_log_target_finer
    }
    
  }
  
  # list(l = l, c = c, max_c_indicator = max_c_indicator, theta = theta_coarse, rho = rho_coarse, grad_log_target = grad_log_target_coarse)
  # list(l = l_finer, c = c + 1, max_c_indicator = max_c_indicator, theta = theta_finer, rho = rho_finer, grad_log_target = grad_log_target_finer, n_evals_ode = n_evals_ode)
  list(l = l, c = c, max_c_indicator = max_c_indicator, theta = theta_finer, rho = rho_finer, grad_log_target = grad_log_target_finer, 
       theta_final_backward = theta_finer_backward, rho_final_backward = rho_finer_backward, grad_log_target_final_backward = grad_log_target_finer_backward,
       n_evals_ode = n_evals_ode,
       gamma_obs = gamma_obs)
  
}


# micro_fun_mrf_error <- function(h, delta, max_c, log_target, grad_log_target, hamiltonian, theta, rho, grad_log_target_initial, norm = "L-infinity") {
#   
#   if (!(norm %in% c("L1", "L2", "L-infinity"))) {
#     stop("Norm must be either: L1, L2 or L-infinity!")
#   } else {
#     
#     if (norm == "L1") {
#       norm_fun <- function(x) {
#         sum(abs(x))
#       }
#     } else if (norm == "L2") {
#       norm_fun <- function(x) {
#         sqrt(sum(x ^ 2))
#       }
#     } else {
#       norm_fun <- function(x) {
#         max(abs(x))
#       }
#     }
#     
#   }
#   
#   theta_initial <- theta
#   rho_initial <- rho
#   if (is.null(grad_log_target_initial)) {
#     grad_log_target_initial <- grad_log_target(theta)
#   }
#   max_c_indicator <- 1 # indicator showing if max_c is reached without achieving the desired criterion
#   n_evals_ode <- 0
#   for (c in 0:max_c) {
#     # print("##########")
#     # print(paste0("c: ", c))
#     l <- 2 ^ c
#     # print(paste0("l: ", l))
#     eps <- h / l # micro step size
#     l_finer <- 2 ^ (c + 1)
#     eps_finer <- h / l_finer
#     
#     if (c == 0) {
#       new_coarse_state <- several_leapfrog_integrator_steps_identity_matrix( # based on the micro step size above, run now l leapfrog steps with this to get one full macro step
#         grad_log_target, 
#         eps,
#         l,
#         theta_initial,
#         rho_initial,
#         grad_log_target_initial
#       )
#       n_evals_ode <- n_evals_ode + l
#       theta_coarse <- new_coarse_state$theta_prime
#       rho_coarse <- new_coarse_state$rho_prime
#       grad_log_target_coarse <- new_coarse_state$grad_log_target_prime
#     }
#     
#     new_finer_state <- several_leapfrog_integrator_steps_identity_matrix( # based on the finer micro step size above, run now l_finer leapfrog steps with this to get one full macro step
#       grad_log_target, 
#       eps_finer,
#       l_finer,
#       theta_initial,
#       rho_initial,
#       grad_log_target_initial
#     )
#     n_evals_ode <- n_evals_ode + l_finer
#     theta_finer <- new_finer_state$theta_prime
#     rho_finer <- new_finer_state$rho_prime
#     grad_log_target_finer <- new_finer_state$grad_log_target_prime
#     
#     forward_error <- norm_fun(c(theta_coarse, rho_coarse) - c(theta_finer, rho_finer)) # check the difference between the two obtained states
#     
#     backward_state_finer <- several_leapfrog_integrator_steps_identity_matrix( # based on the micro step size above, run now l leapfrog steps with this to get one full macro step backward from the forward state obtained from the finer micro step size
#       grad_log_target, 
#       eps,
#       l,
#       theta_finer,
#       -rho_finer, # backward
#       grad_log_target_finer
#     )
#     n_evals_ode <- n_evals_ode + l
#     theta_finer_backward <- backward_state_finer$theta_prime
#     rho_finer_backward <- backward_state_finer$rho_prime
#     
#     backward_error <- norm_fun(c(theta_finer_backward, rho_finer_backward) - c(theta_initial, -rho_initial)) # negative here because backward
#     obs_criterion <- max(backward_error, forward_error)
#     if (c == 0) gamma_obs <- obs_criterion / h ^ 3
#     if (obs_criterion < delta) {
#       max_c_indicator <- 0
#       break
#     } else {
#       theta_coarse <- theta_finer
#       rho_coarse <- rho_finer
#       grad_log_target_coarse <- grad_log_target_finer
#     }
#     
#   }
#   
#   # list(l = l, c = c, max_c_indicator = max_c_indicator, theta = theta_coarse, rho = rho_coarse, grad_log_target = grad_log_target_coarse)
#   # list(l = l_finer, c = c + 1, max_c_indicator = max_c_indicator, theta = theta_finer, rho = rho_finer, grad_log_target = grad_log_target_finer, n_evals_ode = n_evals_ode)
#   list(l = l, c = c, max_c_indicator = max_c_indicator, theta = theta_finer, rho = rho_finer, grad_log_target = grad_log_target_finer, 
#        n_evals_ode = n_evals_ode,
#        gamma_obs = gamma_obs)
#   
# }
