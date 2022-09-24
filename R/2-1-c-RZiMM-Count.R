###########################################################
### A Zero-Inflated Negative Binomial Model

#' A Zero-Inflated Negative Binomial Model
#'   - for RZiMM-Count Model
#'
#' @param y A vector of expression levels
#' @param g_ind A vector of cell-type clusters
#' @param m_ind Sample index, for adjusting within-sample correlation; if set to be NULL,
#'   the within-sample correlation will not be adjusted. 
#' @param g_list A list for cell-type cluster indices
#' @param m_list A list for batch indices
#' @param gm_list A two-way list for cell-type cluster and batch indices
#' @param n_iter Maximum number of total iteration. 
#' @param lambda1 Weight of regularization for group effects. 
#' @param lambda2 Weight of regularization for sample effects. 
#' @param epsilon Iteration stops when the change from last iteration is smaller than this value.
#' @param epsilon0 A small number for calculating the number of parameters. 
#' @param output Whether export a vector or a list
#' @param plot Whether plotting the trajectories for the loss and log-likelihood
#' @param max_l_alpha Upper threshold for l
#' @param min_l_alpha Lower threshold for l
#' @param max_est Upper threshold for estimate
#' @param min_est Lower threshold for estimate
#' @param max_logit_pi Upper threshold for pi
#' @param min_logit_pi Lower threshold for pi
#'
#' @return Returns a vector of parameters.
#'
#' @importFrom methods new
#' @importFrom stats rnorm
#' @importFrom magic adiag
#'
#' @export
zinb_regulated <- function(y, g_ind, m_ind = NULL, 
                           g_list = NULL, m_list = NULL, gm_list = NULL, 
                           lambda1 = 1, lambda2 = 0, n_iter = 10, epsilon = 10**-6, 
                           output = "vec", epsilon0 = 10**-6, plot = FALSE, 
                           max_l_alpha = 10, min_l_alpha = -10, max_est = 100, min_est = -100, 
                           max_logit_pi = 100, min_logit_pi = -100) {
  
  if(is.null(g_list)) {
    g_list <- list()
    gm_list <- list()
    m_list <- list()
    for(iii in 1:length(unique(g_ind))) {
      g_list[[iii]] <- which(g_ind == iii)
      gm_list[[iii]] <- list()
      if(!is.null(m_ind)) {
        for(jjj in 1:(length(unique(m_ind)) - 1)) {
          if(iii == 1) 
            m_list[[jjj]] <- which(m_ind == jjj)
          gm_list[[iii]][[jjj]] <- which(g_ind == iii & m_ind == jjj)
        }
      }
    }
  }
  
  n_ <- length(y)
  n_group <- length(g_list)
  
  pi_est <- max(min(mean(y != 0), 1 - 10**-6), 10**-6)
  logit_pi_est <- log(pi_est/(1 - pi_est))
  alpha_est <- 1
  log_alpha_est <- 0
  m_est <- sapply(g_list, function(x, y) mean(log(y[x][y[x] > 0])), y = y)
  m_est[is.na(m_est)] <- 0
  if(!is.null(m_ind)) {
    e_est <- sapply(m_list, function(x, y) mean(log(y[x][y[x] > 0])), y = y) - mean(m_est)
    e_est[is.na(e_est)] <- 0
  }
  z <- (y != 0)*1
  
  if(!is.null(m_ind)) {
    n_par <- 2 + length(m_est) + length(e_est)
    est_ <- c(logit_pi_est, log_alpha_est, m_est, e_est)
    e_est <- c(e_est, 0)
    lambda_est <- exp(m_est[g_ind] + e_est[m_ind])
  } else {
    n_par <- 2 + length(m_est)
    est_ <- c(logit_pi_est, log_alpha_est, m_est)
    lambda_est <- exp(m_est)[g_ind]
  }
  
  first_deriv <- numeric(n_par)
  second_deriv <- matrix(0, n_par, n_par)
  log_lik <- loss <- numeric(n_iter + 1)
  log_lik[1] <- -log(1 + exp(logit_pi_est))*n_ + 
    sum((1 - z)*log(1 + exp(logit_pi_est - 1/alpha_est*log(1 + alpha_est*lambda_est))) + 
          z*(logit_pi_est + lgamma(y + 1/alpha_est) - lgamma(1/alpha_est) + 
               y*(log(alpha_est) + log(lambda_est)) - (1/alpha_est + y)*log(1 + alpha_est*lambda_est)))
  loss[1] <- -log_lik[1] + 
    lambda1*sum(abs(m_est %*% t(rep(1, n_group)) - rep(1, n_group) %*% t(m_est)))/2
  if(!is.null(m_ind) && lambda2 > 0)
    loss[1] <- loss[1] + lambda2*sum(abs(e_est))
  for(i in 1:n_iter) {
    
    est_0 <- est_
    ttt <- 1 + alpha_est*lambda_est
    aaa <- log(ttt)/alpha_est
    bbb <- lambda_est/ttt
    exp_t_neg_l <- exp(logit_pi_est - aaa)
    expit_exp_t_neg_l <- 1/(1 + exp_t_neg_l)
    expit_exp_t_neg_l <- 1 - expit_exp_t_neg_l
    first_deriv[1] <- -n_*pi_est + sum(z + (1 - z)*expit_exp_t_neg_l)
    first_deriv[2] <- sum(-(1 - z)*expit_exp_t_neg_l*(-aaa + bbb) + 
                            z*(-(digamma(y + 1/alpha_est) - digamma(1/alpha_est))/alpha_est + 
                                 aaa + (y/lambda_est - 1)*bbb))
    first_deriv_l_t_l <- ((y/lambda_est - 1)*z - (1 - z)*expit_exp_t_neg_l)*bbb
    if(!is.null(m_ind)) {
      first_deriv[-(1:2)] <- sapply(append(g_list, m_list), function(x) sum(first_deriv_l_t_l[x]))
    } else {
      first_deriv[-(1:2)] <- sapply(g_list, function(x) sum(first_deriv_l_t_l[x]))
    }
    
    second_deriv[1, 1] <- -pi_est*(1 - pi_est)*n_ + sum((1 - z)*expit_exp_t_neg_l*(1 - expit_exp_t_neg_l))
    second_deriv[2, 2] <- sum(-(1 - z)*expit_exp_t_neg_l*((1 - expit_exp_t_neg_l)*(-aaa + bbb)**2 + 
                                                            expit_exp_t_neg_l*(-bbb + aaa - alpha_est*bbb**2)) + 
                                z*((digamma(y + 1/alpha_est) - digamma(1/alpha_est))/alpha_est + 
                                     (trigamma(y + 1/alpha_est) - trigamma(1/alpha_est))/alpha_est**2 - 
                                     aaa + bbb*(1 - (y - lambda_est)/ttt*alpha_est)))
    second_deriv[1, 2] <- second_deriv[2, 1] <- -sum((1 - z)*expit_exp_t_neg_l*(1 - expit_exp_t_neg_l)*(-aaa + bbb))
    second_deriv_pi_me <- -(1 - z)*expit_exp_t_neg_l*(1 - expit_exp_t_neg_l)*bbb
    second_deriv_al_me <- ((1 - z)*expit_exp_t_neg_l*((1 - expit_exp_t_neg_l)*(-aaa + bbb) + alpha_est*bbb) - 
                             z*(y/lambda_est - 1)*alpha_est*bbb)*bbb
    second_deriv_me_me <- -(z*(1 + alpha_est*y) + (1 - z)*expit_exp_t_neg_l*
                              (1 - (1 - expit_exp_t_neg_l)*lambda_est))*bbb**2/lambda_est
    if(!is.null(m_ind)) {
      second_deriv[1, -(1:2)] <- second_deriv[-(1:2), 1] <- 
        sapply(append(g_list, m_list), function(x) sum(second_deriv_pi_me[x]))
      second_deriv[2, -(1:2)] <- second_deriv[-(1:2), 2] <- 
        sapply(append(g_list, m_list), function(x) sum(second_deriv_al_me[x]))
      diag(second_deriv)[-(1:2)] <- sapply(append(g_list, m_list), function(x) sum(second_deriv_me_me[x]))
      second_deriv[(1+length(m_est)+2):n_par, 2+1:length(m_est)] <- 
        sapply(gm_list, function(x) {sapply(x, function(z) sum(second_deriv_me_me[z]))})
      second_deriv[2+1:length(m_est), (1+length(m_est)+2):n_par] <- 
        t(second_deriv[(1+length(m_est)+2):n_par, 2+1:length(m_est)])
    } else {
      second_deriv[1, -(1:2)] <- second_deriv[-(1:2), 1] <- 
        sapply(g_list, function(x) sum(second_deriv_pi_me[x]))
      second_deriv[2, -(1:2)] <- second_deriv[-(1:2), 2] <- 
        sapply(g_list, function(x) sum(second_deriv_al_me[x]))
      diag(second_deriv)[-(1:2)] <- sapply(g_list, function(x) sum(second_deriv_me_me[x]))
    }
    
    if(lambda1 > 0) {
      b_mat <- abs(m_est %*% t(rep(1, n_group)) - rep(1, n_group) %*% t(m_est)) + epsilon0
      delta_mat <- lambda1/n_group/(n_group-1)*(diag(rowSums(1/b_mat)) - 1/b_mat)
      first_deriv[1:length(m_est)+2] <- first_deriv[1:length(m_est)+2] - as.numeric(delta_mat %*% m_est)
      second_deriv[1:length(m_est)+2, 1:length(m_est)+2] <- 
        second_deriv[1:length(m_est)+2, 1:length(m_est)+2] - delta_mat
    }
    
    if(lambda2 > 0 && !is.null(m_ind)){
      c_mat <- abs(e_est[-length(e_est)]) + epsilon0
      first_deriv[(1+length(m_est)+2):n_par] <- 
        first_deriv[(1+length(m_est)+2):n_par] - 1/c_mat*e_est[-length(e_est)]*lambda2
      diag(second_deriv[(1+length(m_est)+2):n_par, (1+length(m_est)+2):n_par]) <- 
        diag(second_deriv[(1+length(m_est)+2):n_par, (1+length(m_est)+2):n_par]) - 1/c_mat*lambda2
    }
    
    second_deriv_inv <- try(solve(second_deriv), silent = TRUE)
    lll <- 0
    while(class(second_deriv_inv)[1] == "try-error" && lll < 20) {
      second_deriv_inv <- try(solve(second_deriv - epsilon0*diag(n_par)*10**lll))
      lll <- lll + 1
    }
    if(class(second_deriv_inv)[1] == "try-error")
      break()
    est_ <- est_ - as.numeric(second_deriv_inv %*% first_deriv)
    est_[1] <- min(max(est_[1], min_logit_pi), max_logit_pi)
    est_[2] <- min(max(est_[2], min_l_alpha), max_l_alpha)
    est_[-c(1, 2)] <- sapply(est_[-c(1, 2)], function(x) min(max(x, min_est), max_est))
    logit_pi_est <- est_[1]
    log_alpha_est <- est_[2] 
    m_est <- est_[1:length(m_est)+2]
    pi_est <- 1/(1 + exp(-logit_pi_est))
    alpha_est <- exp(log_alpha_est)
    if(!is.null(m_ind)){
      e_est <- c(est_[(1+length(m_est)+2):n_par], 0)
      lambda_est <- exp(m_est[g_ind] + e_est[m_ind])
    } else {
      lambda_est <- exp(m_est)[g_ind]
    }
    log_lik[i + 1] <- -log(1 + exp(logit_pi_est))*n_ + 
      sum((1 - z)*log(1 + exp(logit_pi_est - 1/alpha_est*log(1 + alpha_est*lambda_est))) + 
            z*(logit_pi_est + lgamma(y + 1/alpha_est) - lgamma(1/alpha_est) + 
                 y*(log(alpha_est) + log(lambda_est)) - (1/alpha_est + y)*log(1 + alpha_est*lambda_est)))
    loss[i + 1] <- -log_lik[i + 1] +  
      lambda1*sum(abs(m_est %*% t(rep(1, n_group)) - rep(1, n_group) %*% t(m_est)))/2
    if(lambda2 > 0 && !is.null(m_ind))
      loss[i + 1] <- loss[i + 1] + lambda2*sum(abs(e_est))
    
    if(is.nan(loss[i + 1]) || sum(is.nan(est_)) > 0)
      break()
    if((loss[i] - loss[i + 1] < epsilon) || (sum(abs(est_ - est_0)) < epsilon)) {
      if(loss[i] < loss[i + 1])
        est_ <- est_0
      break()
    }
  }
  
  if(plot) {
    plot(log_lik[1:(i+1)])
    plot(loss[1:(i+1)])
  }
  
  if(output == "vec")
    return(c(1/(1+exp(-est_[1])), exp(est_[2]), est_[-(1:2)]))
  return(list(pi = pi_est, alpha = alpha_est, m = m_est, e = ifelse(is.null(m_ind), NA, e_est), 
              loss = loss[1:(i+1)], log_lik = log_lik[1:(i+1)]))
}

###########################################################
### A Zero-Inflated Poisson Model

#' A Zero-Inflated Poisson Model
#'   - for RZiMM-Count Model
#'
#' @param y A vector of expression levels
#' @param g_ind A vector of cell-type clusters
#' @param m_ind Sample index, for adjusting within-sample correlation; if set to be NULL,
#'   the within-sample correlation will not be adjusted. 
#' @param g_list A list for cell-type cluster indices
#' @param m_list A list for batch indices
#' @param gm_list A two-way list for cell-type cluster and batch indices
#' @param n_iter Maximum number of total iteration. 
#' @param lambda1 Weight of regularization for group effects. 
#' @param lambda2 Weight of regularization for sample effects. 
#' @param epsilon Iteration stops when the change from last iteration is smaller than this value.
#' @param epsilon0 A small number for calculating the number of parameters. 
#' @param output Whether export a vector or a list
#' @param plot Whether plotting the trajectories for the loss and log-likelihood
#' @param max_est Upper threshold for estimate
#' @param min_est Lower threshold for estimate
#' @param max_logit_pi Upper threshold for pi
#' @param min_logit_pi Lower threshold for pi
#'
#' @return Returns a vector of parameters.
#'
#' @importFrom methods new
#' @importFrom stats rnorm
#' @importFrom magic adiag
#'
#' @export
zip_regulated <- function(y, g_ind, m_ind = NULL, 
                          g_list = NULL, m_list = NULL, gm_list = NULL, 
                          lambda1 = 1, lambda2 = 0, n_iter = 10, epsilon = 10**-6, 
                          output = "vec", epsilon0 = 10**-6, plot = FALSE, 
                          max_est = 100, min_est = -100, 
                          max_logit_pi = 100, min_logit_pi = -100) {
  
  if(is.null(g_list)) {
    g_list <- list()
    gm_list <- list()
    m_list <- list()
    for(iii in 1:length(unique(g_ind))) {
      g_list[[iii]] <- which(g_ind == iii)
      gm_list[[iii]] <- list()
      if(!is.null(m_ind)) {
        for(jjj in 1:(length(unique(m_ind)) - 1)) {
          if(iii == 1) 
            m_list[[jjj]] <- which(m_ind == jjj)
          gm_list[[iii]][[jjj]] <- which(g_ind == iii & m_ind == jjj)
        }
      }
    }
  }
  
  n_ <- length(y)
  n_group <- length(g_list)
  
  pi_est <- max(min(mean(y != 0), 1 - 10**-6), 10**-6)
  logit_pi_est <- log(pi_est/(1 - pi_est))
  m_est <- sapply(g_list, function(x, y) mean(log(y[x][y[x] > 0])), y = y)
  m_est[is.na(m_est)] <- 0
  if(!is.null(m_ind)){
    e_est <- sapply(m_list, function(x, y) mean(log(y[x][y[x] > 0])), y = y) - mean(m_est)
    e_est[is.na(e_est)] <- 0
  }
  z <- (y != 0)*1
  
  if(!is.null(m_ind)) {
    n_par <- 1 + length(m_est) + length(e_est)
    est_ <- c(logit_pi_est, m_est, e_est)
    e_est <- c(e_est, 0)
    lambda_est <- exp(m_est[g_ind] + e_est[m_ind])
  } else {
    n_par <- 1 + length(m_est)
    est_ <- c(logit_pi_est, m_est)
    lambda_est <- exp(m_est)[g_ind]
  }
  
  first_deriv <- numeric(n_par)
  second_deriv <- matrix(0, n_par, n_par)
  log_lik <- loss <- numeric(n_iter + 1)
  log_lik[1] <- -log(1 + exp(logit_pi_est))*n_ + sum((1 - z)*log(1 + exp(logit_pi_est - lambda_est)) + 
                                                       z*(logit_pi_est + y*log(lambda_est) - lambda_est))
  loss[1] <- -log_lik[1] + 
    lambda1*sum(abs(m_est %*% t(rep(1, n_group)) - rep(1, n_group) %*% t(m_est)))/2
  if(!is.null(m_ind) && lambda2 > 0)
    loss[1] <- loss[1] + lambda2*sum(abs(e_est))
  for(i in 1:n_iter) {
    
    est_0 <- est_
    exp_t_neg_l <- exp(logit_pi_est - lambda_est)
    expit_exp_t_neg_l <- 1/(1 + exp_t_neg_l)
    expit_exp_t_neg_l <- 1 - expit_exp_t_neg_l
    first_deriv[1] <- -n_*pi_est + sum(z + (1 - z)*expit_exp_t_neg_l)
    first_deriv_l_t_l <- (y - lambda_est)*z - lambda_est*(1 - z)*expit_exp_t_neg_l
    if(!is.null(m_ind)){
      first_deriv[-1] <- sapply(append(g_list, m_list), function(x) sum(first_deriv_l_t_l[x]))
    } else {
      first_deriv[-1] <- sapply(g_list, function(x) sum(first_deriv_l_t_l[x]))
    }
    
    second_deriv[1, 1] <- -pi_est*(1 - pi_est)*n_ + sum((1 - z)*expit_exp_t_neg_l*(1 - expit_exp_t_neg_l))
    second_deriv_pi_me <- -(1 - z)*expit_exp_t_neg_l*(1 - expit_exp_t_neg_l)*lambda_est
    second_deriv_me_me <- -(z + (1 - z)*expit_exp_t_neg_l*(1 - lambda_est*(1 - expit_exp_t_neg_l)))*lambda_est
    if(!is.null(m_ind)) {
      second_deriv[1, -1] <- second_deriv[-1, 1] <- 
        sapply(append(g_list, m_list), function(x) sum(second_deriv_pi_me[x]))
      diag(second_deriv)[-1] <- sapply(append(g_list, m_list), function(x) sum(second_deriv_me_me[x]))
      second_deriv[(1+length(m_est)+1):n_par, 1+1:length(m_est)] <- 
        sapply(gm_list, function(x) {sapply(x, function(z) sum(second_deriv_me_me[z]))})
      second_deriv[1+1:length(m_est), (1+length(m_est)+1):n_par] <- 
        t(second_deriv[(1+length(m_est)+1):n_par, 1+1:length(m_est)])
    } else {
      second_deriv[1, -1] <- second_deriv[-1, 1] <- 
        sapply(g_list, function(x) sum(second_deriv_pi_me[x]))
      diag(second_deriv)[-1] <- sapply(g_list, function(x) sum(second_deriv_me_me[x]))
    }
    
    if(lambda1 > 0) {
      b_mat <- abs(m_est %*% t(rep(1, n_group)) - rep(1, n_group) %*% t(m_est)) + epsilon0
      delta_mat <- lambda1/n_group/(n_group-1)*(diag(rowSums(1/b_mat)) - 1/b_mat)
      first_deriv[1:length(m_est)+1] <- first_deriv[1:length(m_est)+1] - as.numeric(delta_mat %*% m_est)
      second_deriv[1:length(m_est)+1, 1:length(m_est)+1] <- 
        second_deriv[1:length(m_est)+1, 1:length(m_est)+1] - delta_mat
    }
    
    if(lambda2 > 0 && (!is.null(m_ind))){
      c_mat <- abs(e_est[-length(e_est)]) + epsilon0
      first_deriv[(1+length(m_est)+1):n_par] <- first_deriv[(1+length(m_est)+1):n_par] - 1/c_mat*e_est[-length(e_est)]*lambda2
      diag(second_deriv[(1+length(m_est)+1):n_par, (1+length(m_est)+1):n_par]) <- 
        diag(second_deriv[(1+length(m_est)+1):n_par, (1+length(m_est)+1):n_par]) - 1/c_mat*lambda2
    }
    
    second_deriv_inv <- try(solve(second_deriv), silent = TRUE)
    lll <- 0
    while(class(second_deriv_inv)[1] == "try-error" && lll < 20) {
      second_deriv_inv <- try(solve(second_deriv - epsilon0*diag(n_par)*10**lll))
      lll <- lll + 1
    }
    if(class(second_deriv_inv)[1] == "try-error")
      break()
    est_ <- est_ - as.numeric(second_deriv_inv %*% first_deriv)
    est_[1] <- min(max(est_[1], min_logit_pi), max_logit_pi)
    est_[-1] <- sapply(est_[-1], function(x) min(max(x, min_est), max_est))
    logit_pi_est <- est_[1]
    pi_est <- 1/(1 + exp(-logit_pi_est))
    m_est <- est_[1:length(m_est)+1]
    if(!is.null(m_ind)) {
      e_est <- c(est_[(1+length(m_est)+1):n_par], 0)
      lambda_est <- exp(m_est[g_ind] + e_est[m_ind])
    } else {
      lambda_est <- exp(m_est)[g_ind]
    }
    log_lik[i + 1] <- -log(1 + exp(logit_pi_est))*n_ + sum((1 - z)*log(1 + exp(logit_pi_est - lambda_est)) + 
                                                             z*(logit_pi_est + y*log(lambda_est) - lambda_est))
    loss[i + 1] <- -log_lik[i + 1] +  
      lambda1*sum(abs(m_est %*% t(rep(1, n_group)) - rep(1, n_group) %*% t(m_est)))/2
    if(lambda2 > 0 && !is.null(m_ind)) 
      loss[i + 1] <- loss[i + 1] + lambda2*sum(abs(e_est))
    
    if(is.nan(loss[i + 1]) || sum(is.nan(est_)) > 0)
      break()
    if((loss[i] - loss[i + 1] < epsilon) || (sum(abs(est_ - est_0)) < epsilon)) {
      if(loss[i] < loss[i + 1])
        est_ <- est_0
      break()
    }
  }
  
  if(plot) {
    plot(log_lik[1:(i+1)])
    plot(loss[1:(i+1)])
  }
  
  if(output == "vec")
    return(c(1/(1+exp(-est_[1])), est_[-1]))
  return(list(pi = pi_est, m = m_est, e = ifelse(is.null(m_ind), NA, e_est), 
              loss = loss[1:(i+1)], log_lik = log_lik[1:(i+1)]))
}

###########################################################
### RZiMM-Count Model

#' A Regularized Zero-inflated Mixture Model for scRNA Data
#'   - for count data
#'   - Poisson or Negative Binomial
#'
#' @param x A n by p matrix, n cells x p genes. 
#' @param mix_ind Sample index, for adjusting within-sample correlation; if set to be NULL,
#'   the within-sample correlation will not be adjusted. 
#' @param n_group Number of groups for clustering. 
#' @param dist "zinb" for zero-inflated Negative Binomial or "zip" for zero-inflated Poisson. 
#' @param g_ind_init Initial cluster index. 
#' @param n_iter Maximum number of total iteration. 
#' @param n_iter_m Number of iteration for the majorization step. 
#' @param lambda1 Weight of regularization for group effects. 
#' @param lambda2 Weight of regularization for sample effects. 
#' @param epsilon0 A small number to addess matrix singularity. 
#' @param epsilon Iteration stops when the change from last iteration is smaller than this value.
#' @param epsilon_p A small number for calculating the number of parameters. 
#' @param print_e Whether to print the changes compared to the last iteration. 
#' @param ... Other parameters passed to \code{zinb_regulated()} or \code{zip_regulated()}. 
#'
#' @return Returns an \code{RZiMM-class} object.
#'
#' @importFrom methods new
#' @importFrom stats rnorm
#' @importFrom stats na.omit
#' @importFrom magic adiag
#' 
#' @seealso
#' \code{\link{RZiMM-class}}\cr
#'
#' @export
RZiMM_Count <- function(x, mix_ind = NULL, n_group = 4, 
                        dist = c("zinb", "zip")[1], 
                        g_ind_init = NULL, 
                        n_iter = 100, 
                        n_iter_m = 10, 
                        lambda1 = 100, 
                        lambda2 = 0, 
                        epsilon0 = 10**-6, 
                        epsilon = 10**-6, print_e = FALSE, 
                        epsilon_p = 10**-3, ...) {
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(is.null(g_ind_init)) 
    g_ind_init <- sample(n_group, n, replace = TRUE)
  
  g_est <- g_ind_init
  log_l <- rep(NA, n_iter)
  err_l <- rep(NA, n_iter)
  
  if(!is.null(mix_ind)) {
    n_mix <- length(unique(mix_ind))
    if(n_mix == 1)
      mix_ind = NULL
  } else {
    n_mix <- 1
  }
  
  if(dist == "zinb"){
    cov_est <- matrix(0, 1 + n_group + n_mix, p)
    cov_est[2, ] <- 1
  } else {
    cov_est <- matrix(0, n_group + n_mix, p)
  }
  
  z <- (x != 0)*1
  
  for(l in 1:n_iter) {
    
    #### From last iteration
    g_est0 <- g_est
    cov_est0 <- cov_est
    
    #### Update m_est2
    g_list <- list()
    gm_list <- list()
    m_list <- list()
    for(iii in 1:length(unique(g_est))) {
      g_list[[iii]] <- which(g_est == iii)
      gm_list[[iii]] <- list()
      if(!is.null(mix_ind)) {
        for(jjj in 1:(n_mix - 1)) {
          if(iii == 1) 
            m_list[[jjj]] <- which(mix_ind == jjj)
          gm_list[[iii]][[jjj]] <- which(g_est == iii & mix_ind == jjj)
        }
      }
    }
    
    if(dist == "zinb") {
      cov_est <- apply(x, 2, zinb_regulated, g_ind = g_est, m_ind = mix_ind, 
                       g_list = g_list, m_list = m_list, gm_list = gm_list, 
                       lambda1 = lambda1, lambda2 = lambda2, n_iter = n_iter_m, 
                       epsilon = epsilon, epsilon0 = epsilon0, plot = FALSE, ...)
    } else {
      cov_est <- apply(x, 2, zip_regulated, g_ind = g_est, m_ind = mix_ind, 
                       g_list = g_list, m_list = m_list, gm_list = gm_list, 
                       lambda1 = lambda1, lambda2 = lambda2, n_iter = n_iter_m, 
                       epsilon = epsilon, epsilon0 = epsilon0, plot = FALSE, ...)
    }
    
    #### Update g0
    if(dist == "zinb") {
      
      pi_est <- cov_est[1, ]
      alpha_est <- cov_est[2, ]
      log_alpha_est <- log(alpha_est)
      m_est <- cov_est[2 + 1:n_group, ]
      if(!is.null(mix_ind)) {
        
        e_est <- rbind(cov_est[1 + n_group + 2:n_mix, ], rep(0, p))
        p_lik <- apply(cbind(mix_ind, x, z), 1, function(x) {
          y <- x[1 + 1:p]
          z <- x[1 + p + 1:p]
          apply(exp(m_est + rep(1, n_group) %*% t(e_est[x[1], ])), 1, function(mu_) {
            one_plus_alpha_l <- 1 + alpha_est * mu_
            log_one_plus_alpha_l <- log(one_plus_alpha_l)
            aaa <- log_one_plus_alpha_l/alpha_est
            sum((1 - z)*log(1 - pi_est*(1 - exp(-aaa))) + 
                  z*(log(pi_est) + lgamma(y + 1/alpha_est) - lgamma(1/alpha_est) - aaa + 
                       y*(log_alpha_est + log(mu_) - log_one_plus_alpha_l)))
          })
        })
      } else {
        
        p_lik <- apply(cbind(x, z), 1, function(x) {
          y <- x[1:p]
          z <- x[p + 1:p]
          apply(exp(m_est), 1, function(mu_) {
            one_plus_alpha_l <- 1 + alpha_est * mu_
            log_one_plus_alpha_l <- log(one_plus_alpha_l)
            aaa <- log_one_plus_alpha_l/alpha_est
            sum((1 - z)*log(1 - pi_est*(1 - exp(-aaa))) + 
                  z*(log(pi_est) + lgamma(y + 1/alpha_est) - lgamma(1/alpha_est) - aaa + 
                       y*(log_alpha_est + log(mu_) - log_one_plus_alpha_l)))
          })
        })
      }
    } else {
      
      pi_est <- cov_est[1, ]
      m_est <- cov_est[1 + 1:n_group, ]
      if(!is.null(mix_ind)) {
        
        e_est <- rbind(cov_est[n_group + 2:n_mix, ], rep(0, p))
        p_lik <- apply(cbind(mix_ind, x, z), 1, function(x) {
          y <- x[1 + 1:p]
          z <- x[1 + p + 1:p]
          apply(exp(m_est + rep(1, n_group) %*% t(e_est[x[1], ])), 1, function(mu_) {
            sum((1 - z)*log(1 - pi_est*(1 - exp(-mu_))) + z*(log(pi_est) + y*log(mu_) - mu_))
          })
        })
      } else {
        
        p_lik <- apply(cbind(x, z), 1, function(x) {
          y <- x[1:p]
          z <- x[p + 1:p]
          apply(exp(m_est), 1, function(mu_) {
            sum((1 - z)*log(1 - pi_est*(1 - exp(-mu_))) + z*(log(pi_est) + y*log(mu_) - mu_))
          })
        })
      }
    }
    
    g_est <- apply(p_lik, 2, which.max)
    missing_group <- which(sapply(1:n_group, function(x) sum(x == g_est) == 0))
    if(length(missing_group) > 0) {
      g_est_x <- sample(n_group, n, replace = TRUE)
      g_est[g_est_x %in% missing_group] <- g_est_x[g_est_x %in% missing_group]
    }
    
    #### Calculate partial log-likelihood
    if(!is.null(mix_ind)) {
      mu_ <- exp(m_est[g_est, ] + e_est[mix_ind, ])
    } else {
      mu_ <- exp(m_est[g_est, ])
    }
    pi_est_mat <- rep(1, n) %*% t(pi_est)
    if(dist == "zinb") {
      
      alpha_est_mat <- rep(1, n) %*% t(alpha_est)
      one_plus_alpha_l <- 1 + alpha_est_mat * mu_
      log_one_plus_alpha_l <- log(one_plus_alpha_l)
      aaa <- log_one_plus_alpha_l/alpha_est_mat
      log_l[l] <- sum((1 - z)*log(1 - pi_est_mat*(1 - exp(-aaa))) +
                        z*(log(pi_est_mat) + lgamma(x + 1/alpha_est_mat) - lgamma(1/alpha_est_mat) - 
                             aaa + x*(log(alpha_est_mat) + log(mu_) - log_one_plus_alpha_l)))
    } else {
      
      log_l[l] <- sum((1 - z)*log(1 - pi_est_mat*(1 - exp(-mu_))) + 
                        z*(log(pi_est_mat) + x*log(mu_) - mu_))
    }
    
    if(dist == "zinb") {
      err_l[l] <- (sum(abs(cov_est0[-2, ] - cov_est[-2, ])) + 
                     sum(abs(log(cov_est0[2, ]) - log(cov_est[2, ]))))/p
    } else {
      err_l[l] <- (sum(abs(cov_est0 - cov_est)))/p
    }
    
    if(print_e) print(c(l, err_l[l]))
    if(err_l[l] < epsilon)
      break
  }
  
  ind_j <- list()
  for(j in 1:p)
    ind_j[[j]] <- which(x[, j] != 0)
  
  ppx <- 0
  for(j in 1:p)
    ppx <- ppx + (n_clust(m_est[1:n_group, j], epsilon_p) - 1 + n_mix - 1)*log(length(ind_j[[j]]))
  ppx <- ppx - 2*log_l[l]
  
  ppy <- 0
  for(j in 1:p)
    ppy <- ppy + (n_clust(m_est[1:n_group, j], epsilon_p) - 1 + n_mix - 1)*log(length(ind_j[[j]]))*log(log(p))
  ppy <- ppy - 2*log_l[l]
  
  importance_ <- apply(m_est[1:n_group, ], 2, function(x) {
    sum(abs(rep(1, length(x)) %*% t(x) - x %*% t(rep(1, length(x))))/2)
  })
  
  return(methods::new("RZiMM", cluster = g_est,
                      importance = importance_/n_group/(n_group - 1), 
                      param = list(pi = pi_est, m = m_est, 
                                   e = (if(n_mix == 1) {NA} else {e_est}), 
                                   alpha = (if(dist == "zinb") {alpha_est} else {NA})), 
                      info = list(bic = list(bic = ppx, mbic = ppy), 
                                  log_lik = stats::na.omit(log_l), error_traj = stats::na.omit(err_l), 
                                  lambda1 = lambda1, lambda2 = lambda2, n_group = n_group, 
                                  model.type = ifelse(dist == "zinb", "RZiMM-NB", "RZiMM-Pois"))))
}
