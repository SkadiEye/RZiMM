###########################################################
### RZiMM-PF Model

#' A Regularized Zero-inflated Mixture Model for scRNA Data
#'   - with PermFIT-type p-values
#'
#' @param x A n by p matrix, n cells x p genes. 
#' @param mix_ind Sample index, for adjusting within-sample correlation; if set to be NULL,
#'   the within-sample correlation will not be adjusted. 
#' @param n_group Number of groups for clustering. 
#' @param g_ind_init Initial cluster index. 
#' @param n_iter Maximum number of total iteration. 
#' @param n_iter_m Number of iteration for the majorization step. 
#' @param lambda1 Weight of regularization for group effects. 
#' @param lambda2 Weight of regularization for sample effects. 
#' @param epsilon_l1 A small number to addess matrix singularity. 
#' @param epsilon_l2 A small number to addess matrix singularity. 
#' @param epsilon Iteration stops when the change from last iteration is smaller than this value.
#' @param epsilon_p A small number for calculating the number of parameters. 
#' @param print_e Whether to print the changes compared to the last iteration. 
#' @param n_perm Number of permutations for PermFIT p-values
#' @param p_perm Proportion of samples for the training set for PermFIT p-values
#'
#' @return Returns an \code{RZiMM-class} object.
#'
#' @importFrom methods new
#' @importFrom stats rnorm
#' @importFrom stats na.omit
#' @importFrom magic adiag
#' @importFrom stats coefficients
#' @importFrom stats vcov
#' @importFrom stats pnorm
#' @importFrom stats var
#' 
#' @seealso
#' \code{\link{RZiMM-class}}\cr
#'
#' @export
RZiMM_PF <- function(x, mix_ind = NULL, n_group = 4, 
                     g_ind_init = NULL, 
                     n_iter = 100, 
                     n_iter_m = 5, 
                     lambda1 = 100, 
                     lambda2 = 0, 
                     epsilon_l1 = 10**-6, 
                     epsilon_l2 = 10**-6, 
                     epsilon = 10**-6, 
                     epsilon_p = 10**-3, 
                     print_e = FALSE, 
                     n_perm = 100, 
                     p_perm = 0.5) {
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(is.null(mix_ind))
    mix_ind <- rep(1, n)
  n_mix <- length(unique(mix_ind))
  
  if(is.null(g_ind_init)) {
    g_ind_init <- matrix(stats::rnorm(n_group*n), n_group, n)
    g_ind_init <- (g_ind_init == (rep(1, n_group) %*% t(apply(g_ind_init, 2, max))))*1
  }
  
  pi_est <- colMeans(x != 0)
  
  sigma_sq_est <- rep(1, p)
  m_est <- matrix(0, n_group + n_mix - 1, p)
  g_est <- g_ind_init
  log_l <- rep(NA, n_iter)
  err_l <- rep(NA, n_iter)
  ind_j <- list()
  for(j in 1:p)
    ind_j[[j]] <- which(x[, j] != 0)
  
  if(n_mix > 1) {
    
    c_mat <- matrix(0, n_mix-1, n)
    for(j in 1:(n_mix-1))
      c_mat[j, mix_ind == j] <- 1
  }
  
  if(n_mix > 1)
    d_mat <- rbind(g_est, c_mat)
  else 
    d_mat <- g_est
  
  for(l in 1:n_iter) {
    
    #### From last iteration
    g_est0 <- g_est
    m_est0 <- m_est
    
    #### Update m_est2
    for(j in 1:p) {
      
      if(lambda1 == 0 && lambda2 == 0) {
        
        MSM <- d_mat[, ind_j[[j]]] %*% t(d_mat[, ind_j[[j]]])
        if(class(try(solve(MSM), silent = TRUE))[1] == "try-error")
          MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
        m_est[, j] <- solve(MSM) %*% 
          (d_mat[, ind_j[[j]]] %*% (x[ind_j[[j]], j]))
      } else {
        
        for(mmm in 1:n_iter_m) {
          bj <- abs(m_est[1:n_group, j] %*% t(rep(1, n_group)) - 
                      rep(1, n_group) %*% t(m_est[1:n_group, j])) + epsilon_l1
          if(n_mix > 1) {
            
            cj <- abs(m_est[n_group + 1:(n_mix-1), j]) + epsilon_l2
            if(n_mix > 2) {
              MSM <- d_mat[, ind_j[[j]]] %*% t(d_mat[, ind_j[[j]]]) + 
                magic::adiag(lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj), lambda2 * diag(1/cj))
            } else {
              MSM <- d_mat[, ind_j[[j]]] %*% t(d_mat[, ind_j[[j]]]) + 
                magic::adiag(lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj), lambda2 * matrix(1/cj))
            }
            if(class(try(solve(MSM), silent = TRUE))[1] == "try-error")
              MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
            m_est[, j] <- solve(MSM) %*% (d_mat[, ind_j[[j]]] %*% (x[ind_j[[j]], j]))
          } else {
            
            MSM <- d_mat[, ind_j[[j]]] %*% t(d_mat[, ind_j[[j]]]) + 
              lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj)
            if(class(try(solve(MSM), silent = TRUE))[1] == "try-error")
              MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
            m_est[, j] <- solve(MSM) %*% (d_mat[, ind_j[[j]]] %*% (x[ind_j[[j]], j]))
          }
        }
      }
    }
    
    #### Update mu_mat, sigma_sq_est
    for(j in 1:p) 
      sigma_sq_est[j] <- mean(((x[ind_j[[j]], j] - 
                                  as.numeric(t(d_mat[, ind_j[[j]]]) %*% m_est[, j])))**2)
    
    #### Update g0
    for(i in 1:n) {
      mu_group <- m_est[1:n_group, ]
      if(mix_ind[i] < n_mix)
        err <- (rep(1, n_group) %*% t(x[i, ] - m_est[n_group + mix_ind[i], ]) - mu_group)**2 / 
          (rep(1, n_group) %*% t(sigma_sq_est)) / 2
      else 
        err <- (rep(1, n_group) %*% t(x[i, ]) - mu_group)**2 / 
          (rep(1, n_group) %*% t(sigma_sq_est)) / 2
      err <- rowSums(err[, x[i, ] != 0])
      err <- err - min(err)
      g_est[, i] <- 0 
      g_est[which.min(err), i] <- 1 
    }
    
    #### Calculate partial log-likelihood
    if(n_mix > 1)
      d_mat <- rbind(g_est, c_mat)
    else 
      d_mat <- g_est
    
    log_l[l] <- 0
    for(j in 1:p)
      log_l[l] <- log_l[l] - 0.5*log(sigma_sq_est[j])*length(ind_j[[j]]) - 
      sum((x[ind_j[[j]], j] - as.numeric(t(d_mat[, ind_j[[j]]]) %*% m_est[, j]))**2/2/sigma_sq_est[j])
    
    err_l[l] <- mean(apply(g_est, 2, which.max) != apply(g_est0, 2, which.max))
    
    if(print_e) print(c(l, err_l[l]))
    if(err_l[l] < epsilon)
      break
  }
  
  #### p-vals ####
  pval_mat <- matrix(NA, n_group, p)
  coef_mat <- matrix(NA, n_group, p)
  for(j in 1:p) {
    
    xj <- x[ind_j[[j]], j]
    mix_j <- mix_ind[ind_j[[j]]]
    for(k in 1:n_group) {
      
      ind_jk <- g_est[k, ind_j[[j]]]
      if(n_mix > 1) l_mod <- lm(xj ~ ind_jk + factor(mix_j))
      if(n_mix == 1) l_mod <- lm(xj ~ ind_jk)
      coef_mat[k, j] <- coefficients(l_mod)[2]/sqrt(vcov(l_mod)[2, 2])
      pval_mat[k, j] <- pnorm(abs(coefficients(l_mod)[2]), sd = sqrt(vcov(l_mod)[2, 2]), lower.tail = FALSE)*2
    }
  }
  
  #### Sampling 
  training_ind <- sample(n, p_perm*n)
  validate_ind <- (1:n)[-training_ind]
  m0_est <- matrix(0, n_group + n_mix - 1, p)
  mse_ <- list()
  
  #### Update m_est2
  for(j in 1:p) {
    
    ind_j_j <- intersect(training_ind, ind_j[[j]])
    ind_j_j_valid <- intersect(validate_ind, ind_j[[j]])
    if(lambda1 == 0 && lambda2 == 0) {
      
      MSM <- d_mat[, ind_j_j] %*% t(d_mat[, ind_j_j])
      if(class(try(solve(MSM), silent = TRUE))[1] == "try-error")
        MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
      m0_est[, j] <- solve(MSM) %*% 
        (d_mat[, ind_j_j] %*% (x[ind_j_j, j]))
    } else {
      
      for(mmm in 1:n_iter_m) {
        bj <- abs(m0_est[1:n_group, j] %*% t(rep(1, n_group)) - 
                    rep(1, n_group) %*% t(m0_est[1:n_group, j])) + epsilon_l1
        if(n_mix > 1) {
          
          cj <- abs(m0_est[n_group + 1:(n_mix-1), j]) + epsilon_l2
          if(n_mix > 2) {
            MSM <- d_mat[, ind_j_j] %*% t(d_mat[, ind_j_j]) + 
              magic::adiag(lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj), lambda2 * diag(1/cj))
          } else {
            MSM <- d_mat[, ind_j_j] %*% t(d_mat[, ind_j_j]) + 
              magic::adiag(lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj), lambda2 * matrix(1/cj))
          }
          if(class(try(solve(MSM), silent = TRUE))[1] == "try-error")
            MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
          m0_est[, j] <- solve(MSM) %*% (d_mat[, ind_j_j] %*% (x[ind_j_j, j]))
        } else {
          
          MSM <- d_mat[, ind_j_j] %*% t(d_mat[, ind_j_j]) + 
            lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj)
          if(class(try(solve(MSM), silent = TRUE))[1] == "try-error")
            MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
          m0_est[, j] <- solve(MSM) %*% (d_mat[, ind_j_j] %*% (x[ind_j_j, j]))
        }
      }
    }
    
    mse0_j <- (x[ind_j_j_valid, j] - as.numeric(t(d_mat[, ind_j_j_valid]) %*% m0_est[, j]))**2
    mse_j <- matrix(NA, length(ind_j_j_valid), n_perm)
    for(k in 1:n_perm) {
      
      shuffle_k <- sample(ind_j_j_valid)
      d_mat_k <- rbind(d_mat[1:n_group, shuffle_k], d_mat[-(1:n_group), ind_j_j_valid])
      mse_j[, k] <- (x[ind_j_j_valid, j] - as.numeric(t(d_mat_k) %*% m0_est[, j]))**2 - mse0_j
    }
    
    mse_[[j]] <- rowMeans(mse_j)
  }
  
  imp_ <- sapply(mse_, mean)
  imp_sd_ <- sapply(mse_, function(z) sqrt(var(z)/length(z)))
  imp_pval_ <- pnorm(imp_/imp_sd_, lower.tail = FALSE)
  
  overall_imp_0 <- sum(imp_)
  overall_imp_x <- sum(imp_/imp_sd_)
  
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
  
  return(methods::new("RZiMM", cluster = apply(g_est, 2, which.max),
                      importance = importance_/n_group/(n_group - 1), 
                      param = list(pi = pi_est, m = m_est, sigma_sq = sigma_sq_est), 
                      info = list(bic = list(bic = ppx, mbic = ppy), 
                                  log_lik = stats::na.omit(log_l), error_traj = stats::na.omit(err_l), 
                                  lambda1 = lambda1, lambda2 = lambda2, n_group = n_group, 
                                  model.type = "RZiMM-PermFIT", 
                                  importance = data.frame(importance = imp_, 
                                                          importance_sd = imp_sd_, 
                                                          importance_pval = imp_pval_), 
                                  imp_index = list(overall_imp_0 = overall_imp_0, overall_imp_x = overall_imp_x), 
                                  ind_imp = list(z.stat = coef_mat, pval = pval_mat))))
}

