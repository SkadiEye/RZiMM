n_clust <- function(v, cutoff) {
  
  v <- sort(v)
  x <- length(v)
  pre_v <- v[1]
  for(vi in v[-1]) {
    if(abs(vi - pre_v) < cutoff)
      x <- x-1
    pre_v <- vi
  }
  return(x)
}

RZiMM <- function(x, mix_ind = NULL, n_group = 4, 
                  g_ind_init = NULL, 
                  # prob_g = TRUE, 
                  # scale_x = TRUE, 
                  n_iter = 100, 
                  n_iter_m = 5, 
                  lambda1 = 100, 
                  lambda2 = 0, 
                  epsilon_l1 = 10**-6, 
                  epsilon_l2 = 10**-6, 
                  epsilon = 10**-6, 
                  epsilon_p = 10**-3, 
                  print_e = FALSE) {
    
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(is.null(mix_ind))
    mix_ind <- rep(1, n)
  n_mix <- length(unique(mix_ind))
  
  if(is.null(g_ind_init)) {
    g_ind_init <- matrix(rnorm(n_group*n), n_group, n)
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
        if(class(try(solve(MSM), silent = TRUE)) == "try-error")
          MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
        m_est[, j] <- solve(MSM) %*% 
          (d_mat[, ind_j[[j]]] %*% (x[ind_j[[j]], j]))
      }
      else {
        
        for(mmm in 1:n_iter_m) {
          bj <- abs(m_est[1:n_group, j] %*% t(rep(1, n_group)) - 
                      rep(1, n_group) %*% t(m_est[1:n_group, j])) + epsilon_l1
          if(n_mix > 1) {
            
            cj <- abs(m_est[n_group + 1:(n_mix-1), j]) + epsilon_l2
            MSM <- d_mat[, ind_j[[j]]] %*% t(d_mat[, ind_j[[j]]]) + 
              magic::adiag(lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj), lambda2 * diag(1/cj))
            if(class(try(solve(MSM), silent = TRUE)) == "try-error")
              MSM <- MSM + epsilon_l1 * diag(dim(MSM)[1])
            m_est[, j] <- solve(MSM) %*% (d_mat[, ind_j[[j]]] %*% (x[ind_j[[j]], j]))
          } else {
            
            MSM <- d_mat[, ind_j[[j]]] %*% t(d_mat[, ind_j[[j]]]) + 
              lambda1/n_group/(n_group-1) * (diag(rowSums(1/bj)) - 1/bj)
            if(class(try(solve(MSM), silent = TRUE)) == "try-error")
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
    # err_l[l] <- mean((m_est - m_est0)**2)
    
    if(print_e) print(c(l, err_l[l]))
    if(err_l[l] < epsilon)
      break
  }
  
  pp <- 0
  for(j in 1:p)
    pp <- pp + n_clust(m_est[1:n_group, j], epsilon_p) - 1
  
  bic_v00 <- -2*log_l[l] + (pp+p*(n_mix-1))*log(n) + 2*0*lchoose(p*(n_group - 1), pp)
  bic_v05 <- -2*log_l[l] + (pp+p*(n_mix-1))*log(n) + 2*0.5*lchoose(p*(n_group - 1), pp)
  bic_v10 <- -2*log_l[l] + (pp+p*(n_mix-1))*log(n) + 2*1*lchoose(p*(n_group - 1), pp)
  bic_v_chen <- -2*log_l[l] + (pp+p*(n_mix-1))*log(n)*log(log(p))
  
  importance_ <- apply(m_est[1:n_group, ], 2, function(x) {
    sum(abs(rep(1, length(x)) %*% t(x) - x %*% t(rep(1, length(x))))/2)
  })
  
  return(new("ZRiMM", cluster = apply(g_est, 2, which.max),
             importance = importance_/n_group/n_group0, 
             param = list(pi = pi_est, m = m_est, sigma_sq = sigma_sq_est), 
             info = list(bic = list(bic = bic_v00, ebic05 = bic_v05, ebic1 = bic_v10, bic_wang = bic_v_chen), 
                         log_lik = na.omit(log_l), error_traj = na.omit(err_l), 
                         lambda1 = lambda1, lambda2 = lambda2, n_group = n_group, 
                         model.type = ifelse(length(unique(mix_ind)) == 1, "ZRiMM-scRNA", "ZRiMM-Naive"))))
}
