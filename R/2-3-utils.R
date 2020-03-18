optim_ind <- function(ind_pred, ind_true) {
  
  perm <- gtools::permutations(n = length(unique(ind_true)), 
                               r = length(unique(ind_true)), 
                               v = sort(unique(ind_true)))
  optim_acc <- 0
  out_p <- 0
  for(i in 1:dim(perm)[1]) {
    
    acc_i <- mean(perm[i, ind_pred] == ind_true)
    if(acc_i > optim_acc) {
      
      optim_acc <- acc_i
      out_p <- perm[i, ]
    }
  }
  return(list(perm = out_p, acc = optim_acc))
}

step2_sig <- function(x, ind_pred) {
  
  lm_pval <- numeric(dim(x)[2])
  for(i in 1:dim(x)[2]) {
    
    xi_ <- x[, i]
    if(length(unique(ind_pred[xi_ != 0])) > 1) {
      lm_mod <- lm(xi_ ~ factor(ind_pred))
      lm_mod0 <- lm(xi_ ~ 1)
      lm_pval[i] <- anova(lm_mod, lm_mod0)$`Pr(>F)`[2]
    } else {
      lm_pval[i] <- 1
    }
  }
  
  return(lm_pval)
}

perf_stat <- function(ind_est, imp_score, cut_off, true_ind, true_imp) {
  
  acc <- optim_ind(ind_est, group_ind)$acc
  pred_ <- ROCR::prediction(imp_score, true_imp)
  auc_ <- ROCR::performance(pred_, "auc")@y.values[[1]]
  acc_ <- ROCR::performance(pred_, "acc")
  spec_ <- ROCR::performance(pred_, "spec")
  sens_ <- ROCR::performance(pred_, "sens")
  acc_ <- acc_@y.values[[1]][which.min(abs(acc_@x.values[[1]] - cut_off))]
  spec_ <- spec_@y.values[[1]][which.min(abs(spec_@x.values[[1]] - cut_off))]
  sens_ <- sens_@y.values[[1]][which.min(abs(sens_@x.values[[1]] - cut_off))]
  
  return(data.frame(cls_acc = acc, auc = auc_, acc = acc_, 
                    spec = spec_, sens = sens_))
}
