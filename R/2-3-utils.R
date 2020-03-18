###########################################################
### A function to match cluster numbers

#' Match the cluster number of two clustering results, return the 
#'  matching numbers and the matching accuracy. 
#'
#' @param ind_pred Clustering indices from one way.  
#' @param ind_true Clustering indices from another way.  
#'
#' @return Returns a list.
#' 
#' @importFrom gtools permutations
#'
#' @export
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

###########################################################
### Anova test in a two-step algorithm

#' Anova test followed by a clustering algorithm, e.g.
#'   K-Means, H-Clust, etc
#'
#' @param x Data matrix. 
#' @param ind_pred Clustering result. 
#'
#' @return Returns a vector of p-values. 
#'
#' @importFrom stats lm
#' @importFrom stats anova
#' 
#' @export
step2_sig <- function(x, ind_pred) {
  
  lm_pval <- numeric(dim(x)[2])
  for(i in 1:dim(x)[2]) {
    
    xi_ <- x[, i]
    if(length(unique(ind_pred[xi_ != 0])) > 1) {
      lm_mod <- stats::lm(xi_ ~ factor(ind_pred))
      lm_mod0 <- stats::lm(xi_ ~ 1)
      lm_pval[i] <- stats::anova(lm_mod, lm_mod0)$`Pr(>F)`[2]
    } else {
      lm_pval[i] <- 1
    }
  }
  
  return(lm_pval)
}

###########################################################
### Model Performance

#' Model performance in a simulated scenario
#'
#' @param ind_est Clustering result from an algorithm. 
#' @param imp_score Importance scores from an algorithm. 
#' @param cut_off Cutoff for importance score. 
#' @param true_ind True clustering indices. 
#' @param true_imp True importance (0 for unimportant; 1 for important). 
#'
#' @return Returns a list of performance scores. 
#' 
#' @importFrom ROCR prediction
#' @importFrom ROCR performance
#'
#' @export
perf_stat <- function(ind_est, imp_score, cut_off, true_ind, true_imp) {
  
  acc <- optim_ind(ind_est, true_ind)$acc
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

###########################################################
### Plot the Importance

#' Plot the importance from an RZiMM object
#'
#' @param rzimm_mod An RZiMM object. 
#' @param cutoff Cutoff to highlight important genes. 
#' @param top Alternatively, one can choose to highlight top x genes. 
#'
#' @return Returns a ggplot object.
#' 
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#'
#' @seealso
#' \code{\link{RZiMM-class}}\cr
#'
#' @export
importance_plot <- function(rzimm_mod, cutoff = NULL, top = NULL) {
  
  imp <- rzimm_importance(rzimm_mod)
  if(is.null(cutoff)) {
    if(is.null(top)) {
      color <- ifelse(rank(-imp) <= floor(length(imp)/10), "type1", "type2")
    } else {
      color <- ifelse(rank(-imp) <= top, "type1", "type2")
    }
  } else {
    color <- ifelse(imp >= cutoff, "type1", "type2")
  }
  imp_dat <- data.frame(ind = 1:length(imp), 
                        imp = imp, 
                        color = color)
  ggplot2::ggplot(imp_dat, ggplot2::aes(x = ind, y = imp)) + 
    ggplot2::geom_segment(ggplot2::aes(x = ind, xend = ind, y = 0, 
                                       yend = imp, color = color), size = 0.7, alpha = 1) + 
    ggplot2::scale_color_manual(values = c("#E65F00", "#56B4E9")) + 
    ggplot2::theme_light() + ggplot2::ylab("Importance") + ggplot2::xlab("Gene") + 
    ggplot2::theme(legend.position = "none", 
                   axis.ticks.x = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_blank())
}
