###########################################################
### Define generic functions: show, print, '[', '$', plot

#################################
#### RZiMM class

#' @name show
#' @rdname show
#'
#' @title Method show for the package
#'
#' @description The method show for \code{RZiMM} object.
#'
#' @param object A \code{RZiMM} object.
#'
#' @seealso
#' \code{\link{RZiMM-class}}\cr
#'
NULL

#' @rdname show
#' @importFrom methods show
#'
#' @export
setMethod("show",
          "RZiMM",
          function(object) {
            
            cat("A", object@info$model.type, "Object \n")
            n <- length(object@cluster)
            p <- dim(object@importance)
            n_group <- object@info$n_group
            
            cat("Cluster: (first 6 cells) \n")
            print(object@cluster[1:min(6, n)])
            
            cat("Importance: (first 6 genes) \n")
            print(round(object@importance[1:min(6, p)], 3))
            
            cat("Group mean: (first 6 genes) \n")
            print(round(object@param$m[1:n_group, 1:min(6, p)], 3))
          })



