###########################################################
### Generic Functions for an RZiMM object

#' @name rzi
#' @rdname rzi
#'
#' @title RZiMM Functions
#'
#' @description A collection of functions for an RZiMM object
#'
#' @param object An \code{RZiMM} object.
#'
#' @seealso
#' \code{\link{RZiMM-class}}\cr
#' \code{\link{RZiMM}}\cr
NULL

#' @rdname rzi
#' @export
setGeneric("rzimm_importance", function(object) {
  standardGeneric("rzimm_importance")
})

#' @rdname rzi
#' @export
setMethod("rzimm_importance", "RZiMM", function(object) {
  return(object@importance)
})

#' @rdname rzi
#' @export
setGeneric("rzimm_cluster", function(object) {
  standardGeneric("rzimm_cluster")
})

#' @rdname rzi
#' @export
setMethod("rzimm_cluster", "RZiMM", function(object) {
  return(object@cluster)
})

#' @rdname rzi
#' @export
setGeneric("rzimm_info", function(object) {
  standardGeneric("rzimm_info")
})

#' @rdname rzi
#' @export
setMethod("rzimm_info", "RZiMM", function(object) {
  return(object@info)
})

#' @rdname rzi
#' @export
setGeneric("rzimm_param", function(object) {
  standardGeneric("rzimm_param")
})

#' @rdname rzi
#' @export
setMethod("rzimm_param", "RZiMM", function(object) {
  return(object@param)
})

#' @rdname rzi
#' @export
setGeneric("rzimm_bic", function(object) {
  standardGeneric("rzimm_bic")
})

#' @rdname rzi
#' @export
setMethod("rzimm_bic", "RZiMM", function(object) {
  return(object@info$bic)
})

#' @rdname rzi
#' @export
setGeneric("rzimm_nonzero_prob", function(object) {
  standardGeneric("rzimm_nonzero_prob")
})

#' @rdname rzi
#' @export
setMethod("rzimm_nonzero_prob", "RZiMM", function(object) {
  return(object@param$pi)
})

#' @rdname rzi
#' @export
setGeneric("rzimm_m", function(object) {
  standardGeneric("rzimm_m")
})

#' @rdname rzi
#' @export
setMethod("rzimm_m", "RZiMM", function(object) {
  return(object@param$m[1:object@info$n_group, ])
})



