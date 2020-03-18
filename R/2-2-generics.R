setGeneric("rzimm_importance", function(object) {
  standardGeneric("rzimm_importance")
})

setMethod("rzimm_importance", "ZRiMM", function(object) {
  return(object@importance)
})

setGeneric("rzimm_cluster", function(object) {
  standardGeneric("rzimm_cluster")
})

setMethod("rzimm_cluster", "ZRiMM", function(object) {
  return(object@cluster)
})

setGeneric("rzimm_info", function(object) {
  standardGeneric("rzimm_info")
})

setMethod("rzimm_info", "ZRiMM", function(object) {
  return(object@info)
})

setGeneric("rzimm_param", function(object) {
  standardGeneric("rzimm_param")
})

setMethod("rzimm_param", "ZRiMM", function(object) {
  return(object@param)
})

setGeneric("rzimm_bic", function(object) {
  standardGeneric("rzimm_bic")
})

setMethod("rzimm_bic", "ZRiMM", function(object) {
  return(object@info$bic)
})

setGeneric("rzimm_nonzero_prob", function(object) {
  standardGeneric("rzimm_nonzero_prob")
})

setMethod("rzimm_nonzero_prob", "ZRiMM", function(object) {
  return(object@param$pi)
})

setGeneric("rzimm_m", function(object) {
  standardGeneric("rzimm_m")
})

setMethod("rzimm_m", "ZRiMM", function(object) {
  return(object@param$m[1:object@info$n_group, ])
})



