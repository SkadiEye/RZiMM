#' A simulated dataset (normalized)
#'
#' A list containing:
#'   a data matrix: 1000 cells x 1000 genes; 
#'   a vector: cluster index;
#'   a vector: sample index
#'   a vector: indicator for important genes
#'
#' @format A list
#' \describe{
#'   \item{x}{a data matrix: 1000 cells x 1000 genes, normalized, many zeros; 
#'     The first 100 genes are DE genes for the cluster}
#'   \item{cluster}{cluster index}
#'   \item{sample_index}{sample index (batch effect)}
#'   \item{imp_gene_list}{indicator for important genes}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name sim_x 
#' @usage data(sim_x)
#' @format A list of data
NULL

#' A simulated dataset (count)
#'
#' A list containing:
#'   a data matrix: 1000 cells x 1000 genes; 
#'   a vector: cluster index;
#'   a vector: sample index
#'   a vector: indicator for important genes
#'
#' @format A list
#' \describe{
#'   \item{x}{a data matrix: 1000 cells x 1000 genes, normalized, many zeros; 
#'     The first 100 genes are DE genes for the cluster}
#'   \item{cluster}{cluster index}
#'   \item{sample_index}{sample index (batch effect)}
#'   \item{imp_gene_list}{indicator for important genes}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name sim_y 
#' @usage data(sim_y)
#' @format A list of data
NULL
