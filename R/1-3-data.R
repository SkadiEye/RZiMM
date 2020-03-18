#' A simulated dataset 
#'
#' A list containing:
#'   a data matrix: 1000 cells x 1000 genes; 
#'   a vector: cluster index;
#'   a vector: sample index
#'   a vector: indicator for important genes
#'
#' @format A list
#' \describe{
#'   \item{x}{a data matrix: 1000 cells x 1000 genes, normalized, many zeros}
#'   \item{cluster}{cluster index}
#'   \item{sample_index}{sample index}
#'   \item{imp_gene_list}{indicator for important genes}
#' }
"sim_x"

