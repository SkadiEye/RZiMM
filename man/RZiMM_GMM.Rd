% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-1-a-RZiMM-GMM.R
\name{RZiMM_GMM}
\alias{RZiMM_GMM}
\title{A Regularized Zero-inflated Mixture Model for scRNA Data
  - Gaussian Mixture Version}
\usage{
RZiMM_GMM(
  x,
  mix_ind = NULL,
  n_group = 4,
  g_ind_init = NULL,
  n_iter = 100,
  n_iter_m = 5,
  lambda1 = 100,
  lambda2 = 0,
  epsilon_l1 = 10^-6,
  epsilon_l2 = 10^-6,
  epsilon = 10^-6,
  epsilon_p = 10^-3,
  print_e = FALSE
)
}
\arguments{
\item{x}{A n by p matrix, n cells x p genes.}

\item{mix_ind}{Sample index, for adjusting within-sample correlation; if set to be NULL,
the within-sample correlation will not be adjusted.}

\item{n_group}{Number of groups for clustering.}

\item{g_ind_init}{Initial cluster index.}

\item{n_iter}{Maximum number of total iteration.}

\item{n_iter_m}{Number of iteration for the majorization step.}

\item{lambda1}{Weight of regularization for group effects.}

\item{lambda2}{Weight of regularization for sample effects.}

\item{epsilon_l1}{A small number to addess matrix singularity.}

\item{epsilon_l2}{A small number to addess matrix singularity.}

\item{epsilon}{Iteration stops when the change from last iteration is smaller than this value.}

\item{epsilon_p}{A small number for calculating the number of parameters.}

\item{print_e}{Whether to print the changes compared to the last iteration.}
}
\value{
Returns an \code{RZiMM-class} object.
}
\description{
A Regularized Zero-inflated Mixture Model for scRNA Data
  - Gaussian Mixture Version
}
\seealso{
\code{\link{RZiMM-class}}\cr
}
