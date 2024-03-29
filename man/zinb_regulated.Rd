% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2-1-c-RZiMM-Count.R
\name{zinb_regulated}
\alias{zinb_regulated}
\title{A Zero-Inflated Negative Binomial Model
  - for RZiMM-Count Model}
\usage{
zinb_regulated(
  y,
  g_ind,
  m_ind = NULL,
  g_list = NULL,
  m_list = NULL,
  gm_list = NULL,
  lambda1 = 1,
  lambda2 = 0,
  n_iter = 10,
  epsilon = 10^-6,
  output = "vec",
  epsilon0 = 10^-6,
  plot = FALSE,
  max_l_alpha = 10,
  min_l_alpha = -10,
  max_est = 100,
  min_est = -100,
  max_logit_pi = 100,
  min_logit_pi = -100
)
}
\arguments{
\item{y}{A vector of expression levels}

\item{g_ind}{A vector of cell-type clusters}

\item{m_ind}{Sample index, for adjusting within-sample correlation; if set to be NULL,
the within-sample correlation will not be adjusted.}

\item{g_list}{A list for cell-type cluster indices}

\item{m_list}{A list for batch indices}

\item{gm_list}{A two-way list for cell-type cluster and batch indices}

\item{lambda1}{Weight of regularization for group effects.}

\item{lambda2}{Weight of regularization for sample effects.}

\item{n_iter}{Maximum number of total iteration.}

\item{epsilon}{Iteration stops when the change from last iteration is smaller than this value.}

\item{output}{Whether export a vector or a list}

\item{epsilon0}{A small number for calculating the number of parameters.}

\item{plot}{Whether plotting the trajectories for the loss and log-likelihood}

\item{max_l_alpha}{Upper threshold for l}

\item{min_l_alpha}{Lower threshold for l}

\item{max_est}{Upper threshold for estimate}

\item{min_est}{Lower threshold for estimate}

\item{max_logit_pi}{Upper threshold for pi}

\item{min_logit_pi}{Lower threshold for pi}
}
\value{
Returns a vector of parameters.
}
\description{
A Zero-Inflated Negative Binomial Model
  - for RZiMM-Count Model
}
