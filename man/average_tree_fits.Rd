\name{average_tree_fits}
\alias{average_tree_fits}

\title{
Summarize posterior distribution of the average fit of tree ensemble
}
\description{
Computes posterior summaries of average treatment effects based on the fitted BCF model. 
}
\usage{
average_tree_fits(fit,
                  type = c("mu","tau"),
                  X_cont = matrix(0, nrow = 1, ncol = 1),
                  X_cat = matrix(0, nrow = 1, ncol = 1),
                  quantiles = c(0.05, 0.95),
                  weights = NULL,
                  verbose = TRUE, 
                  print_every = floor(max(c(nrow(X_cont), nrow(X_cat)))/10))
}
\arguments{
  \item{fit}{Object returned by \code{\link{flexBCF}}.}
  \item{type}{Character which is equal to either "mu" or "tau". Determines which tree ensemble is used to make predictions before averaging.}
   \item{X_cont}{Matrix of continuous covariates for either the prognostic function mu (\code{type == "mu"}) or treatment effect function tau(\code{type == "tau"}) for the individuals in the group over which you wish to average. Note, predictors must be re-scaled to lie in the interval [-1,1]. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no continuous covariates}
     \item{X_cat}{Integer matrix of categorical covariates for either the prognostic function mu if (\code{type == "mu"}) or treatment effect function tau(\code{type == "tau"}) for the individuals in the group over which you wish to average. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no categorical covariates}
  \item{quantiles}{Vector of quantiles which you wish to compute. Default is \code{c(0.05, 0.95)} so that the function will return 90\% posterior credible intervals for the average effect.}
  \item{weights}{Vector of non-negative weights for each individual. Default is a \code{NULL}, which internally gets converted to a vector of all 1's.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor(max(c(nrow(X_cont), nrow(X_cat)))/10)} so that only 10 messages are printed.}
 
}
\details{
Returns posterior mean and quantiles for the conditional average treatment or prognostic effects averaged over the covariate values contained in \code{X_cont} and \code{X_cat}.
Use the argument \code{weights} to pass observation weights, which may be useful to extrapolate findings to a different population than the one observed at training time.
}
\value{
A list containing
\item{mean}{Posterior mean estimate of the desired average treatment (or prognostic) effect}
\item{quantiles}{Posterior quantiles of the desired average treatment (or prognostic) effect.}
}

}