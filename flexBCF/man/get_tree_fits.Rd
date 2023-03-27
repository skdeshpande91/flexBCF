\name{get_tree_fits}
\alias{get_tree_fits}

\title{
Get fits of regression tree ensembles
}
\description{
Computes posterior samples of evaluations of regression tree ensembles
}
\usage{
get_tree_fits(fit,
              type = c("mu","tau"),
              X_cont = matrix(0, nrow = 1, ncol = 1),
              X_cat = matrix(0, nrow = 1, ncol = 1),
              verbose = TRUE, print_every = floor(max(c(nrow(X_cont), nrow(X_cat)))/10))
}
\arguments{
  \item{fit}{Object returned by \code{flexBCF}.}
  \item{type}{Character which is equal to either "mu" or "tau". Determines which tree ensemble is used to make predictions before averaging. For inference about the ATE, set \code{type == "tau"}.}
  \item{X_cont}{Matrix of continuous covariates for either the prognostic function mu if (\code{type == "mu"}) or treatment effect function tau(\code{type == "tau"}) for the individuals in the group over which you wish to average. Note, predictors must be re-scaled to lie in the interval [-1,1]. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no continuous covariates}
  \item{X_cat}{Integer matrix of categorical covariates for either the prognostic function mu if (\code{type == "mu"}) or treatment effect function tau(\code{type == "tau"}) for the individuals in the group over which you wish to average. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no categorical covariates}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor(max(c(nrow(X_cont), nrow(X_cat)))/10)} so that only 10 messages are printed.}
}
\details{
TO DO: more details about our implementation
}
\value{
A matrix whose rows correspond to MCMC iterations and who columns correspond to individuals
}
\examples{
\dontrun{
1+1
}
}