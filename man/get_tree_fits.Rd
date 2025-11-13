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
              verbose = TRUE, 
              print_every = NULL)
}
\arguments{
  \item{fit}{Object returned by \code{\link{flexBCF}}.}
  \item{type}{Character which is equal to either "mu" or "tau". Determines which tree ensemble is used to make predictions before averaging. For evaluations of the prognostic function (resp. conditional average treatment effect function), set \code{type == "mu"} (resp. \code{type == "tau"}).}
  \item{X_cont}{Matrix of continuous covariates for either the prognostic function \eqn{\mu}{mu} if (\code{type == "mu"}) or CATE function \eqn{\tau}{tau} (\code{type == "tau"}). Note, predictors must be re-scaled to lie in the interval [-1,1]. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no continuous covariates}
  \item{X_cat}{Integer matrix of categorical covariates for either the prognostic function \eqn{\mu}{mu} if (\code{type == "mu"}) or CATE function \eqn{\tau}{tau} (\code{type == "tau"}). Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no categorical covariates}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{A status message is printed every \code{print_every} MCMC iterations. Whens \code{print_every} is \code{NULL}, \code{get_tree_fits} sets the value internally to be 1/10 the number of MCMC iterations.}
}
\details{
To reduce its memory footprint, \code{\link{flexBCF}} only returns posterior samples of the trees used to approximate \eqn{\mu(X)}{mu(X)} and \eqn{\tau(X)}{tau(X)}.
It does \strong{not} return matrices containing posterior samples of these functions evaluated at each training observation.
The function \code{get_tree_fits()} takes an object returned by \code{\link{flexBCF}} and evaluates the trees at the covariate values supplied by \code{X_cont} and \code{X_cat}.
}
\value{
A matrix whose rows correspond to MCMC iterations and who columns correspond to observations (i.e., rows of \code{X_cont} or \code{X_cat}).
}
