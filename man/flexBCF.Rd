\name{flexBCF}
\alias{flexBCF}

\title{
A faster and more flexible Bayesian Causal Forests
}
\description{
Implements a version of Hahn et al. (2020)'s Bayesian Causal Forest (BCF) model that is faster and handles categorical outcomes more flexibly than existing BCF implementations. 
This implementation was created for the 2022 American Causal Inference Conference's Data Challenge.
}
\usage{
flexBCF(Y_train,
        treated,
        X_cont_mu = matrix(0, nrow = 1, ncol = 1),
        X_cat_mu = matrix(0, nrow = 1, ncol = 1),
        X_cont_tau = matrix(0, nrow = 1, ncol = 1),
        X_cat_tau = matrix(0, nrow = 1, ncol = 1),
        unif_cuts_mu = rep(TRUE, times = ncol(X_cont_mu)),
        unif_cuts_tau = rep(TRUE, times = ncol(X_cont_tau)),
        cutpoints_list_mu = NULL,
        cutpoints_list_tau = NULL,
        cat_levels_list_mu = NULL,
        cat_levels_list_tau = NULL,
        sparse = TRUE, 
        M_mu = 50, M_tau = 50,
        nd = 1000, burn = 1000, thin = 1,
        verbose = TRUE, print_every = floor( (nd*thin + burn))/10)
}
\arguments{
  \item{Y_train}{Vector of observed outcomes}
  \item{treated}{Vector of treatment assignments. Assumes that all control subjects are listed first (treated ==0) and then all treated subjects are listed (treated == 1).}
  \item{X_cont_mu}{Matrix of continuous covariates for prognostic function mu. Note, predictors must be re-scaled to lie in the interval [-1,1]. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no continuous covariates for mu.}
  \item{X_cat_mu}{Integer matrix of categorical covariates for prognostic function mu. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no categorical covariates for mu.}
  \item{X_cont_tau}{Matrix of continuous covariates for treatment effect function mu. Note, predictors must be re-scaled to lie in the interval [-1,1]. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no continuous covariates for tau.}
  \item{X_cat_tau}{Integer matrix of categorical covariates for treatment effect function tau. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Assumes that subjects are ordered so that all control subjects are listed before all treated subjects. Default is a 1x1 matrix, which signals that there are no categorical covariates for tau.}
  \item{unif_cuts_mu}{Vector of logical values indicating whether cutpoints for each continuous covariate for mu should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list_mu}. Default is \code{TRUE} for each variable in \code{X_cont_mu}}
   \item{unif_cuts_tau}{Vector of logical values indicating whether cutpoints for each continuous covariate for tau should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list_tau}. Default is \code{TRUE} for each variable in \code{X_cont_mu}}
  \item{cutpoints_list_mu}{List of length \code{ncol(X_cont_mu)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{cutpoints_list_tau}{List of length \code{ncol(X_cont_tau)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{cat_levels_list_mu}{List of length \code{ncol(X_cat_train)} containing a vector of levels for each categorical predictor. If the j-th categorical predictor contains L levels, \code{cat_levels_list[[j]]} should be the vector \code{0:(L-1)}. Default is \code{NULL}, which corresponds to the case that no categorical predictors are available.}
  \item{cat_levels_list_tau}{List of length \code{ncol(X_cat_train)} containing a vector of levels for each categorical predictor. If the j-th categorical predictor contains L levels, \code{cat_levels_list[[j]]} should be the vector \code{0:(L-1)}. Default is \code{NULL}, which corresponds to the case that no categorical predictors are available.}
  \item{sparse}{Logical, indicating whether or not to perform variable selection based on a sparse Dirichlet prior rather than uniform prior; see Linero 2018. Default is \code{FALSE}.}
  \item{M_mu}{Number of trees in the ensemble for prognostic function mu. Default is 50.}
  \item{M_tau}{Number of trees in the ensemble for treatment effect function tau. Default is 50.}
  \item{nd}{Number of posterior draws to return. Default is 1000.}
  \item{burn}{Number of MCMC iterations to be treated as "warmup" or "burn-in". Default is 1000.}
  \item{thin}{Number of post-warmup MCMC iteration by which to thin. Default is 1.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.}
}
\details{
The Bayesian Causal Forest (BCF) model of Hahn et al. (2022) expresses the conditional expectation of a response $Y$ given covariates $X$ and binary treament $Z$ as
\eqn{E[Y | X, Z] = \mu(X) + \tau(X) * Z}{E[Y|X,Z] = mu(X) + tau(X) * Z} and approximates both \eqn{\mu(X)}{mu(X)} and \eqn{\tau(X)}{tau(X)} with ensembles of binary regression trees.
Under standard identifying assumptions, \eqn{\mu}{mu} is the prognostic function \eqn{\mu(X) = E[Y(0) \vert X =x]}{mu(X) = E[Y(0) | X = x]} and \eqn{\tau}{tau} is the conditional average treatment effect function \eqn{\tau(X) = E[Y(1) - Y(0) \vert X]}{tau(X) = E[Y(1) - Y(0)|X]}, where \eqn{Y(0)}{Y(0)} and \eqn{Y(1)}{Y(1)} are the potential outcomes under control and treatment.

Like the original BCF implementation, \code{flex_BCF()} simulates a Markov chain (written in C++) to compute approximate posterior draws of \eqn{\mu(X)}(mu(X)) and \eqn{\tau(X)}{tau(X)}.
By avoiding several redundant internal computations and memory allocations, \code{flex_BCF()} is much faster and scales more gracefully to large datasets than the original BCF implementation.
It is also based on a new prior for regression trees that permits more flexible \dQuote{borrowing of strength} across levels of categorical predictors.
Under the new prior conditional on splitting on a categorical predictor at a particular node in the tree, levels of the predictor are sent to the left and right child uniformly at random. 

To reduce its memory footprint, \code{flex_BCF()} only returns posterior samples of the trees used to approximate \eqn{\mu(X)}{mu(X)} and \eqn{\tau(X)}{tau(X)}.
It does \strong{not} return matrices containing posterior samples of these functions evaluated at each training observation.
Use the functions \code{\link{get_tree_fits}} and \code{\link{average_tree_fits}}, respectively, to obtain posterior samples of these function evaluations and to compute posterior summaries of weighted averages of these evaluations.
}

\value{
A list containing
\item{sigma}{Vector containing ALL samples of the residual standard deviation, including burn-in.}
\item{mu_trees}{A list (of length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees in the ensemble for the prognostic function mu. These strings are parsed by \code{get_tree_fits} and \code{average_tree_fits} to reconstruct the C++ representations of the sampled trees.}
\item{tau_trees}{A list (of length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees in the ensemble for the prognostic function mu. These strings are parsed by \code{get_tree_fits} and \code{average_tree_fits} to reconstruct the C++ representations of the sampled trees.}
\item{y_mean}{Mean of the observed responses. Used by \code{get_tree_fits} and \code{average_tree_fits} to transform predictions to original scale.}
\item{y_sd}{Standard deviation of the observed responses. Used by \code{get_tree_fits} and \code{average_tree_fits} to transform predictions to original scale.}
\item{varcounts_mu}{Matrix that counts the number of times a covariate was used in a decision rule in the ensemble used to approximate treatment effect function tau in each MCMC iteration. Rows correspond to MCMC iterations and columns corresponding to covariates}
\item{varcounts_tau}{Matrix that counts the number of times a covariate was used in a decision rule in the ensemble used to approximate treatment effect function tau in each MCMC iteration. Rows correspond to MCMC iterations and columns corresponding to covariates}
\item{cat_levels_list}{List of containing \code{cat_levels_list_mu} and \code{cat_levels_list_tau}. Used by \code{get_tree_fits} and \code{average_tree_fits}.}
}

\seealso{
\code{\link{get_tree_fits}} to get posterior samples of evaluations of \eqn{\mu}{mu} and \eqn{\tau}{tau} and \code{\link{average_tree_fits}} to compute posterior means and quantiles of weighted averages of such evaluations.
}

\examples{
\dontrun{
p_cont_mu <- 5
p_cat_mu <- 5
p_mu <- p_cont_mu + p_cat_mu

p_cont_tau <- 5
p_cat_tau <- 5
p_tau <- p_cont_tau + p_cat_tau

mu_func <- function(X_cont_mu, X_cat_mu){
  scaled_X_cont <- (X_cont_mu + 1)/2 # moves it to [0,1]
  tmp1 <- sin(pi * scaled_X_cont[,1] * scaled_X_cont[,2])
  tmp2 <- (scaled_X_cont[,3] - 0.5)^2
  tmp3 <- scaled_X_cont[,4]
  tmp4 <- scaled_X_cont[,5]
  return(1 * tmp1 + 
           2 * (X_cat_mu[,1] %in% c(0, 2, 4,6,8)) * tmp2 +
           1 * (X_cat_mu[,2] %in% c(1,2,3,6,7,8)) * tmp3 + 
           0.5 * (X_cat_mu[,1] %in% c(1,2,3,4,5)) * tmp4)
}

tau_func <- function(X_cont_tau, X_cat_tau){
  scaled_X_cont <- (X_cont_tau + 1)/2 # moves it to [0,1]
  z1 <- scaled_X_cont[,1]
  z2 <- 1*(X_cat_tau[,1] %in% c(0,1,2,8,9))
  return(3 * z1 + (2 - 5*z2)*sin(pi * z1) - 2 * (z2))
}

n_control <- 20000
n_treat <- 10000
n <- n_control+ n_treat

treated <- c(rep(0, times = n_control), rep(1, times = n_treat))

# controls first and then treated
set.seed(129)
X_cont_mu <- matrix(runif(n * p_cont_mu, -1, 1), nrow = n, ncol = p_cont_mu)
X_cat_mu <- matrix(NA, nrow = n, ncol = p_cat_mu)
cat_levels_list_mu <- list()
for(j in 1:p_cat_mu){
  X_cat_mu[,j] <- as.integer(sample(0:9), size = n, replace = TRUE)
  cat_levels_list_mu[[j]] <- 0:9
}

# only contains stuff for treated subjects
X_cont_tau <- matrix(runif(n_treat * p_cont_tau, -1, 1), nrow = n_treat, ncol = p_cont_tau)
X_cat_tau <- matrix(NA, nrow = n_treat, ncol = p_cat_tau)
cat_levels_list_tau <- list()
for(j in 1:p_cat_mu){
  X_cat_tau[,j] <- as.integer(sample(0:9), size = n_treat, replace = TRUE)
  cat_levels_list_tau[[j]] <- 0:9
}

mu_true <- mu_func(X_cont_mu, X_cat_mu)
tau_true <- tau_func(X_cont_tau, X_cat_tau)

mean_response <- mu_true + c(rep(0, times = n_control), tau_true)

set.seed(328)
sigma <- 0.1
Y <- mean_response + sigma * rnorm(n = n, mean = 0, sd = 1)



fit <- 
  flexBCF(Y_train = Y, 
          treated = treated,
          X_cont_mu = X_cont_mu,
          X_cat_mu = X_cat_mu,
          X_cont_tau = X_cont_tau,
          X_cat_tau = X_cat_tau,
          cat_levels_list_mu = cat_levels_list_mu,
          cat_levels_list_tau = cat_levels_list_tau,
          sparse = TRUE, M_mu = 50, M_tau = 50)

# Get posterior evaluations of mu

mu_samples <-
  get_tree_fits(fit = fit, type = c("mu"),
                X_cont = X_cont_mu,
                X_cat = X_cat_mu)

plot(mu_true, colMeans(mu_samples), pch = 16, cex = 0.2)

tau_samples <-
  get_tree_fits(fit = fit, type = c("tau"),
                X_cont = X_cont_tau,
                X_cat = X_cat_tau)
plot(tau_true, colMeans(tau_samples), pch = 16, cex = 0.2)

}
}

\references{
Hahn, P.R., Murray, J.S., and Carvalho, C.M. (2020). Bayesian regression tree models for causal inference: Regularization, confounding, and heterogeneous effects (with Discussion). 
\emph{Bayesian Analysis}. \bold{15}(3):965--1056. \doi{doi:10.1214/19-BA1195}.

Kokandakar, A.H., Kang, H., and Deshpande, S.K. (2023). Bayesian Causal Forests & the 2022 ACIC Data Challenge: Scalability and Sensitivity. 
\emph{Observational Studies}. \bold{9}(3):29--41. \doi{doi:10.1353/obs.2023.0024}.
}