\name{flexBCF}
\alias{flexBCF}

\title{
A faster and more flexible Bayesian Causal Forests
}
\description{
Implements a version of Hahn et al. (2020)'s Bayesian Causal Forest (BCF) model that is faster and handles categorical outcomes more flexibly than existing BCF implementations. This implementation was created for the 2022 American Causal Inference Conference's Data Challenge.
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
        M_mu = 200, M_tau = 200,
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
  \item{M_mu}{Number of trees in the ensemble for prognostic function mu. Default is 200.}
  \item{M_tau}{Number of trees in the ensemble for treatment effect function tau. Default is 200.}
  \item{nd}{Number of posterior draws to return. Default is 1000.}
  \item{burn}{Number of MCMC iterations to be treated as "warmup" or "burn-in". Default is 1000.}
  \item{thin}{Number of post-warmup MCMC iteration by which to thin. Default is 1.}
  \item{save_samples}{Logical, indicating whether to return all posterior samples. Default is \code{TRUE}. If \code{FALSE}, only posterior mean is returned.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.}
}
\details{
TO DO: more details about our implementation


Default implementations of Bayesian Additive Regression Trees (BART) represent categorical predictors using several binary indicators, one for each level of each categorical predictor. Axis-aligned decision rules are well-defined with these indicators; they send one level of a categorical predictor to the left and all other levels to the right (or vice versa). Regression trees built with these rules partition the set of all levels of a categorical predictor by recursively removing one level at a time. Unfortunately, most partitions of the levels cannot be built with this ``remove one at a time'' strategy, meaning that default implementations of BART are extremely limited in their ability to ``borrow strength'' across groups of levels.

\code{flexBCF} overcomes this limitation using a new prior on regression trees. Under this new prior, conditional on splitting on a categorical predictor at a particular node in the tree, levels of the predictor are sent to the left and right child uniformly at random. In this way, multiple levels of a categorical predictor are able to be clustered together.
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
\examples{
\dontrun{
1+1
}
}