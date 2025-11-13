flexBCF <- function(Y_train,
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
{
  
  # Standardize the Y's
  y_mean <- mean(Y_train)
  y_sd <- stats::sd(Y_train)
  std_Y <- (Y_train - y_mean)/y_sd
  nu <- 3
  lambda <- stats::qchisq(0.1, df = nu)/nu
  mu0 <- c(0,0)
  tau <- c(1/sqrt(M_mu), 1/sqrt(M_tau))
  
  graph_split_mu <- rep(FALSE, times = ncol(X_cat_mu))
  graph_split_tau <- rep(FALSE, times = ncol(X_cat_tau))
  adj_support_list_mu <- NULL
  adj_support_list_tau <- NULL
  
  fit <- .flexBCF(Y_train = std_Y,
                  treated = treated,
                  tX_cont_mu_train = t(X_cont_mu),
                  tX_cat_mu_train = t(X_cat_mu),
                  tX_cont_tau_train = t(X_cont_tau),
                  tX_cat_tau_train = t(X_cat_tau),
                  unif_cuts_mu = unif_cuts_mu,
                  unif_cuts_tau = unif_cuts_tau,
                  cutpoints_list_mu = cutpoints_list_mu,
                  cutpoints_list_tau = cutpoints_list_tau,
                  cat_levels_list_mu = cat_levels_list_mu,
                  cat_levels_list_tau = cat_levels_list_tau,
                  graph_split_mu = graph_split_mu,
                  graph_split_tau = graph_split_tau,
                  graph_cut_type_mu = 0, graph_cut_type_tau = 0,
                  adj_support_list_mu = adj_support_list_mu,
                  adj_support_list_tau = adj_support_list_tau,
                  sparse = sparse, a_u = 1, b_u = 1,
                  mu0 = mu0, tau = tau,
                  lambda = lambda, nu = nu,
                  M_mu = M_mu, M_tau = M_tau,
                  nd = nd, burn = burn, thin = thin,
                  verbose = verbose, print_every = print_every)
  results <- list()
  results[["sigma"]] <- fit$sigma
  results[["mu_trees"]] <- fit$mu
  results[["tau_trees"]] <- fit$tau
  results[["varcount_mu"]] <- fit$varcount_mu
  results[["varcount_tau"]] <- fit$varcount_tau
  results[["y_mean"]] <- y_mean
  results[["y_sd"]] <- y_sd
  results[["cat_levels_list"]] <- list(mu = cat_levels_list_mu, tau = cat_levels_list_tau)
  return(results)
}