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
                    sparse = TRUE, 
                    M_mu = 200, M_tau = 200,
                    nd = 1000, burn = 1000, thin = 1,
                    verbose = TRUE, print_every = floor( (nd*thin + burn))/10)
{
  
  
  
  
}