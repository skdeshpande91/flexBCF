get_tree_fits <- function(fit, 
                          type = c("mu","tau"),
                          X_cont = matrix(0, nrow = 1, ncol = 1),
                          X_cat = matrix(0, nrow = 1, ncol = 1),
                          verbose = TRUE, 
                          print_every = floor(max(c(nrow(X_cont), nrow(X_cat)))/10))
{
  
  n <- max(c(nrow(X_cont), nrow(X_cat)))
  if(n == 1){
    # Things go a bit haywire if we try to make a prediction with only one subject at a time
    # Best to copy the individual row and make predictions twice. 
    # Will support this case later
    stop("Computing average effects with n = 1 subject is not currently supported")
  } else{
    if(type == "mu"){
      tmp <- .predict_tree_ensemble(tree_draws = fit$mu_trees,
                                    tX_cont = t(X_cont),
                                    tX_cat = t(X_cat),
                                    treat = FALSE, 
                                    y_mean = fit$y_mean,
                                    y_sd = fit$y_sd,
                                    cat_levels_list = fit$cat_levels_list[["mu"]],
                                    verbose = verbose, print_every = print_every)
    } else if(type == "tau"){
      tmp <- .predict_tree_ensemble(tree_draws = fit$tau_trees,
                                    tX_cont = t(X_cont),
                                    tX_cat = t(X_cat),
                                    treat = TRUE, 
                                    y_mean = fit$y_mean,
                                    y_sd = fit$y_sd,
                                    cat_levels_list = fit$cat_levels_list[["tau"]],
                                    verbose = verbose, print_every = print_every)
    } else{
      stop("type must be one of mu or tau")
    }
    return(tmp)
  }

}