#include "funs.h"

// [[Rcpp::export(".predict_tree_ensemble")]]
arma::mat predict_flexBART(Rcpp::List tree_draws,
                           Rcpp::NumericMatrix tX_cont,
                           Rcpp::IntegerMatrix tX_cat,
                           bool treat,
                           double y_mean,
                           double y_sd,
                           Rcpp::Nullable<Rcpp::List> cat_levels_list,
                           bool verbose = true, int print_every = 50)
{

  set_str_conversion set_str;
  
  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  
  parse_training_data(n, p_cont, p_cat, tX_cont, tX_cat);
  int p = p_cont + p_cat;
  
  
  
  // actually we need to look at categorical levels list to get K.
  // maybe we can just read in K?
  
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K;
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      parse_cat_levels(cat_levels, K, p_cat, tmp_cat_levels);
    }
  }
  
  // setup data info object
  data_info di;

  di.n = n;
  if(treat){
    di.p_cont_tau = p_cont;
    di.p_cat_tau = p_cat;
    if(p_cont > 0) di.x_cont_tau = tX_cont.begin();
    if(p_cat > 0)  di.x_cat_tau = tX_cat.begin();
    di.p_tau = p;
  } else{
    di.p_cont_mu = p_cont;
    di.p_cat_mu = p_cat;
    if(p_cont > 0) di.x_cont_mu = tX_cont.begin();
    if(p_cat > 0)  di.x_cat_mu = tX_cat.begin();
    di.p_mu = p;
  }
  
  tree_prior_info tree_pi;
  tree_pi.cat_levels = &cat_levels;
  tree_pi.K = &K;

  int nd = tree_draws.size();
  //int M = tree_draws[0].size();
  
  Rcpp::CharacterVector first_tree_vec = tree_draws[0];
  int M = first_tree_vec.size();

  Rcpp::Rcout << "nd = " << nd << "M = " << M;
  Rcpp::Rcout << " n = " << n << " p_cont = " << p_cont << " p_cat = " << p_cat << std::endl;
  
  std::vector<double> allfit(n);
  //arma::mat pred_out(n, nd);
  arma::mat pred_out(nd,n);

  for(int iter = 0; iter < nd; iter++){
    if( (iter%print_every == 0)){
      Rcpp::Rcout << "  Iteration: " << iter << " of " << nd <<std::endl;
      Rcpp::checkUserInterrupt();
    }
    Rcpp::CharacterVector tmp_string_vec = tree_draws[iter];
    if(tmp_string_vec.size() != M){
      // did we somehow not record enough tree strings?
      // this should really never be hit
      // essentially we're mimicing the R code all(sapply(tree_draws, FUN = length) = length(tree_draws[1]))
      Rcpp::Rcout << "iter = " << iter << " # tree strings = " << tmp_string_vec.size() << std::endl;
      Rcpp::stop("Unexpected number of tree strings!");
    } else{
      std::vector<tree> t_vec(M);
      for(int m = 0; m < M; m++){
        // tmp_string_vec is an Rcpp::CharacterVector
        // let's extract a single element from the CharacterVector and turn it into a std::string
        // that can be passed to read_tree
        std::string tmp_string = Rcpp::as<std::string>(tmp_string_vec[m]); // convert content of the
        read_tree(t_vec[m], tmp_string, tree_pi, set_str);
      }
      if(treat) fit_ensemble_tau(allfit, t_vec, di);
      else fit_ensemble_mu(allfit, t_vec, di);
      //for(int i = 0; i < n; i++) pred_out(i,iter) = allfit[i];
      for(int i = 0; i < n; i++) pred_out(iter,i) = allfit[i];

    } // closes if/else checking that we have M strings for the draw of the ensemble
  } // closes loop over all draws of the ensemble
 
  if(treat){
    // multiply the tau function by y_sd (to undo the standardization)
    pred_out *= y_sd;
  } else{
    // multiply the mu function by y_sd AND add back y_mean (to undo the standardization)
    pred_out *= y_sd;
    pred_out += y_mean;
  }
  
  return pred_out;
}
