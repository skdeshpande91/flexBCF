#include "funs.h"

//

void tree_traversal_mu(suff_stat &ss, tree &t, data_info &di)
{
  double* xx_cont = 0;
  int* xx_cat = 0;
  tree::tree_cp bn;
  int nid;
  ss.clear(); // clear out the sufficient statistic map
  
  tree::npv bnv;
  t.get_bots(bnv);
  suff_stat_it ss_it; // used to look up to which element of ss a particular bottom node corresponds
  
  // add an element to suff stat map for each bottom node
  for(tree::npv_it it = bnv.begin(); it != bnv.end(); ++it){
    nid = (*it)->get_nid(); // bnv is a vector of pointers. it points to elements of bnv so we need (*it) to access the members of these elements
    // add element to sufficient stats map with key = nid, value = empty vector to hold index of observations assigned to this node
    ss.insert(std::pair<int,std::vector<int>>(nid,std::vector<int>()));
  }
  
  for(int i = 0; i < di.n; i++){
    if(di.x_cont_mu != 0) xx_cont = di.x_cont_mu + i * di.p_cont_mu;
    if(di.x_cat_mu != 0) xx_cat = di.x_cat_mu + i * di.p_cat_mu;
    bn = t.get_bn(xx_cont, xx_cat);
    if(bn == 0){
      Rcpp::Rcout << "treated observation i = " << i << std::endl;
      t.print();
      Rcpp::stop("[tree_traversal]: could not find bottom node!");
    }
    else{
      nid = bn->get_nid();
      if(ss.count(nid) != 1) Rcpp::stop("[tree_traversal]: bottom node not included in sufficient statistic map!"); // should never be encountered
      else{
        ss_it = ss.find(nid); // iterator now set to element of ss corresponding to the bottom node holding observation i
        ss_it->second.push_back(i);
      } // closes if/else checking that i-th observation's bottom node was in our map
    } // closes if/else checking that i-th observation maps to valid bottom node
  } // closes loop over all observations
}

void tree_traversal_tau(suff_stat &ss, tree &t, data_info &di)
{
  double* xx_cont = 0;
  int* xx_cat = 0;
  tree::tree_cp bn;
  int nid;
  ss.clear(); // clear out the sufficient statistic map
  
  tree::npv bnv;
  t.get_bots(bnv);
  suff_stat_it ss_it; // used to look up to which element of ss a particular bottom node corresponds
  
  // add an element to suff stat map for each bottom node
  for(tree::npv_it it = bnv.begin(); it != bnv.end(); ++it){
    nid = (*it)->get_nid(); // bnv is a vector of pointers. it points to elements of bnv so we need (*it) to access the members of these elements
    // add element to sufficient stats map with key = nid, value = empty vector to hold index of observations assigned to this node
    ss.insert(std::pair<int,std::vector<int>>(nid,std::vector<int>()));
  }

  // this for a tau tree and we only want it evaluated on the treated subjects
  for(int i = 0; i < di.n_treat; i++){
    if(di.x_cont_tau != 0) xx_cont = di.x_cont_tau + i * di.p_cont_tau;
    if(di.x_cat_tau != 0) xx_cat = di.x_cat_tau + i * di.p_cat_tau;
    bn = t.get_bn(xx_cont, xx_cat);
    if(bn == 0){
      Rcpp::Rcout << "treated observation i = " << i << std::endl;
      t.print();
      Rcpp::stop("[tree_traversal]: could not find bottom node!");
    }
    else{
      nid = bn->get_nid();
      if(ss.count(nid) != 1) Rcpp::stop("[tree_traversal]: bottom node not included in sufficient statistic map!"); // should never be encountered
      else{
        ss_it = ss.find(nid); // iterator now set to element of ss corresponding to the bottom node holding observation i
        ss_it->second.push_back(i);
      } // closes if/else checking that i-th observation's bottom node was in our map
    } // closes if/else checking that i-th observation maps to valid bottom node
  } // closes loop over all observations
}

void fit_ensemble_mu(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di){
  if(fit.size() != di.n) Rcpp::stop("[fit_ensemble]: size of fit must be equal to di.n!"); // honestly should never get triggered
  double* xx_cont = 0;
  int* xx_cat = 0;
  for(int i = 0; i < di.n; i++){
    if(di.x_cont_mu != 0) xx_cont = di.x_cont_mu + i * di.p_cont_mu;
    if(di.x_cat_mu != 0) xx_cat = di.x_cat_mu + i * di.p_cat_mu;
    fit[i] = 0.0;
    for(int m = 0; m < t_vec.size(); m++) fit[i] += t_vec[m].evaluate(xx_cont, xx_cat); // involves a tree traversal
  }
}

void fit_ensemble_tau(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di){
  if(fit.size() != di.n) Rcpp::stop("[fit_ensemble]: size of fit must be equal to di.n!"); // honestly should never get triggered
  double* xx_cont = 0;
  int* xx_cat = 0;
  for(int i = 0; i < di.n; i++){
    if(di.x_cont_tau != 0) xx_cont = di.x_cont_tau + i * di.p_cont_tau;
    if(di.x_cat_tau != 0) xx_cat = di.x_cat_tau + i * di.p_cat_tau;
    fit[i] = 0.0;
    for(int m = 0; m < t_vec.size(); m++) fit[i] += t_vec[m].evaluate(xx_cont, xx_cat); // involves a tree traversal
  }
}


void compute_suff_stat_grow_mu(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, tree &t, data_info &di)
{
  double* xx_cont;
  int* xx_cat;
  int i;
  double tmp_x;
  int l_count;
  int r_count;
  
  // we are growing tree from node nx, which has id of nx_nid
  
  int nxl_nid = 2*nx_nid; // id of proposed left child of nx
  int nxr_nid = 2*nx_nid+1; // id of proposed right child of nx
  
  suff_stat_it nx_it = orig_suff_stat.find(nx_nid); // iterator at element for nx in original sufficient statistic map
  new_suff_stat.clear();
  
  // copy orig_suff_stat into new_suff_stat
  for(suff_stat_it it = orig_suff_stat.begin(); it != orig_suff_stat.end(); ++it){
    new_suff_stat.insert(std::pair<int,std::vector<int>>(it->first, it->second));
  }
  
  // now we manipulate new_suff_stat to drop nx and add nxl and nxr
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxl_nid, std::vector<int>())); // create map element for left child of nx
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxr_nid, std::vector<int>())); // create map element for right child of nx
  new_suff_stat.erase(nx_nid); // remove map element for nx as it is not a bottom leaf node in new tree
  
  suff_stat_it nxl_it = new_suff_stat.find(nxl_nid); // iterator at element for nxl in new sufficient stat map
  suff_stat_it nxr_it = new_suff_stat.find(nxr_nid); // iterator at element for nxr in new sufficient stat map
  
  // loop over all observation that were assigned to nx in original tree
  // note:
  //   nx_it->first is just the node id for nx (nx_nid)
  //   nx_it->second is a vector of integers containing the indicies of observations that land in nx
  // in helper.h we defined int_it as std::vector<int>::iterator
  
  for(int_it it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
    i = *it;
    
    if(di.x_cont_mu != 0) xx_cont = di.x_cont_mu + i * di.p_cont_mu;
    if(di.x_cat_mu != 0) xx_cat = di.x_cat_mu + i * di.p_cat_mu;
    
    if(rule.is_aa && !rule.is_cat){
      // axis-aligned rule
      if(xx_cont[rule.v_aa] < rule.c) nxl_it->second.push_back(i);
      else if(xx_cont[rule.v_aa] >= rule.c) nxr_it->second.push_back(i);
      else Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child");
    } else if(!rule.is_aa && rule.is_cat){
      // categorical rule
      // we need to see whether i-th observation's value of the categorical pred goes to left or right
      // std::set.count returns 1 if the value is in the set and 0 otherwise
      l_count = rule.l_vals.count(xx_cat[rule.v_cat]);
      r_count = rule.r_vals.count(xx_cat[rule.v_cat]);
      if(l_count == 1 && r_count == 0) nxl_it->second.push_back(i);
      else if(l_count == 0 && r_count == 1) nxr_it->second.push_back(i);
      else if(l_count == 1 && r_count == 1) Rcpp::stop("[compute_ss_grow]: observation goes to both left & right child...");
      else Rcpp::stop("[compute_ss_grow]: observation doesn't go to left or right child...");
    } else if(!rule.is_aa && !rule.is_cat){
      // random combination rule
      tmp_x = 0.0;
      for(rc_it rcit = rule.rc_weight.begin(); rcit != rule.rc_weight.end(); ++rcit) tmp_x += (rcit->second) * xx_cont[rcit->first];
      if(tmp_x < rule.c) nxl_it->second.push_back(i);
      else if(tmp_x >= rule.c) nxr_it->second.push_back(i);
      else Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child");
    } else{
      // we should never hit this error
      Rcpp::stop("[compute_ss_grow]: cannot resolve the type of decision rule");
    }
  } // closes loop over all entries in nx
}

void compute_suff_stat_grow_tau(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, tree &t, data_info &di)
{
  double* xx_cont;
  int* xx_cat;
  int i;
  double tmp_x;
  int l_count;
  int r_count;
  
  // we are growing tree from node nx, which has id of nx_nid
  
  int nxl_nid = 2*nx_nid; // id of proposed left child of nx
  int nxr_nid = 2*nx_nid+1; // id of proposed right child of nx
  
  suff_stat_it nx_it = orig_suff_stat.find(nx_nid); // iterator at element for nx in original sufficient statistic map
  new_suff_stat.clear();
  
  // copy orig_suff_stat into new_suff_stat
  for(suff_stat_it it = orig_suff_stat.begin(); it != orig_suff_stat.end(); ++it){
    new_suff_stat.insert(std::pair<int,std::vector<int>>(it->first, it->second));
  }
  
  // now we manipulate new_suff_stat to drop nx and add nxl and nxr
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxl_nid, std::vector<int>())); // create map element for left child of nx
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxr_nid, std::vector<int>())); // create map element for right child of nx
  new_suff_stat.erase(nx_nid); // remove map element for nx as it is not a bottom leaf node in new tree
  
  suff_stat_it nxl_it = new_suff_stat.find(nxl_nid); // iterator at element for nxl in new sufficient stat map
  suff_stat_it nxr_it = new_suff_stat.find(nxr_nid); // iterator at element for nxr in new sufficient stat map
  
  // loop over all observation that were assigned to nx in original tree
  // note:
  //   nx_it->first is just the node id for nx (nx_nid)
  //   nx_it->second is a vector of integers containing the indicies of observations that land in nx
  // in helper.h we defined int_it as std::vector<int>::iterator
  
  for(int_it it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
    i = *it;
    if(di.treated[i + di.n_control] != 1){
      Rcpp::Rcout << "[compute_suff_stat_grow]: trying to update a tree for treatment effect but encountered a control observation" << std::endl;
      Rcpp::Rcout << "  encountered observation " << i + di.n_control << " and treated = " << di.treated[i + di.n_control] << std::endl;
      Rcpp::stop("Trying to update wrong type of tree!");
    } else{
      if(di.x_cont_tau != 0) xx_cont = di.x_cont_tau + i * di.p_cont_tau;
      if(di.x_cat_tau != 0) xx_cat = di.x_cat_tau + i * di.p_cat_tau;
      if(rule.is_aa && !rule.is_cat){
        // axis-aligned rule
        if(xx_cont[rule.v_aa] < rule.c) nxl_it->second.push_back(i);
        else if(xx_cont[rule.v_aa] >= rule.c) nxr_it->second.push_back(i);
        else Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child");
      } else if(!rule.is_aa && rule.is_cat){
        // categorical rule
        // we need to see whether i-th observation's value of the categorical pred goes to left or right
        // std::set.count returns 1 if the value is in the set and 0 otherwise
        l_count = rule.l_vals.count(xx_cat[rule.v_cat]);
        r_count = rule.r_vals.count(xx_cat[rule.v_cat]);
        if(l_count == 1 && r_count == 0) nxl_it->second.push_back(i);
        else if(l_count == 0 && r_count == 1) nxr_it->second.push_back(i);
        else if(l_count == 1 && r_count == 1) Rcpp::stop("[compute_ss_grow]: observation goes to both left & right child...");
        else Rcpp::stop("[compute_ss_grow]: observation doesn't go to left or right child...");
      } else if(!rule.is_aa && !rule.is_cat){
        // random combination rule
        tmp_x = 0.0;
        for(rc_it rcit = rule.rc_weight.begin(); rcit != rule.rc_weight.end(); ++rcit) tmp_x += (rcit->second) * xx_cont[rcit->first];
        if(tmp_x < rule.c) nxl_it->second.push_back(i);
        else if(tmp_x >= rule.c) nxr_it->second.push_back(i);
        else Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child");
      } else{
        // we should never hit this error
        Rcpp::stop("[compute_ss_grow]: cannot resolve the type of decision rule");
      } // closes if/else checking the type of rule
    } // closes if/else checking that the observation is treated
  } // closes loop over all entries in nx
}


void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, tree &t, data_info &di)
{
  //int i;
  if(orig_suff_stat.count(nl_nid) != 1) Rcpp::stop("[compute_ss_prune]: did not find left node in suff stat map");
  if(orig_suff_stat.count(nr_nid) != 1) Rcpp::stop("[compute_ss_prune]: did not find right node in suff stat map");
  
  suff_stat_it nl_it = orig_suff_stat.find(nl_nid); // iterator at element for nl in original suff stat map
  suff_stat_it nr_it = orig_suff_stat.find(nr_nid); // iterator at element for nr in original suff stat map
  
  new_suff_stat.clear();
  // this makes a completely new copy of orig_suff_stat
  for(suff_stat_it ss_it = orig_suff_stat.begin(); ss_it != orig_suff_stat.end(); ++ss_it){
    new_suff_stat.insert(std::pair<int,std::vector<int>>(ss_it->first, ss_it->second));
  }
  new_suff_stat.insert(std::pair<int,std::vector<int>>(np_nid, std::vector<int>())); // add element for np in new suff stat map
  new_suff_stat.erase(nl_nid); // delete element for nl in new suff stat map since nl has been pruned
  new_suff_stat.erase(nr_nid); // delete element for nr in new suff stat map since nr has been pruned
  
  if(new_suff_stat.count(np_nid) != 1) Rcpp::stop("[compute_ss_prune]: didn't create element in new suff stat map for np correctly");
  suff_stat_it np_it = new_suff_stat.find(np_nid); // iterator at element for np in new suff stat map
  
  // time to populate np_it
  // first let's add the elements from nl_it
  for(int_it it = nl_it->second.begin(); it != nl_it->second.end(); ++it) np_it->second.push_back( *it );
  for(int_it it = nr_it->second.begin(); it != nr_it->second.end(); ++it) np_it->second.push_back( *it );
}

double compute_lil(suff_stat &ss, int &nid, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  // reminder posterior of jump mu is N(P^-1 Theta, P^-1)
  if(ss.count(nid) != 1) Rcpp::stop("[compute_lil]: did not find node in suff stat map!");
  suff_stat_it ss_it = ss.find(nid);
  
  double P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
  double Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
  
  for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it]/pow(sigma, 2.0);
  return(-0.5 * log(P) + 0.5 * pow(Theta,2.0) / P);
  
}

double compute_lil_mu(suff_stat &ss, int &nid, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  if(ss.count(nid) != 1) Rcpp::stop("[compute_lil]: did not find node in suff stat map!");
  suff_stat_it ss_it = ss.find(nid);
  
  double P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
  double Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
  
  for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it]/pow(sigma, 2.0);
  return(-0.5 * log(P) + 0.5 * pow(Theta,2.0) / P);
}

double compute_lil_tau(suff_stat &ss, int &nid, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  if(ss.count(nid) != 1) Rcpp::stop("[compute_lil]: did not find node in suff stat map!");
  suff_stat_it ss_it = ss.find(nid);
  
  double P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
  double Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
  // remember that ss_it->second contains the index within the treated set and not the whole dataset
  // we must offset by di.n_control!
  for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it + di.n_control]/pow(sigma, 2.0);
  return(-0.5 * log(P) + 0.5 * pow(Theta,2.0) / P);
}


void draw_mu(tree &t, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  //int i;
  double P;
  double Theta;
  double post_sd;
  double post_mean;
  tree::tree_p bn; // we are modifying bn so we need a pointer not a constant pointer
  
  for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
    bn = t.get_ptr(ss_it->first);
    if(bn == 0) Rcpp::stop("[draw_mu]: could not find node that is in suff stat map in the tree");
    else{
      P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
      Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it]/pow(sigma,2.0);

      post_sd = sqrt(1.0/P);
      post_mean = Theta/P;
      bn->set_mu(gen.normal(post_mean, post_sd));
    }
  }
}

void draw_leaf_mu(tree &t, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  //int i;
  double P;
  double Theta;
  double post_sd;
  double post_mean;
  tree::tree_p bn; // we are modifying bn so we need a pointer not a constant pointer
  
  for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
    bn = t.get_ptr(ss_it->first);
    if(bn == 0) Rcpp::stop("[draw_mu]: could not find node that is in suff stat map in the tree");
    else{
      P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
      Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it]/pow(sigma,2.0);

      post_sd = sqrt(1.0/P);
      post_mean = Theta/P;
      bn->set_mu(gen.normal(post_mean, post_sd));
    }
  }
}

void draw_leaf_tau(tree &t, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  //int i;
  double P;
  double Theta;
  double post_sd;
  double post_mean;
  tree::tree_p bn; // we are modifying bn so we need a pointer not a constant pointer
  
  for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
    bn = t.get_ptr(ss_it->first);
    if(bn == 0) Rcpp::stop("[draw_mu]: could not find node that is in suff stat map in the tree");
    else{
      P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
      Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
      // remember that ss_it->second contains the index within the treated set and not the whole dataset
      // we must offset by di.n_control!
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it + di.n_control]/pow(sigma,2.0);

      post_sd = sqrt(1.0/P);
      post_mean = Theta/P;
      bn->set_mu(gen.normal(post_mean, post_sd));
    }
  }
}

std::string write_tree(tree &t, tree_prior_info &tree_pi, set_str_conversion &set_str)
{
  std::ostringstream os;
  os.precision(32);
  
  tree::cnpv nds;
  rule_t rule;
  t.get_nodes(nds);
  
  for(tree::cnpv_it nd_it = nds.begin(); nd_it != nds.end(); ++nd_it){
    os << (*nd_it)->get_nid() << " ";
    if( (*nd_it)->get_ntype() == 'b'){
      // it's a bottom node
      os << "m " << (*nd_it)->get_mu() << " " << std::endl; // m tells us to expect mu next
    } else if( ((*nd_it)->get_ntype() == 't') && ( !(*nd_it)->l) ){ // because we need to look at left child of a node, make write_tree a friend in tree class
      // it's the top node and it has no children
      os << "m " << (*nd_it)->get_mu() << " " << std::endl; // tree is a stump, m tells us to expect mu next
    } else{
      // we need to print out the rule
      //os << "make rule " << std::endl;
      os << "r "; // node has a decision rule. r tells us to expect a rule next
      rule.clear();
      rule = (*nd_it)->get_rule();
      

      // os << rule.is_cat << " " << rule.is_rc << " ";
      os << rule.is_aa << " " << rule.is_cat << " ";

      if(rule.is_aa && !rule.is_cat){
        // axis-aligned rule
        os << rule.c << " " << rule.v_aa;
      } else if(!rule.is_aa && rule.is_cat){
        // categorical rule
        int K = tree_pi.K->at(rule.v_cat); // how many levels
        os << rule.v_cat << " " << K << " ";
        os << set_str.set_to_hex(K, rule.l_vals) << " ";
        os << set_str.set_to_hex(K, rule.r_vals) << " ";
      } else if(!rule.is_aa && !rule.is_cat){
        // random combination
        os << rule.c << " ";
        for(rc_it rcit = rule.rc_weight.begin(); rcit != rule.rc_weight.end(); ++rcit){
          os << rcit->first << " " << rcit->second << " ";
        }
      } else{
        Rcpp::stop("[write tree]: rule cannot be both axis-aligned and categorical!");
      }
      os << std::endl;
    } // closes if/else checking what type of node we are writing
  } // closes loop over the nodes in t
  
  return os.str();
  
}

void read_tree(tree &t, std::string &tree_string, tree_prior_info &tree_pi, set_str_conversion &set_str)
{
  std::istringstream tree_ss(tree_string); // an in stringstream of the tree's string representation
  std::string node_string; // string for each individual node in the tree
  
  int nid;
  char stream_type; // either 'm' to indicate that next element in stream is mu or 'r' to indicate a rule follows
  
  double tmp_mu; // holds the value of mu for a leaf node

  char aa; // '0' or '1' for rule.is_aa
  char cat; // '0' or '1' for rule.is_cat
  //char rc; // '0' or '1' for rule.is_rc
  rule_t tmp_rule; // temporary rule that gets populated as we read the tree's string/stream
  
  int tmp_v; // for reading in rc weights
  double tmp_phi; // for reading in rc weights
  
  int K; // tells us how many levels there were to the categorical variable
  std::string l_hex; // string representation of the l_vals in a categorical rule
  std::string r_hex; // string representation of the r_vals in a categorical rule
  
  std::map<int, rule_t> decision_nodes;
  std::map<int, double> leaf_nodes;
  
  while(tree_ss){
    std::getline(tree_ss, node_string, '\n');
    if(node_string.size() > 0){
      std::istringstream node_ss(node_string); // in stream for the single node
      node_ss >> nid; // get the node nid
      node_ss >> stream_type;
      
      if(stream_type == 'm'){
        node_ss >> tmp_mu;
        leaf_nodes.insert(std::pair<int,double>(nid, tmp_mu));
      } else if(stream_type == 'r'){
        tmp_rule.clear();
        node_ss >> aa;
        node_ss >> cat;
        
        if(aa == '0') tmp_rule.is_aa = false;
        else tmp_rule.is_aa = true;
        
        if(cat == '0') tmp_rule.is_cat = false;
        else tmp_rule.is_cat = true;
        
        //if(rc == '0') tmp_rule.is_rc = false;
        //else tmp_rule.is_rc = true;
        
        if(tmp_rule.is_aa && !tmp_rule.is_cat){
          // axis-aligned
          node_ss >> tmp_rule.c;
          node_ss >> tmp_rule.v_aa;
        } else if(!tmp_rule.is_aa && tmp_rule.is_cat){
          // categorical rule
          node_ss >> tmp_rule.v_cat; // get the variable index
          node_ss >> K; // we now know how many levels of the categorical variable there were
          node_ss >> l_hex;
          node_ss >> r_hex;
          
          if(K != tree_pi.K->at(tmp_rule.v_cat)){
            Rcpp::Rcout << "v_cat = " << tmp_rule.v_cat << std::endl;
            Rcpp::Rcout << "Read in K = " << K << " total categorical levels" << std::endl;
            Rcpp::Rcout << "Corresponding entry of di.K has length " << tree_pi.K->at(tmp_rule.v_cat) << std::endl;
            Rcpp::stop("mismatch in number of levels recorded in cat_levels & in saved tree!");
          }
          
          tmp_rule.l_vals = set_str.hex_to_set(K, l_hex);
          tmp_rule.r_vals = set_str.hex_to_set(K, r_hex);
        } else if(!tmp_rule.is_aa && !tmp_rule.is_cat){
          // random combination
          node_ss >> tmp_rule.c; // get the cutpoint first
          while(node_ss){
            node_ss >> tmp_v;
            node_ss >> tmp_phi;
            tmp_rule.rc_weight.insert(std::pair<int,double>(tmp_v, tmp_phi));
          }
        } else{
          Rcpp::stop("[read tree]: rule cannot be axis-aligned and categorical");
        }
        decision_nodes.insert(std::pair<int, rule_t>(nid, tmp_rule));
      } // closes if/else checking what type of node we're parsing
    } // closes if checking that we found a valid node in the stream
  } // closes while that parses stream for the tree
  
  // we now have decision_nodes and leaf_nodes and are ready to build up our tree
  t.to_null(); // clear out the tree if there was anything there
  
  // remember std::map is sorted by key.
  // we have always used node id as the key so by iterating over our map, we will *never*
  // attempt to birth from a node that has not already been created.

  for(std::map<int,rule_t>::iterator it = decision_nodes.begin(); it != decision_nodes.end(); ++it){
    t.birth(it->first, it->second); // do the birth.
  }
  
  // since we're messing with private members of tree, do we need to make this function a friend of tree? couldn't hurt.
  tree::npv bnv;
  t.get_bots(bnv); // get the bottom nodes
  std::map<int,double>::iterator leaf_it;

  for(tree::npv_it it = bnv.begin(); it != bnv.end(); ++it){
    leaf_it = leaf_nodes.find( (*it)->get_nid() );
    if(leaf_it == leaf_nodes.end()){
      Rcpp::Rcout << "[read_tree]: we didn't read a leaf node with nid" << (*it)->get_nid() << std::endl;
      Rcpp::stop("mistake in reading in tree");
    } else{
      (*it)->set_mu(leaf_it->second);
    }
  }
  
}



void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen, int &p_cont, int &p_cat){
  rule.clear(); // clear out the rule
  // we now are ready to draw a decision rule

  int rule_counter = 0; // we are allowed multiple tries to draw a valid random combination or categorical rule
  double c_upper = 1.0; // upper bound for range of cutpoints in axis aligned split
  double c_lower = -1.0; // lower bound for range of cutpoints in axis aligned split
  double tmp_weight = 0.0; // weights of random combination
  double c_max = 1.0; // upper bound for absolute value of cutpoint in random combination split
  tree::tree_p nx = t.get_ptr(nid); // at what node are we proposing this rule.
  
  int p = p_cont + p_cat;
  
  double unif = gen.uniform();
  if( (!tree_pi.rc_split) || (tree_pi.rc_split && unif > *tree_pi.prob_rc) ){
    // either we were never allowed to try a rc split OR (1) we are allowed to try rc splits but (2) randomly choose a different type of rule
    
    int v_raw = gen.multinomial(p, tree_pi.theta);
    if(v_raw < p_cont){
      // continuous variable so it's an axis-aligned splitting rule
      rule.is_aa = true;
      rule.is_cat = false;
      rule.v_aa = v_raw;
      if(tree_pi.unif_cuts[rule.v_aa] == 0){
        // draw the cutpoint from the supplied cutpoints
        c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin()); // returns smallest element in set
        c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin()); // reverse iterator, returns largest value in set
        nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
        if(c_lower >= c_upper){
          // this is a weird tree and we'll just propose a trivial split
          c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin());
          c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin());
        }
        std::vector<double> valid_cutpoints;
        if(tree_pi.cutpoints->at(rule.v_aa).count(c_lower) != 1 || tree_pi.cutpoints->at(rule.v_aa).count(c_upper) != 1){
          // c_lower and c_upper were not found in the set of available cutpoints
          Rcpp::Rcout << "[grow tree]: attempting to select a cutpoint from given set" << std::endl;
          Rcpp::Rcout << "  lower bound is: " << c_lower << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_lower) << std::endl;
          Rcpp::Rcout << "  upper bound is: " << c_upper << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_upper) << std::endl;
          //Rcpp::Rcout << "  cutpoints are:";
          //for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it) Rcpp::Rcout << " " << *it;
          //Rcpp::Rcout << std::endl;
          Rcpp::stop("we should never have a c that is outside the pre-defined set of cutpoints!");
        }
        // we want to draw from the cutpoints exclusive of c_lower & c_upper;
        // i.e. we want to start with the one just after c_lower and just before c_upper
        // std::set::lower_bound: iterator at first element that is not considered to come before
        // std::set::upper_bound: iterator at first element considered to come after
        // if value is not in set, lower_bound and upper_bound give same result
        // if value is in set: lower bound returns the value, upper bound returns the next value
        for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).upper_bound(c_lower); it != tree_pi.cutpoints->at(rule.v_aa).lower_bound(c_upper); ++it){
          valid_cutpoints.push_back(*it);
        }
        int num_cutpoints = valid_cutpoints.size();
        if(num_cutpoints < 1){
          // no valid splits are available; we will just pick something, all of the observations will go to one child anyway...
          valid_cutpoints.clear();
          for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it){
            valid_cutpoints.push_back(*it);
          }
          num_cutpoints = valid_cutpoints.size();
        }
        // at this point, valid cutpoints is a vector containing the available cutpoints at this node. we pick one uniformly.
        rule.c = valid_cutpoints[floor(gen.uniform() * num_cutpoints)];
        
      } else{
        // draw cutpoints uniformly
        c_upper = 1.0;
        c_lower = -1.0;
        nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
        if(c_lower >= c_upper){
          c_lower = -1.0;
          c_upper = 1.0;
        }
        rule.c = gen.uniform(c_lower, c_upper);
      }
    } else{
      // categorical decision rule
      rule.is_aa = false;
      rule.is_cat = true;
      rule.v_cat = v_raw - p_cont;

      std::set<int> avail_levels = tree_pi.cat_levels->at(rule.v_cat); // get the full set of levels for this variable
      nx->get_rg_cat(rule.v_cat, avail_levels); // determine the set of levels available at nx.
      // if there is only one level left for this variable at nx, we will just propose a trivial split
      // and will reset the value of avail_levels to be the full set of all levels for the variable
      if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat);
      
      rule.l_vals.clear();
      rule.r_vals.clear();
      
      if(tree_pi.graph_split[rule.v_cat] == 1 && (tree_pi.adj_support->at(rule.v_cat).size() > 0)){
        // if we explicitly say to use the MST to split the variables
        graph_partition(avail_levels, rule.l_vals, rule.r_vals, tree_pi.adj_support->at(rule.v_cat), tree_pi.K->at(rule.v_cat), tree_pi.graph_cut_type, gen);
      } else{
        // otherwise we default to splitting the available levels uniformly at random: prob 0.5 to go to each child
        rule_counter = 0;
        while( ((rule.l_vals.size() == 0) || (rule.r_vals.size() == 0)) && rule_counter < 1000 ){
          rule.l_vals.clear();
          rule.r_vals.clear();
          for(set_it it = avail_levels.begin(); it != avail_levels.end(); ++it){
            if(gen.uniform() <= 0.5) rule.l_vals.insert(*it);
            else rule.r_vals.insert(*it);
          }
          ++(rule_counter);
        }
        if(rule_counter == 1000){
          Rcpp::stop("failed to generate valid categorical split in 1000 attempts"); // this should almost surely not get triggered.
        }
      }
      if( (rule.l_vals.size() == 0) || (rule.r_vals.size() == 0) ) Rcpp::stop("proposed an invalid categorical rule!");
    } // closes if/else determining whether we do an axis-aligned or categorical decision rule
  } else if(tree_pi.rc_split && unif <= *tree_pi.prob_rc){
    // random combination rule
    rule.is_aa = false;
    rule.is_cat = false;

    while( (rule.rc_weight.size() < 2) && (rule_counter < 1000) ){
      rule.rc_weight.clear();
      c_max = 0.0;
      for(int j = 0; j < p_cont; j++){
        if(gen.uniform() < (*tree_pi.theta_rc)){
          tmp_weight = gen.uniform(-1.0,1.0); // Breiman used Uniform(-1,1) weights and so shall we
          rule.rc_weight.insert(std::pair<int,double>(j,tmp_weight));
          c_max += fabs(tmp_weight);
        }
      }
      ++(rule_counter);
    }
    if(rule.rc_weight.size() < 2) Rcpp::stop("[propose_rule]: failed to generate a valid random combination rule in 1000 attempts!");
    else{
      rule.c = gen.uniform(-1.0,1.0) * c_max;
    }
  } else{
    Rcpp::stop("[draw_rule]: unable to draw a rule. check values of rc_split & prob_rc");
  } // closes all if/else determining the type of rule (axis-aligned or categorical) OR random combination
}


// depth-first search, used to find connected components of a graph
void dfs(int i, std::vector<bool> &visited, std::vector<int> &comp, int &n, arma::mat &A)
{
  visited[i] = true; // dfs has reached i for the first time so mark it
  comp.push_back(i); // now that i has been marked, add it to the connected component
  for(int ii = 0; ii < n; ii++){
    if( std::fabs(A(i,ii)) > 1e-16 ){ // in case A is a weighted matrix
      // i is connected to ii
      if(!visited[ii]){
        // somehow ii hasn't been visited before, so we need to continue our dfs from there.
        dfs(ii, visited, comp, n, A);
      }
    }
  }
}

// find connected components of a graph
void find_components(std::vector<std::vector<int> > &components, arma::mat &A)
{
  components.clear(); // clear it out
  int n = A.n_rows;
  std::vector<bool> visited(n,false);
  for(int i = 0; i < n; i++){
    if(!visited[i]){
      std::vector<int> new_comp;
      dfs(i, visited, new_comp, n, A);
      components.push_back(new_comp);
    }
  }
}

// find minimum edge weight in Boruvka's algorithm
std::pair<int,int> find_min_edge_weight(std::vector<int> &components, int &n, arma::mat &W)
{
  if(components.size() == n) Rcpp::stop("[add_min_weight_edge]: component already contains all n nodes!");
  
  std::vector<unsigned int> tmp_in;
  std::vector<unsigned int> tmp_out;
  
  std::vector<int>::iterator find_it;
  for(int i = 0; i < n; i++){
    find_it = std::find(components.begin(), components.end(), i);
    if(find_it != components.end()) tmp_in.push_back( (unsigned int) i);
    else tmp_out.push_back( (unsigned int) i);
  }
  arma::uvec in_index(tmp_in);
  arma::uvec out_index(tmp_out);
  
  arma::mat tmp_W = W(in_index, out_index);
  // check that there are actually edges from component to rest of the graph
  if(arma::all(arma::abs(arma::vectorise(tmp_W)) < 1e-16)) Rcpp::stop("tmpW contains all 0's");
  tmp_W.elem(arma::find(arma::abs(tmp_W) < 1e-16)).ones(); // just to be extra safe, convert all 0's into 1's
  
  
  arma::uword min_index = tmp_W.index_min(); //
  arma::uvec min_sub = arma::ind2sub(arma::size(tmp_W), min_index);
  std::pair<int,int> results( in_index(min_sub(0)), out_index(min_sub(1)) );
  return results;
}
// implement's Boruvka's algorithm
arma::mat boruvka(arma::mat &W)
{
  int n = W.n_rows;
  
  std::vector<std::vector<int> > W_components;
  find_components(W_components, W);
  if(W_components.size() > 1){
    Rcpp::stop("W is not connected!");
  }
  
  arma::mat A_mst = arma::zeros<arma::mat>(n,n); // initialize the
  std::vector<std::vector<int> > mst_components;
  find_components(mst_components, A_mst);
  
  int counter = 0; // failsafe which should *never* be triggered since graph is connected.
  // Boruvka's requires O(log(n)) steps so we have intentionally set the upper bound on the number of iterations much higher.
  while( (mst_components.size() > 1) && (counter < n) ){
    for(std::vector<std::vector<int> >::iterator it = mst_components.begin(); it != mst_components.end(); ++it){
      
      std::pair<int,int> min_edge = find_min_edge_weight(*it,n, W);
      A_mst(min_edge.first, min_edge.second) = 1;
      A_mst(min_edge.second, min_edge.first) = 1;
    }
    find_components(mst_components, A_mst);
    ++counter;
  }
  
  if(mst_components.size() > 1){
    Rcpp::Rcout << "Was not able to find a spanning tree (check connectivity of W!)" << std::endl;
  }
  return A_mst;
}


arma::uword get_cut_edge(const arma::mat &cut_A, const arma::mat &cut_W, const arma::uvec mst_edge_index, const int &graph_cut_type, RNG &gen)
{
  // cut_A is the lower triangular adjacency matrix of the MST
  // cut_W is the lower triangular weighted adjacency matrix of the MST
  // we are going to output a single index which will be cut
  
  int n = cut_A.n_rows;
  
  // n is the number of nodes
  // there should only be n-1 edges in the MST
  std::vector<double> cut_ix_probs(n-1, 1.0/( (double) (n-1) ));
  double tmp_sum;
  arma::uword cut_edge_index;
  if(graph_cut_type == 0){
    // delete an edge uniformly at random
    cut_edge_index = mst_edge_index(gen.multinomial(n-1, cut_ix_probs));
  } else if(graph_cut_type == 1){
    // delete edge w.p. proportional to the weight of the
    tmp_sum = 0.0;
    for(int e_ix = 0; e_ix < n-1; e_ix++){
      cut_ix_probs[e_ix] = cut_W(mst_edge_index(e_ix));
      tmp_sum += cut_ix_probs[e_ix];
    }
    for(int e_ix = 0; e_ix < n-1; e_ix++) cut_ix_probs[e_ix] /= tmp_sum;
    cut_edge_index = mst_edge_index(gen.multinomial(n-1, cut_ix_probs));
    
  } else if(graph_cut_type == 2){
    // size-biased
    tmp_sum = 0.0;
    arma::mat tmp_cut_A = cut_A;
    for(int e_ix = 0; e_ix < n-1; e_ix++){
      tmp_cut_A = cut_A;
      if(tmp_cut_A(mst_edge_index(e_ix)) == 0) Rcpp::stop("trying to delete an edge not in the MST!");
      else{
        tmp_cut_A(mst_edge_index(e_ix)) = 0; // delete the edge!
        arma::mat part_mst_A = arma::symmatl(tmp_cut_A); // adjacency of the partitioned tree
        std::vector<std::vector<int>> cut_components;
        find_components(cut_components, part_mst_A);
        cut_ix_probs[e_ix] = (double) std::min(cut_components[0].size(), cut_components[1].size());
        tmp_sum += cut_ix_probs[e_ix];
      }
    }
    for(int e_ix = 0; e_ix < n-1; e_ix++) cut_ix_probs[e_ix] /= tmp_sum;
    cut_edge_index = mst_edge_index(gen.multinomial(n-1, cut_ix_probs));
  } else if(graph_cut_type == 3){
    // delete the edge with largest weight
    cut_edge_index = cut_W.index_max();
  } else{
    Rcpp::stop("[get_cut_edge]: don't know how to cut an edge from the MST!");
  }
  return cut_edge_index;
}





/*
void get_edge_probs(std::vector<double> &cut_ix_probs, const arma::mat &cut_A, const arma::uvec &mst_index, const int &n)
{
  arma::mat tmp_cut_A = cut_A;
  cut_ix_probs.clear();
  double tmp_sum = 0.0;
  for(int cut_ix = 0; cut_ix < n-1; cut_ix++){
    tmp_cut_A = cut_A;
    arma::uvec cut_sub = arma::ind2sub(arma::size(tmp_cut_A), mst_index(cut_ix));
    if(tmp_cut_A(cut_sub(0), cut_sub(1)) == 0){
      Rcpp::stop("[get_min_component_size]: trying to remove an edge that doesn't exist...");
    } else{
      tmp_cut_A(cut_sub(0), cut_sub(1)) = 0; // delete the edge!
      arma::mat part_mst_A = arma::symmatl(tmp_cut_A); // adjacency matrix of the partitioned tree
      std::vector<std::vector<int> > cut_components;
      find_components(cut_components, part_mst_A);
      cut_ix_probs.push_back( (double) std::min(cut_components[0].size(), cut_components[1].size()));
      tmp_sum += (double) std::min(cut_components[0].size(), cut_components[1].size());
    }
  }
  
  for(int cut_ix = 0; cut_ix < n-1; cut_ix++) cut_ix_probs[cut_ix] /= tmp_sum;
}
*/

// adj_support is only for the lower triangle of the adjacency matrix
void graph_partition(std::set<int> &vals, std::set<int> &l_vals, std::set<int> &r_vals, std::vector<unsigned int> &adj_support, int &K, int &graph_cut_type, RNG &gen)
{
  arma::mat W = arma::zeros<arma::mat>(K,K); // we need to create a weighted adjacency matrix
  for(std::vector<unsigned int>::iterator w_it = adj_support.begin(); w_it != adj_support.end(); ++w_it) W(*w_it) = gen.uniform();
  // note that at this point W is lower triangular
  W = arma::symmatl(W); // make W symmetric
  
  // we need to subset  W to just the rows & columns corresponding to vals
  int n = vals.size(); // how many levels of the variable are there = how many vertices in the induced subgraph
  if(n == 1) Rcpp::stop("[graph_partition]: vals contains only one 1 value; cannot partition it further!");
    
  arma::uword cut_edge_index;
  
  // need to subset A and W to just the subgraph induced by the vertices corresponding to the levels in vals
  std::vector<unsigned int> tmp_in;
  for(set_it it = vals.begin(); it != vals.end(); ++it) tmp_in.push_back( (unsigned int) *it);
  arma::uvec index(tmp_in); // arma is weird with indexing with vectors, need to use a uvec
  arma::mat tmp_W = W(index,index);
  
  std::vector<std::vector<int> > components;
  find_components(components, tmp_W); // how many connected components are in W?
  if(components.size() == 0){
    Rcpp::stop("subgraph induced by current set of levels appears to have no connected components...");
  } else if(components.size() > 1){
    // there are multiple connected components
    // for now: just randomly assign whole components to the left or right uniformly at random
    l_vals.clear();
    r_vals.clear();
    int rule_counter = 0;
    while( ((l_vals.size() == 0 || r_vals.size() == 0)) && rule_counter < 1000 ){
      l_vals.clear();
      r_vals.clear();
      for(int comp_ix = 0; comp_ix < components.size(); comp_ix++){
        if(gen.uniform() <= 0.5){
          // send everything in this component to the left child
          for(int_it it = components[comp_ix].begin(); it != components[comp_ix].end(); ++it) l_vals.insert( (int) tmp_in[*it]);
        } else{
          // send everything in this component to the right child
          for(int_it it = components[comp_ix].begin(); it != components[comp_ix].end(); ++it) r_vals.insert( (int) tmp_in[*it]);
        }
      }
      ++rule_counter;
    }
    if(rule_counter == 1000){
      Rcpp::stop("[graph partition]: graph disconnected. failed to generate a non-trivial partiton of components in 1000 attempts!");
    }
  } else{
    // now that we know levels induce a connected subgraph of the original graph
    // we can run Boruvka's algorithm
    arma::mat A_mst = boruvka(tmp_W); // adjacency matrix of the MST
    arma::mat W_mst = A_mst % tmp_W; // what are the weights
    arma::mat cut_A = arma::trimatl(A_mst); // get the lower triangle of the adjacency matrix of the MST
    arma::mat cut_W = arma::trimatl(W_mst); // get the lower triangle of the weighted adjacency matrix
    arma::uvec mst_edge_index = arma::find(cut_A); // indices of edges in cut_A
    
    cut_edge_index = get_cut_edge(cut_A, cut_W, mst_edge_index, graph_cut_type, gen);
    
    if(cut_A(cut_edge_index) != 1){
      //Rcpp::Rcout << "mst_edge_index = " << std::endl;
      //mst_edge_index.t().print();
      //Rcpp::Rcout << "cut_edge_index = " << cut_edge_index << std::endl;
      //Rcpp::Rcout << cut_A(cut_edge_index) << std::endl;
      Rcpp::stop("[graph_partition]: attempting to delete an edge which does not seem to exist!");
    
    }
    else{
      cut_A(cut_edge_index) = 0; // actually delete the edge
      arma::mat part_mst_A = arma::symmatl(cut_A); // full adjacency matrix of the partitioned tree
      std::vector<std::vector<int>> cut_components;
      find_components(cut_components, part_mst_A);
      if(cut_components.size() != 2){
        // this should never get thrown
        Rcpp::stop("[graph_partition]: we have more than 2 connected components after deleting a single edge from a MST!");
      } else{
        l_vals.clear();
        r_vals.clear();
        // *it is an integer between 0 & n and indexes the temporary edge labels
        // the i-th element of tmp_in corresponds to the i-th edge label
        // we need to look up the correct level (i.e. value in vals) using tmp_in
        for(int_it it = cut_components[0].begin(); it != cut_components[0].end(); ++it) l_vals.insert( (int) tmp_in[*it]);
        for(int_it it = cut_components[1].begin(); it != cut_components[1].end(); ++it) r_vals.insert( (int) tmp_in[*it]);
      }
    } // closes if/else checking that the edge we are trying to delete actually exists
  } // closes if/else checking that levels in vals induced a connected subgraph of the original adjancecy matrix
}

void update_theta_u(std::vector<double> &theta, double &u, std::vector<int> &var_count, int &p, double &a_u, double &b_u, RNG &gen)
{
  if(theta.size() != p){
    Rcpp::Rcout << "theta has size " << theta.size() << "  p = " << p << std::endl;
    Rcpp::stop("theta must have size p!");
  } else{
    double tmp_sum = 0.0;
    double tmp_concentration = 0.0;
    double sum_log_theta = 0.0;
    int v_count;
    std::vector<double> tmp_gamma(p, 0.0);
    
    // update theta first
    double u_orig = u;
    for(int j = 0; j < p; j++){
      v_count = var_count[j];
      tmp_concentration = u_orig/(1.0 - u_orig) + (double) v_count;
      tmp_gamma[j] = gen.gamma(tmp_concentration, 1.0);
      tmp_sum += tmp_gamma[j];
    }
    for(int j = 0; j < p; j++){
      theta[j] = tmp_gamma[j]/tmp_sum;
      sum_log_theta += log(theta[j]);
    }
    
    // we're now ready to update u
    double u_prop = gen.beta(a_u,b_u);
    double log_like_prop = (u_prop)/(1.0 - u_prop) * sum_log_theta;
    double log_like_orig = (u_orig)/(1.0 - u_orig) * sum_log_theta;
    
    log_like_prop += lgamma( (double) p * u_prop/(1.0 - u_prop)) - ((double) p) * lgamma(u_prop/(1.0 - u_prop));
    log_like_orig += lgamma( (double) p * u_orig/(1.0 - u_orig)) - ((double) p) * lgamma(u_orig/(1.0 - u_orig));
    double log_accept = log_like_prop - log_like_orig;
    if(gen.log_uniform() <= log_accept) u = u_prop;
    else u = u_orig;
  }
}
void update_theta_rc(double& theta_rc, int &rc_var_count, int &rc_rule_count, double &a_rc, double &b_rc, int &p_cont, RNG &gen)
{
  // since we disallow rc rules with 0 or 1 variable, the posterior is proportional to
  // theta^(a + rc_var_count - 1) * (1 - theta)^(b + p_cont * rc_rule_count - rc_var_count - 1)/(1 - (1- theta)^p_cont - p_cont * theta * 1 - theta)^(p_cont - 1))
  // note that the denominator here is just P(Bin(p_cont, theta) >= 2) as a function of theta
  
  // use independence MH: transition proposal is Beta(a  + rc_var_count, b + p_cont * rc_rule_count)
  // acceptance ratio turns out to:
  // P( Bin(p_cont, theta_orig) >= 2) / P(Bin(p_cont, theta_prop))
  
  double a_post = a_rc + (double) rc_var_count;
  double b_post = b_rc + ((double) p_cont) * rc_rule_count - (double) rc_var_count;
  
  double theta_orig = theta_rc;
  double theta_prop = gen.beta(a_post, b_post);
  
  double log_post_prop = log(1.0 - pow(1.0 - theta_prop, p_cont) - (double) p_cont * theta_prop * pow(1.0 - theta_prop, p_cont-1));
  double log_post_orig = log(1.0 - pow(1.0 - theta_orig, p_cont) - (double) p_cont * theta_orig * pow(1.0 - theta_orig, p_cont-1));
 
  double log_alpha = log_post_orig - log_post_prop;
  if(gen.log_uniform() < log_alpha) theta_rc = theta_prop;
  else theta_rc = theta_orig;
}

/* eventually we will add this functionality back in.


void update_theta_u_cat(std::vector<double> &theta_cat, std::vector<int> &cat_var_count, double &u_cat, double& a_cat, double& b_cat, int &p_cat, RNG &gen)
{
  // stuff for updating theta
  double tmp_sum = 0.0;
  int v_count = 0;
  double tmp_concentration = 0.0;
  std::vector<double> tmp_gamma(p_cat);
  
  // stuff for updating u
  double u_prop = 0.0;
  double u_orig = 0.0;
  double sum_log_theta = 0.0;
  double log_like_prop = 0.0;
  double log_like_orig = 0.0;
  double log_accept = 0.0;
  
  u_orig = u_cat;
  for(int j = 0; j < p_cat; j++){
    v_count = cat_var_count[j];
    tmp_concentration = u_orig/(1.0 - u_orig) + (double) v_count;
    tmp_gamma[j] = gen.gamma(tmp_concentration, 1.0);
    tmp_sum += tmp_gamma[j];
  }
  for(int j = 0; j < p_cat; j++){
    theta_cat[j] = tmp_gamma[j]/tmp_sum;
    sum_log_theta += log(theta_cat[j]);
  }
  
  u_prop = gen.beta(a_cat, b_cat);
  log_like_prop = (u_prop)/(1.0 - u_prop) * sum_log_theta;
  log_like_prop += lgamma( ((double) p_cat) * u_prop/(1.0 - u_prop)) - ((double) p_cat) * lgamma(u_prop/(1.0 - u_prop));
  
  log_like_orig = (u_orig)/(1.0 - u_orig) * sum_log_theta;
  log_like_orig += lgamma( ((double) p_cat) * u_orig/(1.0 - u_orig)) - ( (double) p_cat) * lgamma(u_orig/(1.0 - u_orig));
  
  log_accept = log_like_prop - log_like_orig;
  if(log_accept >= 0.0) log_accept = 0.0;
  if(gen.log_uniform() <= log_accept) u_cat = u_prop;
  else u_cat = u_orig;
}
 
 
 
 void update_theta_cont(std::vector<double> &theta_cont, std::vector<int> &cont_var_count, int &cont_rule_count, double &a_cont, double &b_cont, int &p_cont, RNG &gen)
 {
   if(theta_cont.size() != p_cont) Rcpp::stop("[update_theta_cont]: theta_cont must have size p_cont");
   double a_post = a_cont;
   double b_post = b_cont;
   for(int j = 0; j < p_cont; j++){
     a_post = a_cont + (double) cont_var_count[j];
     b_post = b_cont + (double)(cont_rule_count - cont_var_count[j]);
     theta_cont[j] = gen.beta(a_post, b_post);
   }
   
 }
*/


