#ifndef GUARD_update_tree_h
#define GUARD_update_tree_h

#include "funs.h"

void grow_tree_mu(tree &t, suff_stat &ss_train, int &accept, double &sigma, data_info &di_train, tree_prior_info &tree_pi, RNG &gen);
void grow_tree_tau(tree &t, suff_stat &ss_train, int &accept, double &sigma, data_info &di_train, tree_prior_info &tree_pi, RNG &gen);


void prune_tree_mu(tree &t, suff_stat &ss_train, int &accept, double &sigma, data_info &di_train, tree_prior_info &tree_pi, RNG &gen);
void prune_tree_tau(tree &t, suff_stat &ss_train, int &accept, double &sigma, data_info &di_train, tree_prior_info &tree_pi, RNG &gen);

void update_tree_mu(tree &t, suff_stat &ss_train, int &accept, double &sigma, data_info &di_train, tree_prior_info &tree_pi, RNG &gen);
void update_tree_tau(tree &t, suff_stat &ss_train, int &accept, double &sigma, data_info &di_train, tree_prior_info &tree_pi, RNG &gen);



#endif /* update_tree_h */
