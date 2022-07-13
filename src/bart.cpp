#define _USE_MATH_DEFINES
#include <Rcpp.h>
#include<RcppEigen.h>
#include <vector>
#include<cmath>
#include <math.h>
#include "aux_functions.h"
#include "tree.h"

// [[Rcpp::depends(RcppEigen)]]
// The line above (depends) it will make all the dependcies be included on the file
using namespace Rcpp;
using namespace RcppEigen;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Upper;


// Nan error checker
// void contain_nan(Eigen::VectorXd vec){
//
//   for(int i = 0; i<vec.size();i++){
//     if(isnan(vec.coeff(i))==1){
//       cout << "This vector contain nan" << endl;
//     }
//   }
//
// }




// [[Rcpp::depends(RcppEigen)]]
Tree grow(Tree curr_tree,
          Eigen::MatrixXd x,
          int n_min_size){

  // Declaring variables
  bool bad_tree=true;
  int bad_tree_counter = 0;

  // Getting information from x
  int p = x.cols();
  int n_terminals;
  int n_nodes = curr_tree.list_node.size();
  int g_node, g_var, n;
  double g_rule;

  // Getting observations that are on the left and the ones that are in the right
  vector<int> new_left_index;
  vector<int> new_right_index;
  vector<int> curr_obs; // Observations that belong to that terminal node

  // Getting terminal nodes
  vector<node> t_nodes = curr_tree.getTerminals();
  n_terminals = t_nodes.size();

  while(bad_tree){

    // Sampling one of the nodes (this function will sample the node index from
    // the vector)
    g_node = sample_int(n_terminals);
    g_var = sample_int(p);

    // Number of observations in the node that will grow
    curr_obs = t_nodes[g_node].obs; // Current observations on that node
    n = t_nodes[g_node].obs.size();

    // Do not let split if is not possible to split
    if( floor(n/2) <= n_min_size){
      break;
    }

    g_rule = sample_double(x.col(g_var),n_min_size);

    // Iterating over the observations
    for(int i = 0; i<n; i++){
      if(x.coeff(curr_obs[i],g_var)<=g_rule){
        new_left_index.push_back(curr_obs[i]);
      } else {
        new_right_index.push_back(curr_obs[i]);
      }
    }

    if(new_right_index.size()>=n_min_size && new_left_index.size()>=n_min_size){

      // Getting the new nodes
      node new_node_left(n_nodes,new_left_index,-1,-1,t_nodes[g_node].depth+1,
                 g_var,g_rule,t_nodes[g_node].mu);

      node new_node_right(n_nodes+1,new_right_index,-1,-1,t_nodes[g_node].depth+1,
                         g_var,g_rule,t_nodes[g_node].mu);

      // Adding the new nodes
      curr_tree.list_node.push_back(new_node_left);
      curr_tree.list_node.push_back(new_node_right);

      // Transforming the growned node into non-terminal
      curr_tree.list_node[t_nodes[g_node].index].left = n_nodes;
      curr_tree.list_node[t_nodes[g_node].index].right = n_nodes+1;


      // Success to create a good tree;
      bad_tree = false;

    } else {

      // Cleaning the vectors
      new_left_index.clear();
      new_right_index.clear();
      curr_obs.clear();

      bad_tree_counter++;
      // Avoid infinite loop
      if(bad_tree_counter==2){
        bad_tree = false;
      }
    }


  }

  return curr_tree;


}

// Prune a tree
// [[Rcpp::depends(RcppEigen)]]
Tree prune(Tree curr_tree){

  // New list of nodes
  int n_nodes = curr_tree.list_node.size();
  vector<node> new_nodes;

  // Can't prune a root
  if(curr_tree.list_node.size()==1){
   return curr_tree;
  }

  // Getting the parents of terminal nodes
  vector<node> parent_left_right;
  for(int i=0;i<n_nodes;i++){
    if(!curr_tree.list_node[i].isTerminal()){
      if(curr_tree.list_node[curr_tree.list_node[i].left].isTerminal()){
        if(curr_tree.list_node[curr_tree.list_node[i].right].isTerminal()){
            parent_left_right.push_back(curr_tree.list_node[i]); // Adding the parent
          }
      }
    }
  }


  // Getting the number of internal and the node to be pruned
  int n_parent = parent_left_right.size();
  int p_node = sample_int(n_parent);


  // Returning to be a terminal node
  int left_index = curr_tree.list_node[parent_left_right[p_node].index].left;
  int right_index = curr_tree.list_node[parent_left_right[p_node].index].right;
  curr_tree.list_node[parent_left_right[p_node].index].left=-1;
  curr_tree.list_node[parent_left_right[p_node].index].right=-1;



  // Adding new trees

  for(int k = 0;k<curr_tree.list_node.size();k++){
    if((curr_tree.list_node[k].index != left_index) && (curr_tree.list_node[k].index != right_index)){


          // Doing for trees greater than left index
          if(curr_tree.list_node[k].index>left_index){ /// COULD ALSO USE IF curr_tree.list_node[k].depth >depth_prune_node
             curr_tree.list_node[k].index = curr_tree.list_node[k].index-2;
          }
          // Checking if the left index is greater than it should be
          if(curr_tree.list_node[k].left>left_index){
            curr_tree.list_node[k].left = curr_tree.list_node[k].left-2;
          }
          // Checking if the right index is greater than it should be
          if(curr_tree.list_node[k].right>left_index){
            curr_tree.list_node[k].right = curr_tree.list_node[k].right-2;
          }

      new_nodes.push_back(curr_tree.list_node[k]);
    }
  }


  // Updating the new nodes
  curr_tree.list_node = new_nodes;
  return curr_tree;
}

// Change a tree
// [[Rcpp::depends(RcppEigen)]]
Tree change(Tree curr_tree,
            Eigen::MatrixXd x,
            int n_min_size){

  // Declaring important values and variables
  int n_nodes = curr_tree.list_node.size();
  int p = x.cols(),n;
  int c_var,c_node;
  double c_rule;
  int bad_tree_counter=0;
  bool bad_tree = true;

  // Getting observations that are on the left and the ones that are in the right
  vector<int> new_left_index;
  vector<int> new_right_index;
  vector<int> curr_obs; // Observations that belong to that terminal node


  // Can't prune a root
  if(curr_tree.list_node.size()==1){
    return curr_tree;
  }

  // Getting the parents of terminal nodes (NOG)
  vector<node> parent_left_right;
  for(int i=0;i<n_nodes;i++){
    if(curr_tree.list_node[i].isTerminal()==0){
      if(curr_tree.list_node[curr_tree.list_node[i].left].isTerminal()==1 && curr_tree.list_node[curr_tree.list_node[i].right].isTerminal()==1){
        parent_left_right.push_back(curr_tree.list_node[i]); // Adding the parent
      }
    }
  }


  // Getting the number of parent of terminals and selecting the one to be changed
  int n_parent = parent_left_right.size();

  while(bad_tree){

    // Selecting the node and var
    c_node = sample_int(n_parent);
    c_var = sample_int(p);
    c_rule = sample_double(x.col(c_var),n_min_size);

    // Number of observations that will be changed
    curr_obs = parent_left_right[c_node].obs;
    n = parent_left_right[c_node].obs.size();

    // Iterating over the observations
    for(int i = 0; i<n; i++){
      if(x.coeff(curr_obs[i],c_var)<=c_rule){
        new_left_index.push_back(curr_obs[i]);
      } else {
        new_right_index.push_back(curr_obs[i]);
      }
    }

    // Verifying if is the correct min node size
    if(new_right_index.size()>=n_min_size && new_left_index.size()>=n_min_size){

      // Returning the nodes that will be changed
      int left_index = curr_tree.list_node[parent_left_right[c_node].index].left;
      int right_index = curr_tree.list_node[parent_left_right[c_node].index].right;


      // Changing the rules and observations that belong to each one of
      // that nodes
      curr_tree.list_node[left_index].var = c_var;
      curr_tree.list_node[right_index].var = c_var;
      curr_tree.list_node[left_index].var_split = c_rule;
      curr_tree.list_node[right_index].var_split = c_rule;
      curr_tree.list_node[left_index].obs = new_left_index;
      curr_tree.list_node[right_index].obs = new_right_index;

      // Success to create a good tree
      bad_tree = false;

    } else {

      // Cleaning the current vectors
      new_left_index.clear();
      new_right_index.clear();
      curr_obs.clear();

      bad_tree_counter++;

      // Avoid infinite while
      if(bad_tree_counter==2){
        bad_tree = false;
      }

    }

  }

  return curr_tree;

};

// Get the verbs
Tree swap(Tree curr_tree,
           Eigen::MatrixXd x,
           int n_min_size){

  // Testing if there are at least enough terminal nodes
  int n_int = curr_tree.n_internal();
  int branch_counter;
  // Declaring important values and variables
  int n_nodes = curr_tree.list_node.size();
  int n;

  int accept_new_tree = 0;
  Tree new_tree = curr_tree;


  // Getting observations that are on the left and the ones that are in the right
  vector<int> new_branch_left_index;
  vector<int> new_branch_right_index;

  vector<int> new_leaf_left_index;
  vector<int> new_leaf_right_index;

  vector<int> curr_branch_obs; // Observations that belong to that terminal node
  vector<int> curr_leaf_obs; // Observations that belong to that terminal node

  // Returning if there are not enough internal nodes
  if(n_int<3){
    return curr_tree;
  }

  // Getting the parents of terminal nodes (NOG)
  vector<node> swap_branches;

  for(int i=0;i<n_nodes;i++){
    if(curr_tree.list_node[i].isTerminal()==0){
      if((curr_tree.list_node[curr_tree.list_node[i].left].isTerminal()==1) || (curr_tree.list_node[curr_tree.list_node[i].right].isTerminal()==1)){
        branch_counter = 0;
        // Checking if the left child-branch is a parent of non-terminal nodes
        if((curr_tree.list_node[curr_tree.list_node[curr_tree.list_node[i].left].left].isTerminal()==1) && (curr_tree.list_node[curr_tree.list_node[curr_tree.list_node[i].left].right].isTerminal()==1)){
          swap_branches.push_back(curr_tree.list_node[i]);
          branch_counter++;
        }

        if ((curr_tree.list_node[curr_tree.list_node[curr_tree.list_node[i].right].left].isTerminal()==1) && (curr_tree.list_node[curr_tree.list_node[curr_tree.list_node[i].right].right].isTerminal()==1)){
          swap_branches.push_back(curr_tree.list_node[i]);
          branch_counter++;
        }


        // Getting all possible candidates to swap
        if(branch_counter==2){
          // cout << " BRANC COUNTER "<< branch_counter << endl;
          swap_branches.pop_back(); // Removing the branch node added twice
        }
      }
    }
  }


  // // Getting the possible number of parents to be swapped
  int n_swap = swap_branches.size();
  int swap_node = sample_int(n_swap);
  int child_swap_node_index = 0;
  double leaf_bottom_var;
  double leaf_bottom_rule;

  double branch_bottom_var;
  double branch_bottom_rule;

  // No internal with the proper conditions
  if(n_swap == 0){
    return curr_tree;
  }



  // Getting the values
  branch_bottom_var = curr_tree.list_node[swap_branches[swap_node].left].var;
  branch_bottom_rule = curr_tree.list_node[swap_branches[swap_node].left].var_split;


  // Get the child-branch
  if(curr_tree.list_node[swap_branches[swap_node].left].isTerminal()==1){
    child_swap_node_index = curr_tree.list_node[swap_branches[swap_node].right].index;
    leaf_bottom_var = curr_tree.list_node[curr_tree.list_node[swap_branches[swap_node].right].left].var;
    leaf_bottom_rule = curr_tree.list_node[curr_tree.list_node[swap_branches[swap_node].right].left].var_split;


  } else if(curr_tree.list_node[swap_branches[swap_node].right].isTerminal()==1){
    child_swap_node_index = curr_tree.list_node[swap_branches[swap_node].left].index;
    leaf_bottom_var = curr_tree.list_node[curr_tree.list_node[swap_branches[swap_node].left].left].var;
    leaf_bottom_rule = curr_tree.list_node[curr_tree.list_node[swap_branches[swap_node].left].left].var_split;

  }

  if(child_swap_node_index == 0 ){
    return curr_tree;
  }


  // Getting the values
  for(int l = 0;l<n_nodes;l++){

    if(curr_tree.list_node[l].index == swap_branches[swap_node].index){


      curr_branch_obs = curr_tree.list_node[l].obs;
      n = curr_tree.list_node[l].obs.size();

      for(int i=0;i<n;i++){
        if(x.coeff(curr_branch_obs[i],leaf_bottom_var)<=leaf_bottom_rule){
          new_branch_left_index.push_back(curr_branch_obs[i]);
        } else {
          new_branch_right_index.push_back(curr_branch_obs[i]);
        }
      }

      // Certifying that I will have only nodes with enough observations
      if(new_branch_right_index.size()>=n_min_size && new_branch_left_index.size()>=n_min_size){
        new_tree.list_node[new_tree.list_node[l].left].var = leaf_bottom_var;
        new_tree.list_node[new_tree.list_node[l].left].var_split = leaf_bottom_rule;
        new_tree.list_node[new_tree.list_node[l].right].var = leaf_bottom_var;
        new_tree.list_node[new_tree.list_node[l].right].var_split = leaf_bottom_rule;

        // Replacing the left and right index
        new_tree.list_node[new_tree.list_node[l].left].obs = new_branch_left_index;
        new_tree.list_node[new_tree.list_node[l].right].obs = new_branch_right_index;


        accept_new_tree++;
        // Cleaning the auxiliary vectors

        // new_left_index.clear();
        // new_right_index.clear();
        // curr_obs.clear();

      } else {
        // cout << "NO VALID TREE FIRST" <<  endl;
        return curr_tree; // Finish the VERB
      }

      // Skipping for the bottom leaves nodes
    }

  }

  for(int l = 0; l<n_nodes;l++){

    if(curr_tree.list_node[l].index==child_swap_node_index){

      curr_leaf_obs = new_tree.list_node[l].obs;
      n = new_tree.list_node[l].obs.size();

      for(int i=0;i<n;i++){
        if(x.coeff(curr_leaf_obs[i],branch_bottom_var)<=branch_bottom_rule){
          new_leaf_left_index.push_back(curr_leaf_obs[i]);
        } else {
          new_leaf_right_index.push_back(curr_leaf_obs[i]);
        }
      }

      // Certifying that I will have only nodes with enough observations
      if(new_leaf_right_index.size()>=n_min_size && new_leaf_left_index.size()>=n_min_size){
        new_tree.list_node[new_tree.list_node[l].left].var = branch_bottom_var;
        new_tree.list_node[new_tree.list_node[l].left].var_split = branch_bottom_rule;
        new_tree.list_node[new_tree.list_node[l].right].var = branch_bottom_var;
        new_tree.list_node[new_tree.list_node[l].right].var_split = branch_bottom_rule;

        // Replacing the left and right index
        new_tree.list_node[new_tree.list_node[l].left].obs = new_leaf_left_index;
        new_tree.list_node[new_tree.list_node[l].right].obs = new_leaf_right_index;

        accept_new_tree++;
      } else {
        // cout << "NO VALID TREE SECOND" <<  endl;

        return curr_tree; // Finish the VERB
      }

    }
  }

  if(accept_new_tree==2){
    // cout <<  "   " << endl;
    // cout << " THE NODE WAS SWAPPPEED" << endl;
    // cout <<  "   " << endl;
    return new_tree;
    // return curr_tree;
  } else {
    return curr_tree;
  }
}


// [[Rcpp::depends(RcppEigen)]]
double node_loglikelihood(Eigen::VectorXd residuals,
                          node current_node,
                          double tau,
                          double tau_mu) {

  // Decarling the quantities
  int n_size = current_node.obs.size();
  double sum_sq_r = 0 , sum_r = 0;

  for(int i = 0;i<n_size;i++){
    sum_sq_r+=residuals.coeff(i)*residuals.coeff(i);
    sum_r+=residuals.coeff(i);

  }

  return -0.5*tau*sum_sq_r-0.5*tau*tau*sum_r*sum_r/(tau*n_size+tau_mu)-0.5*log(tau*n_size+tau_mu);
}

// [[Rcpp::depends(RcppEigen)]]
double tree_loglikelihood(Eigen::VectorXd residuals,
                          Tree curr_tree,
                          double tau,
                          double tau_mu
) {

  // Declaring important quantities
  vector<node> terminal_nodes = curr_tree.getTerminals();
  int number_nodes = terminal_nodes.size();
  double loglike_sum = 0;


  for(int i = 0; i<number_nodes; i++) {
      loglike_sum+=node_loglikelihood(residuals,
                                      terminal_nodes[i],tau,tau_mu);
  }

  return loglike_sum;
}

// Updating the \mu
// [[Rcpp::depends(RcppEigen)]]
Tree update_mu(Eigen::VectorXd residuals,
               Tree curr_tree,
               double tau,
               double tau_mu){

  // Declaring important quantities
  vector<node> terminal_nodes = curr_tree.getTerminals();
  int number_nodes = terminal_nodes.size();
  double sum_residuals;

  // Iterating over terminal nodes
  for(int i = 0;i<number_nodes;i++){
    sum_residuals = 0;
    for(int j = 0;j<terminal_nodes[i].obs.size();j++){
      sum_residuals+=residuals.coeff(terminal_nodes[i].obs[j]);
    }
    curr_tree.list_node[terminal_nodes[i].index].mu = R::rnorm( ((tau)/(terminal_nodes[i].obs.size()*tau+tau_mu))*sum_residuals,sqrt(1/(terminal_nodes[i].obs.size()*tau+tau_mu)) );
  }

  return curr_tree;
}

// Calculate the density of a half cauchy location 0 and sigma
//[[Rcpp::export]]
double dhcauchy(double x, double location, double sigma){

  if( x>location) {
    return (1/(M_PI_2*sigma))*(1/(1+((x-location)*(x-location))/(sigma*sigma)));
  } else {
    return 0.0;
  }
}

double update_tau_old(Eigen::VectorXd y,
                  Eigen::VectorXd y_hat,
                  double a_tau,
                  double d_tau){

  // Function used in the development of the package where I checked
  // contain_nan(y_hat);

  int n = y.size();
  double sum_sq_res=0;
  for(int i=0;i<n; i++){
      sum_sq_res += (y.coeff(i)-y_hat.coeff(i))*(y.coeff(i)-y_hat.coeff(i));
  }
  return R::rgamma((0.5*n+a_tau),1/(0.5*sum_sq_res+d_tau));
}

double update_tau(Eigen::VectorXd y,
                  Eigen::VectorXd y_hat,
                  double naive_sigma,
                  double curr_tau){

      int n=y.size();
      double curr_sigma, proposal_tau,proposal_sigma, acceptance;

      curr_sigma = 1/(sqrt(curr_tau));

      // Getting the sum squared of residuals
      double sum_sq_res=0;
        for(int i=0;i<n; i++){
            sum_sq_res += ((y.coeff(i)-y_hat.coeff(i)))*((y.coeff(i)-y_hat.coeff(i)));
      }

      // Getting a proposal sigma
      proposal_tau = R::rgamma(0.5*n+1,1/(0.5*sum_sq_res));
      proposal_sigma  = 1/(sqrt(proposal_tau));

      acceptance = exp(log(dhcauchy(proposal_sigma,0,naive_sigma))+3*log(proposal_sigma)-log(dhcauchy(curr_sigma,0,naive_sigma))-3*log(curr_sigma));

      if(R::runif(0,1)<=acceptance){
        return 1/(proposal_sigma*proposal_sigma);
      } else{
        return curr_tau;
      }
}





// [[Rcpp::depends(RcppEigen)]]
Eigen::VectorXd get_prediction_tree(Tree curr_tree,
                             Eigen::MatrixXd x,
                             bool is_train){

  // Creating the prediction vector
  int n = x.rows();
  Eigen::VectorXd prediction(n);
  prediction.setZero(n);
  vector<int> curr_obs;
  vector<int> new_left_index, new_right_index;

  // Getting terminal nodes
  vector<node> terminal_nodes = curr_tree.getTerminals();
  int n_terminals = terminal_nodes.size();

  // If x is the training data
  if(is_train){
      // Iterating to get the predictions for a training test
      for(int i = 0;i<n_terminals;i++){
        for(int j = 0; j<terminal_nodes[i].obs.size();j++){
          prediction[terminal_nodes[i].obs[j]] = terminal_nodes[i].mu;
        }
      }
  } else {

    curr_tree.list_node[0].obs = seq_along(n);

    for(int i = 0; i<curr_tree.list_node.size();i++){
      curr_obs = curr_tree.list_node[i].obs;

      // If the node is terminal update the new predictions with the mu values
      //from the terminal node;

      if(curr_tree.list_node[i].isTerminal()==1){
        for(int j=0;j<curr_obs.size();j++){
          prediction(curr_obs[j]) = curr_tree.list_node[i].mu;
        }
      } else {

        for(int j=0; j<curr_obs.size();j++){

          if(x.coeff(j,curr_tree.list_node[curr_tree.list_node[i].left].var)<=curr_tree.list_node[curr_tree.list_node[i].left].var_split){
              new_left_index.push_back(curr_obs[j]);
          } else {
              new_right_index.push_back(curr_obs[j]);
          }
        }

        // Storing the new ones
        curr_tree.list_node[curr_tree.list_node[i].left].obs=new_left_index;
        curr_tree.list_node[curr_tree.list_node[i].right].obs=new_right_index;

        // Clearning the aux vectors afterwards
        new_left_index.clear();
        new_right_index.clear();

      }

    }

  }

  return prediction;
}

// [[Rcpp::depends(RcppEigen)]]
double tree_log_prior(Tree curr_tree,
                  double alpha,
                  double beta) {

  double log_tree_p = 0;

  for(int i=0;i<curr_tree.list_node.size();i++){
    if(curr_tree.list_node[i].isTerminal()==1){
      log_tree_p += log(1-alpha/pow((1+curr_tree.list_node[i].depth),beta));
    } else {
      log_tree_p += log(alpha)-beta*log(1+curr_tree.list_node[i].depth);
    }
  }

  return log_tree_p;
}

// Get the log transition probability
// [[Rcpp::depends(RcppEigen)]]
double log_transition_prob(Tree curr_tree,
                           Tree new_tree,
                           double verb){

  // Getting the probability
  double log_prob = 0;
  // In case of Grow: (Prob from Grew to Current)/(Current to Grow)
  if(verb < 0.3){
    log_prob = log(0.3/new_tree.n_nog())-log(0.3/curr_tree.n_terminal());
  } else if (verb < 0.6) { // In case of Prune: (Prob from the Pruned to the current)/(Prob to the current to the prune)
    log_prob = log(0.3/new_tree.n_terminal())-log(0.3/curr_tree.n_nog());
  }; // In case of change log_prob = 0; it's already the actual value

  return log_prob;

}

//[[Rcpp::export]]
List bart(Eigen::MatrixXd x,
          Eigen::VectorXd y,
          int n_tree,
          int n_mcmc,
          int n_burn,
          int n_min_size,
          double tau, double mu,
          double tau_mu, double naive_sigma,
          double a_tau, double d_tau,
          double alpha, double beta){

  // Declaring common variables
  double verb;
  double acceptance;
  double log_transition_prob_obj;
  double past_tau;
  int post_counter = 0;

  // Getting the number of observations
  int n = x.rows();

  // Creating the variables
  int n_post = (n_mcmc-n_burn);
  Eigen::MatrixXd y_hat_post(n_post,n);
  NumericVector tau_post;


  // Creating a vector of multiple trees
  vector<Tree> current_trees(n_tree,n);
  Tree new_tree(n);


  // Creating a matrix of zeros of y_hat
  y_hat_post.setZero(n_post,n);

  // Creating the partial residuals and partial predictions
  Eigen::VectorXd partial_pred,partial_residuals;

  // Initializing zero values
  partial_pred.setZero(n);
  partial_residuals.setZero(n);


  // Creating loglikelhood objects
  double log_like_old,log_like_new;

  // Iterating over all MCMC samples
  for(int i=0;i<n_mcmc;i++){

      // Iterating over the trees
      for(int t = 0; t<n_tree;t++){


        // Updating the residuals
        partial_residuals = y - (partial_pred - get_prediction_tree(current_trees[t],x,true));

        // Sampling a verb. In this case I will sample a double between 0 and 1
        // Grow: 0-0.3
        // Prune 0.3 - 0.6
        // Change: 0.6 - 1.0
        // Swap: Not in this current implementation
        verb = R::runif(0,1);

        // Forcing the first trees to grow
        if(current_trees[t].list_node.size()==1 ){
          verb = 0.1;
        }

        // Proposing a new tree
        if(verb<=0.25){
              new_tree = grow(current_trees[t], x, n_min_size);
        } else if ( verb <= 0.5){
              if(current_trees[t].list_node.size()>1){
                  new_tree = prune(current_trees[t]);
                }
        } else if (verb <= 0.9){
               new_tree = change(current_trees[t],x,n_min_size);
        } else {
               new_tree = swap(current_trees[t],x,n_min_size);
        }

        // No new tree is proposed (Jump log calculations)
        if( (verb <=0.6) && (current_trees[t].list_node.size()==new_tree.list_node.size())){
          log_transition_prob_obj = 0;
        } else {
          log_transition_prob_obj = log_transition_prob(current_trees[t],new_tree,verb);

        }

        // Calculating all log_likelihoods
        log_like_old = tree_loglikelihood(partial_residuals,current_trees[t],tau,tau_mu) + tree_log_prior(current_trees[t],alpha,beta);
        // log_like_old =  tree_log_prior(current_trees[t],alpha,beta);

        log_like_new = tree_loglikelihood(partial_residuals,new_tree,tau,tau_mu) + tree_log_prior(new_tree,alpha,beta);
        // log_like_new =  tree_log_prior(new_tree,alpha,beta);

        // acceptance = log_like_new-log_like_old + log_transition_prob_obj;
        acceptance = 0.1;

        // Testing if will acceptance or not
        if( (acceptance > 0) || (acceptance > -R::rexp(1))){
          current_trees[t] = new_tree;
        }

        // Generating new \mu values for the accepted (or not tree)
        current_trees[t] = update_mu(partial_residuals,current_trees[t],tau,tau_mu);

        partial_pred = y + get_prediction_tree(current_trees[t],x,true) - partial_residuals;

      }

      // Updating tau
      past_tau = tau;
      tau = update_tau(y,partial_pred,naive_sigma,past_tau);
      // tau = update_tau_old(y,partial_pred,a_tau,d_tau);

      // Updating the posterior matrix
      if(i>=n_burn){
        y_hat_post.row(post_counter) += partial_pred;
        //Updating tau
        tau_post.push_back(tau);
        post_counter++;
      }

  }

  return Rcpp::List::create(_["y_hat_post"] = y_hat_post,_["tau_post"] = tau_post);

}


//[[Rcpp::export]]
int initialize_test(Eigen::MatrixXd x,Eigen::MatrixXd x_new,
                    Eigen::VectorXd residuals,
                    double tau, double tau_mu){

  int n_obs(x.rows());
  Tree tree1(n_obs);
  Tree new_tree = grow(tree1,x,2);
  Tree new_tree_two = grow(new_tree,x,2);
  Tree new_tree_three = grow(new_tree_two,x,2);
  Tree new_tree_four = grow(new_tree_three,x,2);
  Tree new_tree_five = grow(new_tree_four,x,2);
  Tree new_tree_six = grow(new_tree_five,x,2);
  new_tree_four.DisplayNodes();
  // // Tree new_three_tree = grow(new_tree_two,x,2);

  // cout << "======= SWAPP LINE =======" << endl;
  Tree swap_tree(n_obs);
  // get_prediction_tree(new_tree_five,x,true);
  // for(int k = 0;k<1000;k++){
    // cout << k << endl;
    swap_tree = swap(new_tree_four,x,2);
    swap_tree.DisplayNodes();
  // }
  // Tree update_mu_tree = update_mu(residuals,change_tree_two,tau,tau_mu);
  // Eigen::VectorXd pred_vec  = get_prediction_tree(update_mu_tree,x_new,false);
  // // get_prediction_tree(change_tree_two,x,true);
  // cout << "Pred x: ";
  // for(int i = 0; i<pred_vec.size();i++){
  //   cout  << pred_vec.coeff(i)<< " " ;
  // }

  // cout << "Tree prior" << tree_log_prior(new_tree_two,0.95,2) << endl;
  // vector<int> p_index = new_tree_two.getParentTerminals();

  // Tree prune_tree = prune(new_three_tree);
  // prune_tree.DisplayNodes();

  return 0;

}


