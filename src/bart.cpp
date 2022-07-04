#include <Rcpp.h>
#include<RcppEigen.h>
#include <vector>
#include <math.h>
#include "aux_functions.h"
#include "tree.h"

// [[Rcpp::depends(RcppEigen)]]
// The line above (depends) it will make all the dependcies be included on the file
using namespace Rcpp;
using namespace RcppEigen;
using namespace std;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Upper;


// Nan error checker
void contain_nan(VectorXd vec){

  for(int i = 0; i<vec.size();i++){
    if(isnan(vec.coeff(i))==1){
      cout << "This vector contain nan" << endl;
    }
  }

}

void print_vector(VectorXd vector) {

  for(int k=0;k<vector.size();k++){
    cout << vector[k] << " ";
  }
  cout << endl;
}

// [[Rcpp::depends(RcppEigen)]]
Tree grow(Tree curr_tree,
          Eigen::MatrixXd x,
          int n_min_size){

  // Declaring variables
  bool bad_tree=true;
  int bad_tree_counter = 0;

  // Getting information from x
  int p(x.cols());
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

    g_rule = sample_double(x.col(g_var));

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
      if(bad_tree_counter == 2){
        bad_tree = false;
      }
    }


  }

  return curr_tree;


}

// Prune a tree
// [[Rcpp::depends(RcppEigen)]]
Tree prune(Tree curr_tree){

  int n_nodes = curr_tree.list_node.size();

  // Can't prune a root
  if(curr_tree.list_node.size()==1){
    return curr_tree;
  }

  // Getting the parents of terminal nodes
  vector<node> parent_left_right;
  for(int i=0;i<n_nodes;i++){
    if(curr_tree.list_node[i].isTerminal()==0){
      if(curr_tree.list_node[curr_tree.list_node[i].left].isTerminal()==1 && curr_tree.list_node[curr_tree.list_node[i].right].isTerminal()==1){
            parent_left_right.push_back(curr_tree.list_node[i]); // Adding the parent
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

  // Removing the two tree elements ( The plust 1 inn the right because erase([,)]))
  curr_tree.list_node.erase(curr_tree.list_node.begin()+left_index,curr_tree.list_node.begin()+right_index+1);

  return curr_tree;
}

// Change a tree
// [[Rcpp::depends(RcppEigen)]]
Tree change(Tree curr_tree,
            MatrixXd x,
            int n_min_size){

  // Declaring important values and variables
  int n_nodes = curr_tree.list_node.size();
  int p(x.cols()),n;
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
    c_rule = sample_double(x.col(c_var));

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

// // Get the verbs
// Tree swap(Tree curr_tree,
//            MatrixXd x,
//            int n_min_size){
//
//   // Testing if there are at least enough terminal nodes
//   int n_int = curr_tree.n_internal();
//
//   // Returning if there are not enough internal nodes
//   if(n_int<2){
//     return curr_tree;
//   }
//
//   return curr_tree;
// }

// [[Rcpp::depends(RcppEigen)]]
double node_loglikelihood(VectorXd residuals,
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
double tree_loglikelihood(VectorXd residuals,
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
Tree update_mu(VectorXd residuals,
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
    curr_tree.list_node[terminal_nodes[i].index].mu = R::rnorm( ((tau)/(terminal_nodes[i].obs.size()*tau+tau_mu))*sum_residuals,sqrt(1/((terminal_nodes[i].obs.size()*tau+tau_mu))) );
  }

  return curr_tree;
}

// [[Rcpp::depends(RcppEigen)]]
double update_tau(VectorXd y,
                  VectorXd y_hat,
                  double a_tau,
                  double d_tau){

  contain_nan(y_hat);

  int n = y.size();
  double sum_sq_res=0;
  sum_sq_res = ((y-y_hat)*(y-y_hat)).sum();
  return R::rgamma((0.5*n+a_tau),1/(0.5*sum_sq_res+d_tau));
}

// [[Rcpp::depends(RcppEigen)]]
VectorXd get_prediction_tree(Tree curr_tree,
                             MatrixXd x,
                             bool is_train){

  // Creating the prediction vector
  int n(x.rows());
  VectorXd prediction(n);
  prediction.setZero(n);
  vector<int> curr_obs;
  vector<int> new_left_index, new_right_index;

  // Getting terminal nodes
  vector<node> terminal_nodes = curr_tree.getTerminals();
  int n_terminals = terminal_nodes.size();

  // If x is the training data
  if(is_train==1){
      // Iterating to get the predictions for a training test
      for(int i = 0;i<n_terminals;i++){
        for(int j = 0; j<terminal_nodes[i].obs.size();j++){
          prediction(terminal_nodes[i].obs[j]) = terminal_nodes[i].mu;
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

          if(x(j,curr_tree.list_node[curr_tree.list_node[i].left].var)<=curr_tree.list_node[curr_tree.list_node[i].left].var_split){
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

  contain_nan(prediction);
  return prediction;
}

// [[Rcpp::depends(RcppEigen)]]
double tree_log_prior(Tree curr_tree,
                  double alpha,
                  double beta) {

  double log_tree_p = 0;

  for(int i=0;i<curr_tree.list_node.size();i++){
    if(curr_tree.list_node[i].isTerminal()==1){
      // cout << log(1-alpha/pow((1+curr_tree.list_node[i].depth),beta)) << endl;
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
    log_prob = log(0.3/new_tree.n_nog())-log(0.3/new_tree.n_terminal());
  } else if (verb < 0.6) { // In case of Prune: (Prob from the Pruned to the current)/(Prob to the current to the prune)
    log_prob = log(0.3/new_tree.n_terminal())-log(0.3/new_tree.n_nog());
  }; // In case of change log_prob = 0; it's already the actual value

  return log_prob;

}






//[[Rcpp::export]]
MatrixXd bart(MatrixXd x,
          VectorXd y,
          int n_tree,
          int n_mcmc,
          int n_burn,
          int n_min_size,
          double tau, double mu,
          double a_tau, double d_tau, double tau_mu,
          double alpha, double beta){

  // Declaring common variables
  double verb;
  double acceptance;
  int post_counter = 0, a_counter = 0;

  // Getting the number of observations
  int n(x.rows());

  // Creating the variables
  int n_post = (n_mcmc-n_burn);
  MatrixXd y_hat_post;

  // Creating a vector of multiple trees
  vector<Tree> current_trees;
  Tree new_tree(n);

  // Creating a matrix of zeros of y_hat
  y_hat_post.setZero(n_post,n);

  // Creating the partial residuals and partial predictions
  VectorXd partial_pred,partial_residuals;

  // Initializing zero values
  partial_pred.setZero(n);
  partial_residuals.setZero(n);

  // Initializing the "n_trees"
  for(int i=0; i<n_tree;i++){
    current_trees.push_back(Tree(n));

  }

  // Creating loglikelhood objects
  double log_like_old,log_like_new;

  // Iterating over all MCMC samples
  for(int i=0;i<n_mcmc;i++){



      // Iterating over the trees
      for(int t = 0; t<n_tree;t++){

        // cout << "Y:" << endl;
        // print_vector(get_prediction_tree(current_trees[t],x,true));

        // Updating the residuals
        partial_residuals = y - (partial_pred - get_prediction_tree(current_trees[t],x,true));

        // Sampling a verb. In this case I will sample a double between 0 and 1
        // Grow: 0-0.3
        // Prune 0.3 - 0.6
        // Change: 0.6 - 1.0
        // Swap: Not in this current implementation
        verb = R::runif(0,1);

        // Proposing a new tree
        if(verb<=0.3){
          new_tree = grow(current_trees[t], x, n_min_size);
        } else if ( verb <= 0.6){
          new_tree = prune(current_trees[t]);
        } else {
          new_tree = change(current_trees[t],x,n_min_size);
        }

        // Calculating all log_likelihoods
        log_like_old = tree_loglikelihood(partial_residuals,current_trees[t],tau,tau_mu) + tree_log_prior(current_trees[t],alpha,beta);

        log_like_new = tree_loglikelihood(partial_residuals,new_tree,tau,tau_mu) + tree_log_prior(new_tree,alpha,beta);

        acceptance = log_like_new-log_like_old + log_transition_prob(current_trees[t],new_tree,verb);

        // Testing if will acceptance or not
        if( (acceptance > 0) || (acceptance > -R::rexp(1))){
          // cout << a_counter << endl;
          current_trees[t] = new_tree;
          // a_counter++;
        }

        // Generating new \mu values for the accepted (or not tree)
        current_trees[t] = update_mu(partial_residuals,current_trees[t],tau,tau_mu);

        partial_pred = y + get_prediction_tree(current_trees[t],x,true) - partial_residuals;

      }

      // Updating the posterior matrix
      if(n_mcmc>n_burn){
        y_hat_post.row(post_counter) = partial_pred;
        post_counter++;
      }

      //Updating tau
      // cout << update_tau(y,partial_pred,a_tau,d_tau) << endl;
      // tau = update_tau(y,partial_pred,a_tau,d_tau);
  }

  // Cleaning some memory;
  current_trees.clear();

  return y_hat_post;
}

//[[Rcpp::export]]
int initialize_test(MatrixXd x,MatrixXd x_new,
                    VectorXd residuals,
                    double tau, double tau_mu){

  // int n_obs(x.rows());
  // vector<int> vec;
  // for(int i = 0; i<3;i++){
  //   vec[i] = i;
  // }
  //
  // Tree tree1(n_obs);
  // node node1(1,vec,1,1,1,1,0.1,0.1);
  // //node1.DisplayNode();
  // Tree new_tree = grow(tree1,x,2);
  // // new_tree.DisplayNodes();
  // Tree new_tree_two = grow(new_tree,x,2);
  // // new_tree_two.DisplayNodes();
  // Tree change_tree_two = change(new_tree_two,x,2);
  // // change_tree_two.DisplayNodes();
  // // Tree new_three_tree = grow(new_tree_two,x,2);
  // // new_three_tree.DisplayNodes();
  // Tree update_mu_tree = update_mu(residuals,change_tree_two,tau,tau_mu);
  // VectorXd pred_vec  = get_prediction_tree(update_mu_tree,x_new,false);
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


