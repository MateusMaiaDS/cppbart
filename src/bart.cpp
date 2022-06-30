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
using Eigen::LLT;
using Eigen::Lower;
using Eigen::Upper;

using namespace Rcpp;



//[[Rcpp::export]]
Tree grow(Tree curr_tree,
          Eigen::MatrixXd x,
          int n_min_size){

  // Declaring variables
  bool bad_trees=true;
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

  // Getting terminal nodes
  vector<node> t_nodes = curr_tree.getTerminals();
  n_terminals = t_nodes.size();

  // Defining new tree
  Tree new_tree = curr_tree;


  while(bad_trees){

    // Sampling one of the nodes (this function will sample the node index from
    // the vector)
    g_node = sample_int(n_terminals);
    g_var = sample_int(p);

    // Number of observations in the node that will grow
    vector<int> curr_obs = t_nodes[g_node].obs;
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

    if(new_right_index.size()>=n_min_size & new_left_index.size()>=n_min_size){

      // Getting the new nodes
      node new_node_left(n_nodes+1,new_left_index,-1,-1,t_nodes[g_node].depth+1,
                 g_var,g_rule,t_nodes[g_node].mu);

      node new_node_right(n_nodes+2,new_right_index,-1,-1,t_nodes[g_node].depth+1,
                         g_var,g_rule,t_nodes[g_node].mu);

      // Adding the new nodes
      new_tree.list_node.push_back(new_node_left);
      new_tree.list_node.push_back(new_node_right);

      // Transforming the growned node into non-terminal
      new_tree.list_node[g_node].left = n_nodes+1;
      new_tree.list_node[g_node].right = n_nodes+2;


      // Success to create a good tree;
      bad_trees = false;

    } else {

      // Cleaning the vectors
      new_left_index.clear();
      new_right_index.clear();

      bad_tree_counter++;
      // Avoid infinite loop
      if(bad_tree_counter == 2){
        bad_trees = false;
      }
    }


  }

  return new_tree;


}


//[[Rcpp::export]]
int initialize_test(MatrixXd x){

  vector<int> vec;
  for(int i = 0; i<3;i++){
    vec[i] = i;
  }

  Tree tree1(10);
  node node1(1,vec,1,1,1,1,0.1,0.1);
  //node1.DisplayNode();
  Tree new_tree = grow(tree1,x,3);
  //new_tree.DisplayNodes();
  Tree new_tree_two = grow(new_tree,x,18);
  new_tree_two.DisplayNodes();

  return 0;

}
