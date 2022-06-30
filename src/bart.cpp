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

    if(new_right_index.size()>=n_min_size && new_left_index.size()>=n_min_size){

      // Getting the new nodes
      node new_node_left(n_nodes,new_left_index,-1,-1,t_nodes[g_node].depth+1,
                 g_var,g_rule,t_nodes[g_node].mu);

      node new_node_right(n_nodes+1,new_right_index,-1,-1,t_nodes[g_node].depth+1,
                         g_var,g_rule,t_nodes[g_node].mu);

      // Adding the new nodes
      new_tree.list_node.push_back(new_node_left);
      new_tree.list_node.push_back(new_node_right);

      // Transforming the growned node into non-terminal
      new_tree.list_node[t_nodes[g_node].index].left = n_nodes;
      new_tree.list_node[t_nodes[g_node].index].right = n_nodes+1;


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
  cout << "#========#" << endl;

  return new_tree;


}

// Prune a tree
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
        cout<< " The parent node is given by "<<i<<endl;
      }
    }
  }


  // Transforming back the parent and left nodes
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




//[[Rcpp::export]]
int initialize_test(MatrixXd x){

  vector<int> vec;
  for(int i = 0; i<3;i++){
    vec[i] = i;
  }

  Tree tree1(10);
  node node1(1,vec,1,1,1,1,0.1,0.1);
  //node1.DisplayNode();
  Tree new_tree = grow(tree1,x,2);
  // new_tree.DisplayNodes();
  Tree new_tree_two = grow(new_tree,x,2);
  // new_tree_two.DisplayNodes();
  Tree new_three_tree = grow(new_tree_two,x,2);
  new_three_tree.DisplayNodes();

  // vector<int> p_index = new_tree_two.getParentTerminals();

  Tree prune_tree = prune(new_three_tree);
  prune_tree.DisplayNodes();

  return 0;

}
