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
  cout << "#========#" << endl;

  return curr_tree;


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

  // Getting the parents of terminal nodes
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
  new_tree_two.DisplayNodes();
  Tree change_tree_two = change(new_tree_two,x,2);
  change_tree_two.DisplayNodes();
  // Tree new_three_tree = grow(new_tree_two,x,2);
  // new_three_tree.DisplayNodes();

  // vector<int> p_index = new_tree_two.getParentTerminals();

  // Tree prune_tree = prune(new_three_tree);
  // prune_tree.DisplayNodes();

  return 0;

}
