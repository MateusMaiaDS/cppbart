#include<iostream>
#include <vector>
#include <math.h>
#include<RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>

using Eigen::VectorXd;
using namespace std;

// Function to sample an integer from a sequence from
// 0:n
int sample_int(int n){
  return floor(R::runif(0,n));
}

// Sample an uniform value from a double vector
double sample_double(VectorXd vec, int n_min_size){

  // Getting the range of the vector
  std::sort(vec.data(),vec.data()+vec.size() );
  return R::runif(vec[(n_min_size-1)],vec[(vec.size()-(n_min_size+1))]);

}



