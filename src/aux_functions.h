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
double sample_double(VectorXd vec){

  // Getting the range of the vector
  return R::runif(vec.minCoeff(),vec.maxCoeff());
}



