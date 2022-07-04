#' @useDynLib cppbart
#' @importFrom Rcpp sourceCpp

# Calling the function in the R version
#' @export
r_bart <- function(x,
                   y,
                   n_tree,
                   n_mcmc,
                   n_burn,
                   n_min_size,
                   tau, mu,
                   a_tau,d_tau, tau_mu,
                   alpha, beta){

  # Getting the main function
  bart_obj <- bart(x = x,
                   y = y,
                   n_tree = n_tree,
                   n_mcmc = n_mcmc,
                   n_burn = n_burn,
                   n_min_size = n_min_size,
                   tau = tau, mu = mu,
                   a_tau = a_tau,d_tau = d_tau, tau_mu = tau_mu,
                   alpha = alpha, beta = beta)
  return(bart_obj)
}
