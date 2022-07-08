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
                   alpha, beta,
                   scale_boolean = TRUE,
                   K_bart = 2){


  # Saving a_min and b_max
  a_min <- NULL
  b_max <- NULL

  a_min <- min(y)
  b_max <- max(y)

  # Scale values
  if(scale_boolean) {
    # Normalizing y
    y_scale <- normalize_bart(y = y)

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu <- (4 * n_tree * (K_bart^2))

    nsigma <- naive_sigma(x = x,y = y_scale)

  } else {

    # Not scaling the y
    y_scale <- y

    # Calculating \tau_{\mu} based on the scale of y
    # Need to change this value in case of non-scaling
    tau_mu <- (4 * n_tree * (K_bart^2))

    nsigma <- naive_sigma(x = x,y = y_scale)


  }

  # Getting the main function
  bart_obj <- bart(x = x,
                   y = y_scale,
                   n_tree = n_tree,
                   n_mcmc = n_mcmc,
                   n_burn = n_burn,
                   n_min_size = n_min_size,
                   tau = tau, mu = mu,
                   tau_mu = tau_mu,
                   naive_sigma = nsigma, # naive sigma value
                   alpha = alpha, beta = beta)

  # Returning to the normal scale
  if(scale_boolean){
    bart_obj$y_hat_post <- unnormalize_bart(bart_obj$y_hat_post,
                                            a = a_min, b = b_max)
    bart_obj$tau_post <- bart_obj$tau_post/((b_max-a_min)^2)
  }



  return(bart_obj)
}


zero_tau_prob_squared <- function(x, naive_tau_value, prob, shape) {

  # Find the zero to the function P(tau < tau_ols) = 0.1, for a defined
  return((stats::pgamma(naive_tau_value,
                        shape = shape,
                        scale = x) - (1 - prob))^2)
}

# Naive sigma_estimation
naive_sigma <- function(x,y){

  # Getting the valus from n and p
  n <- length(y)

  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))

  # Adjusting the df
  df <- data.frame(x,y)
  colnames(df)<- c(colnames(x),"y")

  # Naive lm_mod
  lm_mod <- stats::lm(formula = y ~ ., data =  df)

  # Getting sigma
  sigma <- stats::sigma(lm_mod)


  # sigma <- sqrt(sum((lm_mod$residuals)^2)/(n - p))
  # sigma <- stats::sd(y)
  sigma <- stats::sigma(lm_mod)

  return(sigma)
}


# Normalize BART function (Same way as theOdds code)
normalize_bart <- function(y) {

  # Defining the a and b
  a <- min(y)
  b <- max(y)

  # This will normalize y between -0.5 and 0.5
  y  <- (y - a)/(b - a) - 0.5
  return(y)
}

# Now a function to return everything back to the normal scale

unnormalize_bart <- function(z, a, b) {
  # Just getting back to the regular BART
  y <- (b - a) * (z + 0.5) + a
  return(y)
}
