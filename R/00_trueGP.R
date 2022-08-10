#' Simulates a Gaussian process with a given covariance function.
#'
#' This function simulates a Gaussian process over a fine-grid,
#'  provided with a given covariance function.
#'
#' @param from numeric value for the left limit of the interval of interest
#' @param to numeric value for the right limit of the interval of interest
#' @param start numerical of the starting value of the GP, default to zero.
#' @param K A covariance function of two arguments that is PD. Default being Brownian motion.
#' @param m An integer value of number of points to simulate.
#' @return A dataframe with time index and value of the sample path.
#' @export
#' @examples
#' gaussprocess(from = 0, to = 1, K = function(s, t) {min(s, t)}, start = 0, m = 10)
gaussprocess <- function(from = 0, to = 1, K = function(s, t) {min(s, t)},
                         start = 0, m = 1000) {
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })

  path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
  path <- path - path[1] + start

  return(data.frame("t" = t, "xt" = path))
}



#' Generate the true covariance function of the sGP, given its periodicity parameter
#' and its SD parameter.
#'
#' This function takes the periodicity parameter and SD parameter to generate a covariance function
#' that correspond to the specified sGP.
#'
#' @param sigma The value of the SD parameter.
#' @param alpha The value of the periodicity parameter.
#' @return A PD covariance function that takes two arguments.
#' @export
#' @examples
#' fun <- generate_K_true(1,1)
#' fun(1,2)
generate_K_true <- function(sigma,alpha){
  K_true <- function(s,t){
    ((sigma/alpha)^2)*((min(s,t)/2)*cos(alpha*(abs(s-t))) - (cos(alpha*max(c(s,t)))*sin(alpha*min(c(s,t))))/(2*alpha))
  }
  K_true
}



#' Compute covariance matrix given covariance function.
#'
#' This function computes the covariance matrix of the GP given its covariance function,
#' over a specified region.
#'
#' @param from numeric value for the left limit of the interval of interest
#' @param to numeric value for the right limit of the interval of interest
#' @param m An integer value of number of points to simulate.
#' @param K A covariance function of two arguments that is PD.
#' @return A m by m covariance matrix
#' @export
#' @examples
#' fun <- generate_K_true(1,1)
#' compute_matrix_given_cov(from = 0, to = 1, m = 5, K = fun)
compute_matrix_given_cov <- function(from, to, m, K){
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  Sigma
}
