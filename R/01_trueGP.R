#' Generate the true covariance function of the sGP, given its periodicity parameter
#' and its SD parameter.
#'
#' This function takes the periodicity parameter and SD parameter to generate a covariance function
#' that correspond to the specified sGP.
#'
#' @param sigma The value of the SD parameter.
#' @param alpha The value of the periodicity parameter.
#' @return A Positive Definite covariance function that takes two arguments.
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
#' @param grid A vector of locations to be provided as the grid.
#' @param K A covariance function of two arguments that is PD.
#' @return A m by m covariance matrix, with m being length of grid.
#' @export
#' @examples
#' fun <- generate_K_true(1,1)
#' compute_matrix_given_cov(grid = c(0.1,0.2,0.3,0.4), K = fun)
compute_matrix_given_cov <- function(grid, K){
  t <- sort(grid)
  Cov <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  Cov
}


#' Simulates a Gaussian process with a given covariance function.
#'
#' This function simulates a Gaussian process over a fine-grid,
#'  provided with a given covariance function.
#'
#' @param grid A vector of locations to be provided as the grid.
#' @param K A covariance function of two arguments that is PD. Default being Brownian motion.
#' @return A dataframe with time index and value of the sample path.
#' @export
#' @examples
#' sim_GP(grid = c(0.1,0.2,0.3), K = function(s, t) {min(s, t)})
sim_GP <- function(grid, K = function(s, t) {min(s, t)}) {
  x <- sort(grid)
  Cov <- compute_matrix_given_cov(grid = x, K = K)
  sim_path <- MASS::mvrnorm(mu = rep(0, times = length(x)), Sigma = Cov)
  return(data.frame("x" = x, "f" = sim_path))
}

