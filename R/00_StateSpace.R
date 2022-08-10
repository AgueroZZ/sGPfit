#' Construct the joint precision matrix of the sGP process and its derivative.
#'
#' This function constructs the joint precision matrix of the sGP process and its
#' derivative (second entry) given the sGP's periodicity and standard deviation,
#' at a given set of locations.
#'
#' @param t_vec A positive numeric vector denotes the locations to evaluate the GP.
#' @param a A positive scalar represents the periodicity parameter.
#' @param sd A positive scalar represents the standard deviation parameter.
#' @return A sparse precision matrix for the sGP (and derivative) at the locations.
#' @importFrom methods as
#' @export
#' @examples
#' joint_prec_construct(t_vec = c(1,2,3), a = 1, sd = 1)
joint_prec_construct <- function(t_vec, a, sd){
  n <- length(t_vec)
  Blist <- list()
  AClist <- list()
  Clist <- list()

  ### Construct transition matrix:
  M_construct <- function(t0, t1, a){
    M <- matrix(nrow = 2, ncol = 2)
    d <- t1 - t0
    M[1,1] <- cos(a*d)
    M[1,2] <- (1/a)*sin(a*d)
    M[2,1] <- (-a)*sin(a*d)
    M[2,2] <- cos(a*d)
    M
  }

  ### Construct Sig matrix:
  Sig_construct <- function(t0, t1, a, sd){
    Sig <- matrix(nrow = 2, ncol = 2)
    d <- t1 - t0
    Sig[1,1] <- (d/(2*(a^2))) - sin(2*a*d)/(4*(a^3))
    Sig[1,2] <- (sin(a*d)^2)/(2*(a^2))
    Sig[2,1] <- Sig[1,2]
    Sig[2,2] <- (d/2) + (sin(2*a*d)/(4*a))
    (sd^2)*Sig
  }
  Compute_Ci <- function(t_vec, i){
    Ci <- Matrix::forceSymmetric(solve(Sig_construct(t0 = t_vec[i], t1 = t_vec[i+1], a, sd)))
    Ci
  }

  Compute_Ai <- function(t_vec, i, Ci) {
    Ti <- M_construct(t0 = t_vec[i], t1 = t_vec[i+1], a)
    t(Ti) %*% Ci %*% Ti
  }

  Compute_Bi <- function(t_vec, i, Ci) {
    Ti <- M_construct(t0 = t_vec[i], t1 = t_vec[i+1], a)
    -t(Ti) %*% Ci
  }

  for (i in 1:(n - 1)) {
    Clist[[i]] <- Compute_Ci(t_vec = t_vec, i = i)
  }

  for (i in 2:(n - 1)) {
    AClist[[i]] <- Compute_Ai(t_vec = t_vec, i = i, Ci = Clist[[i]]) + Clist[[i-1]]
  }

  AClist[[1]] <- Compute_Ai(t_vec = t_vec, i = 1, Ci = Clist[[1]]) + Compute_Ci(t_vec = c(0,t_vec[1]), i = 1)
  AClist[[n]] <- Compute_Ci(t_vec = t_vec, i = (n - 1))

  for (i in 1:(n - 1)) {
    Blist[[i]] <- Compute_Bi(t_vec = t_vec, i = i, Ci = Clist[[i]])
  }

  Qlist <- list()
  Q <- matrix(0, nrow = 0, ncol = n*2)
  for (i in 1:(n-1)) {
    Qlist[[i]] <- cbind(matrix(0,nrow = 2, ncol = 2 * (i-1)), AClist[[i]], Blist[[i]], matrix(0,nrow = 2, ncol = (2 * (n-i-1))) )
    Q <- rbind(Q,Qlist[[i]])
  }
  Q <- rbind(Q,cbind(matrix(0,nrow = 2, ncol = 2 * (n-1)), AClist[[n]]))
  Q <- Matrix::forceSymmetric(Q)
  as(as.matrix(Q), "dgTMatrix")
}


#' Simulate the exact sGP and its derivative at a fine mesh, using the state-space
#' representation.
#'
#' This function simulates sample path of the exact sGP with its derivatives at a fine
#' mesh, using the state-space representation.
#'
#' @param t A vector of locations to evaluate the sGP and its derivative, by default being NULL.
#' @param mesh_size A scalar denotes the size of each mesh to create, if t is NULL
#' @param max_t A scalar denotes the maximum of location to consider if t is NULL.
#' @param a A positive scalar represents the periodicity parameter.
#' @param sd A positive scalar represents the standard deviation parameter.
#' @return A dataframe that contains the sampled sGP and its derivative at the vector t.
#' @export
sim_sGP_Var <- function(t = NULL, mesh_size = 0.01, max_t = 10, a, sd){
  if(is.null(t)){
    t <- seq(0, max_t, by = mesh_size)
  }
  ### Construct transition matrix:
  M_construct <- function(t0, t1, a){
    M <- matrix(nrow = 2, ncol = 2)
    d <- t1 - t0
    M[1,1] <- cos(a*d)
    M[1,2] <- (1/a)*sin(a*d)
    M[2,1] <- (-a)*sin(a*d)
    M[2,2] <- cos(a*d)
    M
  }
  ### Construct Sig matrix:
  Sig_construct <- function(t0, t1, a, sd){
    Sig <- matrix(nrow = 2, ncol = 2)
    d <- t1 - t0
    Sig[1,1] <- (d/(2*(a^2))) - sin(2*a*d)/(4*(a^3))
    Sig[1,2] <- (sin(a*d)^2)/(2*(a^2))
    Sig[2,1] <- Sig[1,2]
    Sig[2,2] <- (d/2) + (sin(2*a*d)/(4*a))
    (sd^2)*Sig
  }
  n <- length(t) - 1
  matT <- M_construct(t0 = t[1], t1 = t[2], a = a)
  Sig <- Sig_construct(t0 = t[1], t1 = t[2], a = a, sd = sd)
  sample_path <- tsDyn::VAR.sim(B = matT, lag = 1, n = n, starting = NULL, varcov = Sig, include = "none", returnStarting = TRUE)
  result <- data.frame(t = t, sample_path)
  names(result)[2] <- "function"
  names(result)[3] <- "1stDeriv"
  result
}





