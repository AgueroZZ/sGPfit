#' Create the Precision matrix of the BT approximation
#'
#' This function outputs the precision matrix of the BT approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the BT basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param accuracy A positive integer represents the integration size used to compute the inner product. Smaller value
#' means more accurate inner product with longer computational time.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the precision matrix of the weights.
#' @export
Compute_Q_Bt <- function(a,k,region, accuracy = 0.01, boundary = TRUE){
  ss <- function(M) {forceSymmetric(M + t(M))}
  x <- seq(min(region),max(region),by = accuracy)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4,
                                                     dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4))
  }
  Bmatrix <- eval.basis(x, B_basis, Lfdobj=0, returnMatrix=T)
  B1matrix <-  eval.basis(x, B_basis, Lfdobj=1, returnMatrix=T)
  B2matrix <-  eval.basis(x, B_basis, Lfdobj=2, returnMatrix=T)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- as(apply(Bmatrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B1cos <- as(apply(B1matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  B2cos <- as(apply(B2matrix, 2, function(x) x*cos_matrix), "dgCMatrix")
  Bsin <- as(apply(Bmatrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B1sin <- as(apply(B1matrix, 2, function(x) x*sin_matrix), "dgCMatrix")
  B2sin <- as(apply(B2matrix, 2, function(x) x*sin_matrix), "dgCMatrix")

  ### Compute I, L, T:
  Numerical_I <- as(diag(c(diff(c(0,x)))), "dgCMatrix")

  ### T
  T00 <- t(Bcos) %*% Numerical_I %*% Bcos
  T10 <- t(B1cos) %*% Numerical_I %*% Bcos
  T11 <- t(B1cos) %*% Numerical_I %*% B1cos
  T20 <- t(B2cos) %*% Numerical_I %*% Bcos
  T21 <- t(B2cos) %*% Numerical_I %*% B1cos
  T22 <- t(B2cos) %*% Numerical_I %*% B2cos

  ### L
  L00 <- t(Bsin) %*% Numerical_I %*% Bsin
  L10 <- t(B1sin) %*% Numerical_I %*% Bsin
  L11 <- t(B1sin) %*% Numerical_I %*% B1sin
  L20 <- t(B2sin) %*% Numerical_I %*% Bsin
  L21 <- t(B2sin) %*% Numerical_I %*% B1sin
  L22 <- t(B2sin) %*% Numerical_I %*% B2sin

  ### I
  I00 <- t(Bsin) %*% Numerical_I %*% Bcos
  I10 <- t(B1sin) %*% Numerical_I %*% Bcos
  I11 <- t(B1sin) %*% Numerical_I %*% B1cos
  I20 <- t(B2sin) %*% Numerical_I %*% Bcos
  I21 <- t(B2sin) %*% Numerical_I %*% B1cos
  I22 <- t(B2sin) %*% Numerical_I %*% B2cos

  ## G = <phi,phj>
  G <- rbind(cbind(T00, t(I00)), cbind(I00,L00))

  ## C = <D^2phi,D^2phj>
  C11 <- T22 - 2*a*ss(I21) - (a^2)*ss(T20) + 2*(a^3)*ss(I10) + 4 * (a^2) * L11 + (a^4)*T00
  C22 <- L22 + 2*a*ss(I21) - (a^2)*ss(L20) - 2*(a^3)*ss(I10) + 4 * (a^2) * T11 + (a^4)*L00
  C12 <- I22 + 2*a*T21 - (a^2)* ss(I20) - 2*a*t(L21) - 4*(a^2)*I11 + 2*(a^3)*L10 - 2*(a^3)*t(T10) + (a^4)*I00
  C <- rbind(cbind(C11,C12), cbind(t(C12), C22))

  ## M = <phi,D^2phj>
  M11 <- t(T20) - (2*a)*t(I10) - (a^2)*T00
  M12 <- t(I20) + (2*a)*t(T10) - (a^2)*I00
  M21 <- t(I20) - (2*a)*t(L10) - (a^2)*I00
  M22 <- t(L20) + (2*a)*t(I10) - (a^2)*L00
  M <- rbind(cbind(M11,M12), cbind(M21, M22))


  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)
  forceSymmetric(Q)
}


#' Create the Design matrix of the BT approximation
#'
#' This function outputs the design matrix of the BT approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the BT basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param accuracy A positive integer represents the integration size used to compute the inner product. Smaller value
#' means more accurate inner product with longer computational time.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the design matrix of the weights.
#' @export
Compute_B_Bt <- function(x, a, k, region, boundary = T){
  if(boundary){
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4,
                                                     dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4))
  }
  Bmatrix <- eval.basis(x, B_basis, Lfdobj=0, returnMatrix=T)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- apply(Bmatrix, 2, function(x) x*cos_matrix)
  Bsin <- apply(Bmatrix, 2, function(x) x*sin_matrix)
  cbind(Bcos, Bsin)
}



#' Sample from the approximate prior using BT splines
#'
#' This function samples sample paths from the approximate prior using BT splines.
#'
#' @param x A vector that specifies where to evaluate the sample path.
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the BT basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param accuracy A positive integer represents the integration size used to compute the inner product. Smaller value
#' means more accurate inner product with longer computational time.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @param n The number of samples to draw.
#' @return A matrix with each column denote a sample path f(x).
#' @export
sampling_from_weights <- function(x, a, k, region, boundary = T, n = 1){
  Prec <- Compute_Q_Bt(a, k, region, boundary = boundary)
  B <- Compute_B_Bt(x, a, k, region, boundary = boundary)
  if(boundary){
    coefs_samps <- rmvnp(n = n, mu = rep(0,(2*(k-2))), Omega = as.matrix(Prec))
  }
  else{
    coefs_samps <- rmvnp(n = n, mu = rep(0,(2*k)), Omega = as.matrix(Prec))
  }
  splfd <- B %*% t(coefs_samps)
  (splfd)
}













