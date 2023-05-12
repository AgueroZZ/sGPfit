#' Create the Precision matrix of the sB approximation
#'
#' This function outputs the precision matrix of the sB approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the sB basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param accuracy A positive integer represents the integration size used to compute the inner product. Smaller value
#' means more accurate inner product with longer computational time.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the precision matrix of the weights.
#' @importFrom methods as
#' @import Matrix
#' @export
#' @examples
#' Compute_Q_sB(a = 1, k = 5, region = c(0,1))
Compute_Q_sB <- function(a,k,region, accuracy = 0.01, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + t(M))}
  x <- seq(min(region),max(region),by = accuracy)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  B1matrix <-  fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE)
  B2matrix <-  fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE)
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


  ### Inner product involving B spline:
  Bmatrix <- as(Bmatrix, "dgCMatrix")
  B1matrix <- as(B1matrix, "dgCMatrix")
  B2matrix <- as(B2matrix, "dgCMatrix")

  BB <- t(Bmatrix) %*% Numerical_I %*% Bmatrix
  B2B2 <- t(B2matrix) %*% Numerical_I %*% B2matrix
  BB2 <- t(Bmatrix) %*% Numerical_I %*% B2matrix
  BS <- t(Bmatrix) %*% Numerical_I %*% Bsin
  BC <- t(Bmatrix) %*% Numerical_I %*% Bcos
  BS1 <- t(Bmatrix) %*% Numerical_I %*% B1sin
  BC1 <- t(Bmatrix) %*% Numerical_I %*% B1cos
  BS2 <- t(Bmatrix) %*% Numerical_I %*% B2sin
  BC2 <- t(Bmatrix) %*% Numerical_I %*% B2cos
  B2S <- t(B2matrix) %*% Numerical_I %*% Bsin
  B2C <- t(B2matrix) %*% Numerical_I %*% Bcos
  B2S1 <- t(B2matrix) %*% Numerical_I %*% B1sin
  B2C1 <- t(B2matrix) %*% Numerical_I %*% B1cos
  B2S2 <- t(B2matrix) %*% Numerical_I %*% B2sin
  B2C2 <- t(B2matrix) %*% Numerical_I %*% B2cos

  ## G = <phi,phj>
  G <- rbind(cbind(T00, t(I00), t(BC)), cbind(I00,L00, t(BS)), cbind(BC,BS,BB))

  ## C = <D^2phi,D^2phj>
  C11 <- T22 - 2*a*ss(I21) - (a^2)*ss(T20) + 2*(a^3)*ss(I10) + 4 * (a^2) * L11 + (a^4)*T00
  C22 <- L22 + 2*a*ss(I21) - (a^2)*ss(L20) - 2*(a^3)*ss(I10) + 4 * (a^2) * T11 + (a^4)*L00
  C12 <- I22 + 2*a*T21 - (a^2)* ss(I20) - 2*a*t(L21) - 4*(a^2)*I11 + 2*(a^3)*L10 - 2*(a^3)*t(T10) + (a^4)*I00

  C13 <- t(B2C2) - 2*a*t(B2S1) - (a^2)*t(B2C)
  C23 <- t(B2S2) + 2*a*t(B2C1) - (a^2)*t(B2S)
  C33 <- B2B2
  C <- rbind(cbind(C11,C12,C13), cbind(t(C12), C22, C23), cbind(t(C13), t(C23), C33))


  ## M = <phi,D^2phj>
  M11 <- t(T20) - (2*a)*t(I10) - (a^2)*T00
  M12 <- t(I20) + (2*a)*t(T10) - (a^2)*I00
  M21 <- t(I20) - (2*a)*t(L10) - (a^2)*I00
  M22 <- t(L20) + (2*a)*t(I10) - (a^2)*L00

  M13 <- t(B2C)
  M23 <- t(B2S)
  M31 <- BC2 - (2*a)*BS1 - (a^2)*BC
  M32 <- BS2 + (2*a)*BC1 - (a^2)*BS
  M33 <- BB2

  M <- rbind(cbind(M11,M12,M13), cbind(M21,M22,M23), cbind(M31,M32,M33))


  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)
  Matrix::forceSymmetric(Q)
}




#' Create the Precision matrix of the sB approximation (without the B spline)
#'
#' This function outputs the precision matrix of the sB approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the sB basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param accuracy A positive integer represents the integration size used to compute the inner product. Smaller value
#' means more accurate inner product with longer computational time.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the precision matrix of the weights.
#' @importFrom methods as
#' @import Matrix
#' @export
#' @examples
#' Compute_Q_sB_old(a = 1, k = 5, region = c(0,1))
Compute_Q_sB_old <- function(a,k,region, accuracy = 0.01, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + t(M))}
  x <- seq(min(region),max(region),by = accuracy)
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4,
                                                     dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  B1matrix <-  fda::eval.basis(x, B_basis, Lfdobj=1, returnMatrix=TRUE)
  B2matrix <-  fda::eval.basis(x, B_basis, Lfdobj=2, returnMatrix=TRUE)
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
  Matrix::forceSymmetric(Q)
}



#' Create the Precision matrix of the B spline approximation, using Richardson Method
#'
#' This function outputs the precision matrix of the B spline approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the B spline basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the precision matrix of the weights.
#' @importFrom methods as
#' @import Matrix
#' @export
#' @examples
#' Compute_Q_BR(a = 1, k = 5, region = c(0,1))
Compute_Q_BR <- function(a,k,region, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + t(M))}
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }

  ### Compute G,C,M:
  G <- fda::inprod(fdobj1 = B_basis, fdobj2 = B_basis, Lfdobj1 = 0, Lfdobj2 = 0)
  C <- fda::inprod(fdobj1 = B_basis, fdobj2 = B_basis, Lfdobj1 = 2, Lfdobj2 = 2)
  M <- fda::inprod(fdobj1 = B_basis, fdobj2 = B_basis, Lfdobj1 = 2, Lfdobj2 = 0)

  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)
  Matrix::forceSymmetric(Q)
}


#' Create the Precision matrix of the B spline approximation with grid approximation
#'
#' This function outputs the precision matrix of the B spline approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the B spline basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param accuracy A positive integer represents the integration size used to compute the inner product. Smaller value
#' means more accurate inner product with longer computational time.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the precision matrix of the weights.
#' @importFrom methods as
#' @import Matrix
#' @export
#' @examples
#' Compute_Q_B_grid(a = 1, k = 5, region = c(0,1))
Compute_Q_B_grid <- function(a,k,region, boundary = TRUE, accuracy = 0.01){
  x <- seq(min(region),max(region),by = accuracy)
  ss <- function(M) {Matrix::forceSymmetric(M + t(M))}
  if(boundary == TRUE){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }

  Bmatrix <- eval.basis(x, B_basis, Lfdobj=0, returnMatrix=T)
  B2matrix <-  eval.basis(x, B_basis, Lfdobj=2, returnMatrix=T)
  Numerical_I <- as(diag(c(diff(c(0,x)))), "dgCMatrix")


  ### Compute G,C,M:
  G <- t(Bmatrix) %*% Numerical_I %*% Bmatrix
  C <- t(B2matrix) %*% Numerical_I %*% B2matrix
  M <- t(Bmatrix) %*% Numerical_I %*% B2matrix

  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)
  Matrix::forceSymmetric(Q)
}




#' Create the Precision matrix of the B spline approximation
#'
#' This function outputs the precision matrix of the B spline approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the B spline basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the precision matrix of the weights.
#' @importFrom methods as
#' @import Matrix
#' @export
#' @examples
#' Compute_Q_B(a = 1, k = 5, region = c(0,1))
Compute_Q_B <- function(a,k,region, boundary = TRUE){
  ss <- function(M) {Matrix::forceSymmetric(M + t(M))}
  B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                        nbasis = k,
                                                        norder = 4))
  B_basis_fd <- fd(coef = diag(1,k), basisobj = B_basis)
  ### Compute G,C,M:
  G <- fda::inprod.bspline(fdobj1 = B_basis_fd, fdobj2 = B_basis_fd, nderiv1 = 0, nderiv2 = 0)
  C <- fda::inprod.bspline(fdobj1 = B_basis_fd, fdobj2 = B_basis_fd, nderiv1 = 2, nderiv2 = 2)
  M <- fda::inprod.bspline(fdobj1 = B_basis_fd, fdobj2 = B_basis_fd, nderiv1 = 2, nderiv2 = 0)

  ### Compute the final precision matrix: Q
  Q <- (a^4)*G + C + (a^2)*ss(M)

  if(boundary == TRUE){
    Q <- Q[-c(1,2),-c(1,2)]
  }
  Q
}



#' Create the Design matrix of the sB approximation
#'
#' This function outputs the design matrix of the sB approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param x A vector of numeric values that denotes the evaluation locations.
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the sB basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the design matrix of the weights.
#' @export
#' @examples
#' Compute_B_sB(x = c(1,2,3), a = 1, k = 5, region = c(0,3))
Compute_B_sB <- function(x, a, k, region, boundary = TRUE){
  if(boundary){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- apply(Bmatrix, 2, function(x) x*cos_matrix)
  Bsin <- apply(Bmatrix, 2, function(x) x*sin_matrix)
  cbind(Bcos, Bsin,Bmatrix)
}




#' Create the Design matrix of the sB approximation (without the B spline)
#'
#' This function outputs the design matrix of the sB approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param x A vector of numeric values that denotes the evaluation locations.
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the sB basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the design matrix of the weights.
#' @export
#' @examples
#' Compute_B_sB_old(x = c(1,2,3), a = 1, k = 5, region = c(0,3))
Compute_B_sB_old <- function(x, a, k, region, boundary = TRUE){
  if(boundary){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4,
                                                     dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                     nbasis = k,
                                                     norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  cos_matrix <- cos(a*x)
  sin_matrix <- sin(a*x)
  Bcos <- apply(Bmatrix, 2, function(x) x*cos_matrix)
  Bsin <- apply(Bmatrix, 2, function(x) x*sin_matrix)
  cbind(Bcos, Bsin)
}




#' Create the Design matrix of the B spline approximation
#'
#' This function outputs the design matrix of the B spline approximation as a sparse matrix, given the number of knots,
#' the periodicity parameter a, and the region of interest.
#'
#' @param x A vector of numeric values that denotes the evaluation locations.
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the B spline basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @return A sparse matrix as the design matrix of the weights.
#' @export
#' @examples
#' Compute_B_B(x = c(1,2,3), a = 1, k = 5, region = c(0,3))
Compute_B_B <- function(x, a, k, region, boundary = TRUE){
  if(boundary){
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4,
                                                          dropind = c(1,2)))
  }
  else{
    B_basis <- suppressWarnings(fda::create.bspline.basis(rangeval = c(min(region),max(region)),
                                                          nbasis = k,
                                                          norder = 4))
  }
  Bmatrix <- fda::eval.basis(x, B_basis, Lfdobj=0, returnMatrix=TRUE)
  return(Bmatrix)
}





#' Sample from the approximate prior using sB splines
#'
#' This function samples sample paths from the approximate prior using sB splines.
#'
#' @param x A vector that specifies where to evaluate the sample path.
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the sB basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @param n The number of samples to draw.
#' @return A matrix with each column denote a sample path f(x).
#' @export
#' @examples
#' sampling_from_sB(x = c(1,2,3), a = 1, k = 5, region = c(0,5))
sampling_from_sB <- function(x, a, k, region, boundary = TRUE, n = 1){
  Prec <- Compute_Q_sB(a, k, region, boundary = boundary)
  B <- Compute_B_sB(x, a, k, region, boundary = boundary)
  if(boundary){
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,(2*(k-2))), Omega = as.matrix(Prec))
  }
  else{
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,(2*k)), Omega = as.matrix(Prec))
  }
  splfd <- B %*% t(coefs_samps)
  (splfd)
}



#' Sample from the approximate prior using B splines
#'
#' This function samples sample paths from the approximate prior using B splines.
#'
#' @param x A vector that specifies where to evaluate the sample path.
#' @param a A positive scalar represents the periodicity parameter.
#' @param k A positive integer represents the number of knots used to define the B basis. The number of
#' basis functions equals to 2 times k or 2 times (k-1) if boundary is TRUE.
#' @param region A vector of size 2 that denotes the upper and lower interval limit of the region of interest.
#' @param boundary A logical value to denote whether to consider the boundary conditions.
#' @param n The number of samples to draw.
#' @return A matrix with each column denote a sample path f(x).
#' @export
#' @examples
#' sampling_from_B(x = c(1,2,3), a = 1, k = 5, region = c(0,5))
sampling_from_B <- function(x, a, k, region, boundary = TRUE, n = 1){
  Prec <- Compute_Q_B(a, k, region, boundary = boundary)
  B <- Compute_B_B(x, a, k, region, boundary = boundary)
  if(boundary){
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,((k-2))), Omega = as.matrix(Prec))
  }
  else{
    coefs_samps <- LaplacesDemon::rmvnp(n = n, mu = rep(0,(k)), Omega = as.matrix(Prec))
  }
  splfd <- B %*% t(coefs_samps)
  (splfd)
}













