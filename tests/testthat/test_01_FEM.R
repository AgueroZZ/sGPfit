## Testing for sGP approximation with sB splines:
library(sGPfit)

## Construct precision matrix:
Q <- Compute_Q_sB(a = 1, k = 5, region = c(0,10), boundary = TRUE)
Q2 <- Compute_Q_sB(a = 1, k = 5, region = c(0,1), boundary = FALSE)
## Construct design matrix:
B <- Compute_B_sB(x = c(0.1,0.5,0.8), a = 1, k = 5, region = c(0,1))
B2 <- Compute_B_sB(x = c(0.1,0.5,0.8), a = 1, k = 5, region = c(0,1), boundary = FALSE)

test_that("Dimensions of created matrices are correct",{
  expect_equal(dim(Q), c(6,6))
  expect_equal(dim(Q2), c(10,10))
  expect_equal(dim(B), c(3,6))
  expect_equal(dim(B2), c(3,10))
}
)

test_that("Classes of created matrices are correct",{
  expect_equal(class(Q)[[1]], "dsCMatrix")
  expect_equal(class(Q2)[[1]], "dsCMatrix")
  expect_equal(class(B)[[1]], "matrix")
  expect_equal(class(B2)[[1]], "matrix")
}
)

test_that("Precision matrices have correct eigenvalues",{
  expect_equal(min(eigen(Q, only.values = T)$values)>=0, TRUE)
  expect_equal(min(eigen(Q2, only.values = T)$values), 0, tolerance = 0.0001)
}
)

