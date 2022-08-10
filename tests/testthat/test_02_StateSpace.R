## Testing for sGP computing with State Space:
library(sGPfit)

## Construct precision matrix:
Q <- joint_prec_construct(t_vec = c(0.1,0.2,0.3), a = 1, sd = 1)
Exact <- compute_matrix_given_cov(from = 0.1, to = 0.3, m = 3, K = generate_K_true(1,1))

test_that("Dimensions of created matrices are correct",{
  expect_equal(dim(Q), c(6,6))
}
)

test_that("Precision matrices have correct eigenvalues",{
  expect_equal(min(eigen(Q, only.values = T)$values)>=0, TRUE)
  expect_equal(min(eigen(Q[c(1,3,5), c(1,3,5)], only.values = T)$values)>=0, TRUE)
  expect_equal(min(eigen(Q[c(2,4,6), c(2,4,6)], only.values = T)$values)>=0, TRUE)
}
)

test_that("Covariance matrices match the exact one",{
  expect_equal(max(solve(Q)[c(1,3,5), c(1,3,5)] - Exact), 0, tolerance = 0.001)
}
)
