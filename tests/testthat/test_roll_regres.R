context("Testing `roll_regres` and `roll_regres.fit`")

# set.seed(23457799)
# n <- 25
# p <- 2
# X <- cbind(1, matrix(rnorm(p * n), ncol = 2))
# dimnames(X) <- list(paste0("I", 1:n), c("(Intercept)", paste0("X", 1:p)))
# y <- drop(X %*% c(0, runif(p))) + rnorm(n)
# dput(round(X, 3))
# dput(round(y, 3))
X <- structure(c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, -0.55, 1.11, 0.26, 1.591, 0.762, 0.082,
  0.136, 0.176, -0.9, 0.489, -0.938, -0.679, -0.516, 0, -0.513,
  0.463, -2.902, 1.237, 0.896, 0.255, -0.148, 1.378, -1.957, 0.026,
  0.443, 0.503, 0.54, -1.559, -1.086, -0.795, -0.242, -0.44, -2.424,
  -1.485, -1.345, 0.855, 0.084, 1.784, 1.46, 0.207, -0.974, -0.877,
  -0.588, 0.551, 1.332, 0.641, 0.434, -0.349, -1.097, -1.832),
  .Dim = c(25L, 3L), .Dimnames = list(c(
    "I1", "I2", "I3", "I4", "I5", "I6", "I7",
    "I8", "I9", "I10", "I11", "I12", "I13", "I14", "I15", "I16",
    "I17", "I18", "I19", "I20", "I21", "I22", "I23", "I24", "I25"),
    c("(Intercept)", "X1", "X2")))
y <- c(I1 = 1.644, I2 = -0.511, I3 = -0.636, I4 = 0.596, I5 = -0.694,
       I6 = 0.568, I7 = -1.929, I8 = 0.433, I9 = -1.124, I10 = -0.094,
       I11 = -2.489, I12 = -0.571, I13 = 1.336, I14 = -1.34, I15 = 1.463,
       I16 = -0.744, I17 = -1.409, I18 = 0.429, I19 = 0.387, I20 = -0.571,
       I21 = 1.838, I22 = 1.602, I23 = -1.261, I24 = -1.073, I25 = -0.735)

test_that("`roll_regres.fit` gives the same as `roll_cpp` but with dimnames", {
  out <- roll_regres.fit(X, y, 20L)
  expect_equal(out, roll_cpp(X = X, Y = y, 20), check.attributes = FALSE)

  expect_equal(dimnames(out), dimnames(X))
})

test_that("`roll_regres` gives the same as `roll_regres.fit`", {
  df <- cbind(y = y, data.frame(X[, -1L, drop = FALSE]))
  expect_equal(roll_regres.fit(X, y, 20L), roll_regres(y ~ X1 + X2, df, 20L),
               check.attributes = FALSE)
})