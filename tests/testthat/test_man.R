context("Testing man files gives same results")

test_that("'roll_regres' gives the same", {
  set.seed(29132867)
  n <- 50
  p <- 2
  X <- cbind(1, matrix(rnorm(p * n), ncol = p))
  y <- drop(X %*% c(1, -1, 1)) + rnorm(n)
  df <- data.frame(y, X[, -1])

  expect_known_value(
    roll_regres(y ~ X1 + X2, df, width = 45L), file = "roll_regres_man_1.RDS")

  out <- roll_regres(
    y ~ X1 + X2, df, width = 45L,
    do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
  expect_known_value(lapply(out, tail), file = "roll_regres_man_2.RDS")
})

test_that("'roll_regres.fit' gives the same", {
  set.seed(9623556)
  n <- 50
  p <- 2
  X <- cbind(1, matrix(rnorm(p * n), ncol = p))
  y <- drop(X %*% c(1, -1, 1)) + rnorm(n)

  expect_known_value(
    roll_regres.fit(x = X, y = y, width = 45L),
    file = "roll_regres.fit_man_1.RDS")
})
