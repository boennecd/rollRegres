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
  expect_equal(
    out, roll_cpp(X = X, Y = y, 20, FALSE, FALSE, FALSE, 1:nrow(X), FALSE,
                  TRUE),
    check.attributes = FALSE)

  expect_equal(dimnames(out$coefs), dimnames(X))
})

test_that("`roll_regres.fit` gives the same as `roll_cpp` but with dimnames when `do_compute` arg is used", {
  out <- roll_regres.fit(X, y, 20L, do_compute = "sigmas")
  expect_equal(
    out$sigmas,
    drop(roll_cpp(
      X = X, Y = y, 20, do_compute_sigmas = TRUE,
      do_compute_R_sqs = FALSE, do_1_step_forecasts = FALSE,
      grp = 1:nrow(X), use_grp = FALSE, do_downdates = TRUE)$sigmas),
    check.attributes = FALSE)
  expect_equal(names(out$sigmas), row.names(X))

  out <- roll_regres.fit(X, y, 20L, do_compute = "r.squareds")
  expect_equal(
    out$r.squareds,
    drop(roll_cpp(
      X = X, Y = y, 20, do_compute_sigmas = FALSE,
      do_compute_R_sqs = TRUE, do_1_step_forecasts = FALSE,
      grp = 1:nrow(X), use_grp = FALSE, do_downdates = TRUE)$r.squareds),
    check.attributes = FALSE)
  expect_equal(names(out$r.squareds), row.names(X))

  out <- roll_regres.fit(X, y, 20L, do_compute = "1_step_forecasts")
  expect_equal(
    out$one_step_forecasts,
    drop(roll_cpp(
      X = X, Y = y, 20, do_compute_sigmas = FALSE,
      do_compute_R_sqs = FALSE, do_1_step_forecasts = TRUE,
      grp = 1:nrow(X), use_grp = FALSE,
      do_downdates = TRUE)$one_step_forecasts),
    check.attributes = FALSE)
  expect_equal(names(out$one_step_forecasts), row.names(X))
})

test_that("`roll_regres` gives the same as `roll_regres.fit`", {
  df <- cbind(y = y, data.frame(X[, -1L, drop = FALSE]))
  expect_equal(roll_regres.fit(X, y, 20L), roll_regres(y ~ X1 + X2, df, 20L),
               check.attributes = FALSE)
})

test_that("`roll_regres` gives the same as `roll_regres.fit` with `do_compute` arg", {
  df <- cbind(y = y, data.frame(X[, -1L, drop = FALSE]))
  do_comp <- c("sigmas", "r.squareds", "1_step_forecasts")
  expect_equal(
    roll_regres.fit(X, y, 20L, do_compute = do_comp),
    roll_regres(y ~ X1 + X2, df, 20L, do_compute = do_comp),
    check.attributes = FALSE)

  # see github.com/boennecd/rollRegres/issues/1
  expect_equal(
    roll_regres.fit(X, y, 20L, do_compute = do_comp, do_downdates = FALSE),
    roll_regres(y ~ X1 + X2, df, 20L, do_compute = do_comp,
                do_downdates = FALSE),
    check.attributes = FALSE)
})

test_that("`roll_regres.fit` post warning when low p compared to n", {
  # see github.com/boennecd/rollRegres/issues/2
  x <- matrix(rnorm(4 * 20), ncol = 4)
  y <- rnorm(20)

  expect_warning(
    roll_regres.fit(x, y, width = 11L),
    "low sample size relative to number of parameters")

  # should have no warning
  roll_regres.fit(x, y, width = 12L)
})


test_that("`roll_regres.fit` works as expected with 'min_obs'", {
  roll_regress_R_for_loop <- function(X, y, width, grp, downdate, min_obs){
    grp <- grp + 1L - min(grp)
    u_grp = unique(grp)
    n <- nrow(X)
    p <- ncol(X)
    out <- matrix(NA_real_, n, p)
    sigmas             <- rep(NA_real_, n)
    r.squared          <- rep(NA_real_, n)

    start_val <- max(which(u_grp <= width))
    for(g in u_grp[start_val:length(u_grp)]){
      idx <-
        if(downdate)
          which(grp %in% (g - width + 1L):g) else
            which(grp %in% 1:g)
      i <- which(grp == g)
      if(length(idx) < min_obs)
        next
      fit <- lm(y[idx] ~ -1 + X[idx, , drop = FALSE])
      out[i, ] <- sapply(fit$coefficients, rep, times = length(i))

      su <- summary(fit)
      sigmas[i] <- su$sigma

      ss1 <- sum((y[idx] - mean(y[idx]))^2)
      ss2 <- sum(fit$residuals^2)
      r.squared[i] <- 1 - ss2 / ss1
    }

    list(coef = out, sigmas = sigmas, r.squared = r.squared,
         one_step_forecasts = NULL)
  }

  #####
  # simulate complete data
  n   <- 2L * 12L * 21L # x years w/ 12 months of 21 trading days
  mth <- (seq_len(n) - 1L) %/% 21L + 1L # group by months
  set.seed(29478439)
  X <- matrix(rnorm(n * 4L), ncol = 4L,
              dimnames = list(1:n, paste0("X", 1:4)))
  y <- rnorm(n)

  for(i in 1:20){
    keep <- seq_along(y) %in% sample.int(nrow(X), as.integer(n * .5))
    x1   <- X  [keep, ]
    y1   <- y  [keep]
    mth1 <- mth[keep]

    o1 <- roll_regress_R_for_loop(X = x1, y = y1, grp = mth1, width = 6L,
                                  downdate = TRUE, min_obs = 63L)

    o2 <- roll_regres.fit(
      x = x1, y = y1, width = 6L, grp = mth1, do_downdates = TRUE,
      min_obs = 63L,
      do_compute = c("sigmas", "r.squareds"))
    expect_equal(o1, o2, check.attributes = FALSE)
  }
})

