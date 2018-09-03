context("Testing `roll_cpp`")

test_that("`roll_cpp` works in 1 length block case with and without downdating", {
  # set.seed(101)
  # n <- 100
  # p <- 3
  # X <- matrix(rnorm(p * n), n, p)
  # y <- drop(X %*% runif(p)) + rnorm(n)
  # dput(round(X, 3))
  # dput(round(y, 3))
  wdth <- 40
  X <- structure(c(
    -0.326, 0.552, -0.675, 0.214, 0.311, 1.174, 0.619,
    -0.113, 0.917, -0.223, 0.526, -0.795, 1.428, -1.467, -0.237,
    -0.193, -0.85, 0.058, -0.818, -2.05, -0.164, 0.709, -0.268, -1.464,
    0.744, -1.41, 0.467, -0.119, 0.467, 0.498, 0.895, 0.279, 1.008,
    -2.073, 1.19, -0.724, 0.168, 0.92, -1.672, 0.448, 0.482, 0.758,
    -2.319, -0.46, -1.105, 0.403, 0.569, -0.706, -0.29, -1.484, -1.15,
    -0.274, 0.578, -1.397, 0.749, -1.051, 0.165, 1.13, 1.174, -0.428,
    -0.26, -1.411, -0.641, 0.112, 0.423, 0.387, -0.688, 0.149, -0.058,
    -0.075, 1.51, 1.62, 1.153, -0.078, -1.819, -1.037, 0.302, -1.278,
    0.138, -0.051, 1.852, 1.112, -0.511, -0.544, -1.729, 0.471, 0.005,
    1.348, 0.724, 1.553, 1.325, -0.034, -0.361, -0.72, 0.282, -0.791,
    -0.445, 1.365, 0.497, -0.814, 0.268, -0.592, 2.133, 1.173, 0.747,
    -0.231, 0.088, -2.184, -0.467, 1.686, -0.568, -0.047, -0.157,
    1.602, 0.769, -0.772, -0.631, -0.83, -0.591, 0.981, -0.662, -0.772,
    -2.018, -0.534, 0.435, -0.771, -0.754, -0.299, 1.664, -1.244,
    -0.783, 0.245, -0.144, -1.609, 0.952, -1.819, 1.784, 1.887, 1.491,
    -0.381, -0.909, -0.338, -1.412, 0.218, 0.67, -0.288, 0.469, -0.47,
    -0.239, -0.447, -0.619, 0.253, -0.753, 0.732, -0.403, -2.823,
    0.463, 2.133, -0.27, 0.249, 0.038, 0.394, -1.504, -1.587, -0.927,
    0.776, -0.781, -1.279, -0.001, -1.851, 0.452, -0.433, 0.714,
    0.961, 0.382, 1.218, -0.017, -0.038, 1.244, -0.956, 0.915, -0.939,
    0.112, 0.553, 0.532, -0.874, -0.187, -0.214, -0.204, 1.72, 0.202,
    0.513, 1.452, 0.364, -0.876, -0.015, -0.724, 1.969, -0.536, -0.026,
    -0.164, -1.383, 0.424, -0.79, 1.21, 0.895, -0.101, 0.297, 0.197,
    -0.157, 1.537, -2.168, 0.598, 0.043, 1.295, 0.706, 0.346, -0.08,
    0.455, 1.276, 1.265, 0.269, -0.121, 0.795, -0.514, -0.407, 1.22,
    0.084, 0.59, -0.517, 0.769, 0.802, -0.697, 1.178, 0.586, -0.467,
    0.386, -0.535, 1.057, -0.206, 0.607, -0.548, -2.1, 0.251, -0.055,
    -0.66, -1.456, 0.024, 0.548, -0.809, -0.239, -0.353, 0.82, -0.345,
    -1.183, -1.031, -0.075, 0.829, -1.036, -0.147, -0.281, -1.374,
    1.558, -0.575, 2.187, 0.794, 0.203, -0.022, 0.305, -1.109, 0.765,
    -0.022, -0.904, 0.4, -1.149, 0.188, 0.217, -1.049, -0.075, -1.568,
    -1.174, -2.148, 0.342, 0.905, 1.096, -1.471, -0.281, 0.846, -1.286,
    -0.312, -0.363, 1.412, -0.255, 0.387, 0.525, 0.615, -0.333, 0.745,
    -0.311, -0.234), .Dim = c(100L, 3L))
  y <- c(
    -0.162, -1.562, 0.805, -0.611, 0.787, -1.146, -0.569, -0.914,
    0.372, 0.398, 1.033, 0.009, 0.439, 0.967, 0.249, -0.978, -0.853,
    -2.209, -0.721, 1.521, 2.914, -2.086, -4.518, 0.653, 0.802, -1.308,
    -0.759, -1.235, 2.379, -0.188, 1.192, 0.342, 2.513, -0.145, 0.874,
    -0.492, 0.423, 1.248, 0.445, -0.459, 0.19, 1.269, -0.611, -0.167,
    1.112, -0.422, 1.795, -1.979, -1.497, -1.009, -2.166, -0.41,
    0.65, 2.016, 1.24, -0.518, -0.564, 0.123, -0.849, -1.199, -1.076,
    0.539, -0.558, -1.474, 0.691, 1.205, 0.859, -0.564, 0.152, -2.522,
    1.242, 0.891, -0.545, 0.936, 0.286, 0.032, 0.748, -1.623, -0.145,
    -2.159, 0.752, -2.578, 1.022, 1.759, 0.59, -2.389, -0.953, 2.007,
    1.262, -0.16, -0.349, 0.617, 0.445, -1.548, -1.48, -1.067, -0.739,
    1.373, -0.834, -0.656)
  wdth = 40

  roll_regress_R_for_loop <- function(X, y, width, downdate){
    n <- nrow(X)
    p <- ncol(X)
    out <- matrix(NA_real_, n, p)
    sigmas             <- rep(NA_real_, n)
    r.squared          <- rep(NA_real_, n)
    one_step_forecasts <- rep(NA_real_, n)

    for(i in width:n){
      idx <- if(downdate) (i - width + 1L):i else 1:i
      fit <- lm(y[idx] ~ -1 + X[idx, , drop = FALSE])
      out[i, ] <- fit$coefficients

      su <- summary(fit)
      sigmas[i] <- su$sigma

      ss1 <- sum((y[idx] - mean(y[idx]))^2)
      ss2 <- sum(fit$residuals^2)
      r.squared[i] <- 1 - ss2 / ss1

      if(i < n){
        next_i <- i + 1L
        one_step_forecasts[next_i] <- fit$coefficients %*% X[next_i, ]
      }
    }

    list(coef = out, sigmas = sigmas, r.squared = r.squared,
         one_step_forecasts = one_step_forecasts)
  }

  . <- function(
    do_compute_R_sqs, do_compute_sigmas, do_1_step_forecasts, use_grp,
    do_downdates)
    substitute({
      cpp_out <- roll_cpp(Y = y, X = X, window = wdth,
               do_compute_R_sqs = do_compute_R_sqs,
               do_compute_sigmas = do_compute_sigmas,
               do_1_step_forecasts = do_1_step_forecasts, grp = 1:nrow(X),
               use_grp = use_grp, do_downdates = do_downdates)
      expect_equal(r_out$coef, cpp_out$coefs)
      t1
      t2
      t3
    }, list(
      t1 = if(do_compute_sigmas)
        quote(expect_equal(r_out$sigmas, drop(cpp_out$sigmas))) else
           quote(expect_null(cpp_out$sigmas)),
      t2 = if(do_compute_R_sqs)
        quote(expect_equal(r_out$r.squared, drop(cpp_out$r.squareds))) else
          quote(expect_null(cpp_out$r.squareds)),
      t3 = if(do_1_step_forecasts)
        quote(expect_equal(
          r_out$one_step_forecasts, drop(cpp_out$one_step_forecasts))) else
            quote(expect_null(cpp_out$one_step_forecasts)),
      do_compute_R_sqs = do_compute_R_sqs,
      do_compute_sigmas = do_compute_sigmas,
      do_1_step_forecasts = do_1_step_forecasts,
      use_grp = use_grp, do_downdates = do_downdates))

  vals <- expand.grid(
    do_compute_R_sqs    = c(T, F),
    do_compute_sigmas   = c(T, F),
    do_1_step_forecasts = c(T, F),
    use_grp             = c(T, F),
    do_downdates        =   T)

  r_out <- roll_regress_R_for_loop(X, y, wdth, downdate = TRUE)
  for(i in 1:nrow(vals))
    eval(do.call(., as.list(vals[i, ])))

  vals$do_downdates <- FALSE
  r_out <- roll_regress_R_for_loop(X, y, wdth, downdate = FALSE)
  for(i in 1:nrow(vals))
    eval(do.call(., as.list(vals[i, ])))
})

test_that("`roll_cpp` works in n > 1 length block case with and without downdating", {
  # set.seed(71336382)
  # week <- as.integer(gl(25, 5))
  # week <- week[!week %in% c(3, 10:12, 18)] # miss some weeks
  # n <- length(week)
  # p <- 2
  # X <- matrix(rnorm(p * n), n, p)
  # y <- drop(X %*% runif(p)) + rnorm(n)
  # dput(week)
  # dput(round(X, 3))
  # dput(round(y, 3))
  week <- c(
    1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 4L, 4L, 4L, 4L, 4L,
    5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 8L,
    8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L, 13L, 13L, 13L, 13L, 13L,
    14L, 14L, 14L, 14L, 14L, 15L, 15L, 15L, 15L, 15L, 16L, 16L, 16L,
    16L, 16L, 17L, 17L, 17L, 17L, 17L, 19L, 19L, 19L, 19L, 19L, 20L,
    20L, 20L, 20L, 20L, 21L, 21L, 21L, 21L, 21L, 22L, 22L, 22L, 22L,
    22L, 23L, 23L, 23L, 23L, 23L, 24L, 24L, 24L, 24L, 24L, 25L, 25L,
    25L, 25L, 25L)
  X <- structure(c(
    -0.462, 0.247, -0.181, -0.66, 0.01, -1.595, -0.936,
    1.845, 0.054, 1.41, -0.094, -0.895, 0.564, 0.568, 1.929, 1.153,
    0.147, -1.213, -1.042, 0.074, -0.992, -0.261, -0.199, -0.574,
    -0.419, 0.362, 2.554, -0.611, 0.571, 0.872, 0.405, -0.133, 0.526,
    -2.09, -1.467, -1.533, 0.894, 0.645, -0.011, 1.98, -0.399, -0.977,
    -0.419, -1.165, -0.381, -0.247, -0.941, -0.613, 0.831, 1.814,
    0.672, -1.285, 0.013, -0.728, 1.329, 0.042, 0.067, -0.132, -2.167,
    0.075, 0.148, -0.122, -0.44, 0.931, -0.199, 0.406, 0.599, -0.385,
    0.582, -2.315, 0.512, 0.817, 0.041, -0.798, -0.555, -0.13, 0.164,
    0.366, -0.035, 0.21, -1.349, -0.479, 1.536, -1.405, 0.508, 0.768,
    -0.662, 1.515, 0.008, -0.374, 0.342, 0.31, 0.402, -0.556, -1.977,
    1.239, 0.164, 0.655, -0.838, -0.087, -0.72, 0.048, 0.319, 1.711,
    -1.013, 1.737, 0.359, 1.205, -0.04, -0.308, -0.669, 0.42, 0.115,
    -0.169, 0.205, -2.057, -1.172, -1.229, -1.139, 0.842, 1.43, 1.038,
    -2.857, 0.039, 1.257, -0.198, 0.872, -0.283, -1.407, -0.29, -1.724,
    -0.081, -0.978, 0.223, 0.312, 2.32, 0.18, 1.344, -0.025, -0.044,
    0.573, -0.153, -0.618, -1.644, 0.754, 1.541, -1.089, 0.69, 0.498,
    0.686, 0.57, -1.445, -0.042, -1.054, 0.489, -1.188, 0.911, 1.579,
    1.666, -0.276, 0.331, -3.046, -0.017, 0.435, 0.748, -0.952, -0.188,
    0.362, 0.434, 0.966, -0.729, 1.168, 0.595, 1.109, 0.682, -0.74,
    0.455, 0.744, 0.091, -1.292, 0.286, -0.806, -1.724, -0.245, 1.054,
    -1.315, -0.707, 0.593, 0.613, -0.895, -1.54, 0.52, -1.006, 0.101,
    -0.278, -1.573, -0.363, -0.427, 1.654, -0.452), .Dim = c(100L, 2L))
  y <- c(
    -1.053, -1.138, -1.87, -0.279, 0.419, 0.612, 0.746, -0.163,
    -0.641, 1.512, -1.285, -0.206, 0.127, 1.439, 3.36, -0.293, 0.344,
    -0.551, -2.987, -1.526, 0.607, -0.247, -2.639, -0.573, 2.146,
    1.873, 4.497, -0.999, -0.55, 0.994, 1.116, -0.074, -0.261, -2.843,
    0.12, -2.176, 1.362, 1.009, -0.354, 3.655, 1.15, -0.321, -0.074,
    -1.488, -0.593, 1.93, -2.265, -0.249, 1.21, 1.099, 0.531, -1.684,
    -0.456, -1.444, 1.665, -1.835, 1.783, 1.622, -0.592, 0.69, 0.115,
    -0.943, 1.96, 0.031, -0.523, -0.363, 3.047, -2, 1.448, -2.132,
    1.201, -0.844, 0.932, -0.722, 0.024, -3.273, 0.804, -0.051, -0.483,
    0.863, -2.713, -0.429, 0.865, -0.688, 0.498, 0.607, 0.068, 1.671,
    -0.928, 0.925, -0.417, 0.436, 0.554, -2.642, -0.141, 0.486, 0.51,
    -0.914, -1.104, -0.114)
  wdth = 10L

  roll_regress_R_for_loop <- function(X, y, width, grp, downdate){
    u_grp = unique(grp)
    n <- nrow(X)
    p <- ncol(X)
    out <- matrix(NA_real_, n, p)
    sigmas             <- rep(NA_real_, n)
    r.squared          <- rep(NA_real_, n)
    one_step_forecasts <- rep(NA_real_, n)

    start_val <- max(which(u_grp <= width))
    for(g in u_grp[start_val:length(u_grp)]){
      idx <-
        if(downdate)
          which(grp %in% (g - width + 1L):g) else
            which(grp %in% 1:g)
      i <- which(grp == g)
      fit <- lm(y[idx] ~ -1 + X[idx, , drop = FALSE])
      out[i, ] <- sapply(fit$coefficients, rep, times = length(i))

      su <- summary(fit)
      sigmas[i] <- su$sigma

      ss1 <- sum((y[idx] - mean(y[idx]))^2)
      ss2 <- sum(fit$residuals^2)
      r.squared[i] <- 1 - ss2 / ss1

      if(g != max(grp)){
        next_g <- grp[min(which(grp > g))]
        next_g <- which(grp == next_g)
        one_step_forecasts[next_g] <- X[next_g, ] %*% fit$coefficients
      }
    }

    list(coef = out, sigmas = sigmas, r.squared = r.squared,
         one_step_forecasts = one_step_forecasts)
  }

  . <- function(
    do_compute_R_sqs, do_compute_sigmas, do_1_step_forecasts, use_grp,
    do_downdates)
    substitute({
      cpp_out <- roll_cpp(Y = y, X = X, window = wdth,
                          do_compute_R_sqs = do_compute_R_sqs,
                          do_compute_sigmas = do_compute_sigmas,
                          do_1_step_forecasts = do_1_step_forecasts,
                          grp = week, use_grp = TRUE,
                          do_downdates = do_downdates)
      expect_equal(r_out$coef, cpp_out$coefs)
      t1
      t2
      t3
    }, list(
      t1 = if(do_compute_sigmas)
        quote(expect_equal(r_out$sigmas, drop(cpp_out$sigmas))) else
          quote(expect_null(cpp_out$sigmas)),
      t2 = if(do_compute_R_sqs)
        quote(expect_equal(r_out$r.squared, drop(cpp_out$r.squareds))) else
          quote(expect_null(cpp_out$r.squareds)),
      t3 = if(do_1_step_forecasts)
        quote(expect_equal(
          r_out$one_step_forecasts, drop(cpp_out$one_step_forecasts))) else
            quote(expect_null(cpp_out$one_step_forecasts)),
      do_compute_R_sqs = do_compute_R_sqs,
      do_compute_sigmas = do_compute_sigmas,
      do_1_step_forecasts = do_1_step_forecasts, do_downdates = do_downdates))

  vals <- expand.grid(
    do_compute_R_sqs    = c(T, F),
    do_compute_sigmas   = c(T, F),
    do_1_step_forecasts = c(T, F),
    do_downdates        =   T)

  r_out <- roll_regress_R_for_loop(X, y, wdth, grp = week, downdate = TRUE)
  for(i in 1:nrow(vals))
    eval(do.call(., as.list(vals[i, ])))

  r_out <- roll_regress_R_for_loop(X, y, wdth, grp = week, downdate = FALSE)
  vals$do_downdates <- FALSE
  for(i in 1:nrow(vals))
    eval(do.call(., as.list(vals[i, ])))
})

test_that("'.find_grps' gives the same as a simpler R-version", {
  # Idea: write a simpler R function where we can check the result of
  #       View(do.call(cbind, R_func(...)$tmp)) and then test cpp function
  #       against this R function
  R_func <- function(grp, width, min_obs){
    out <- istart <- integer(length(grp))
    for(i in 1:length(grp)){
      idx <- which(grp == grp[i])
      keep <- grp > grp[i] - width & grp <= grp[i]
      out[idx] <- sum(keep)
      istart[idx] <- min(which(keep))
    }

    can_comp <- out >= min_obs

    # we first want to set values after a full window
    can_comp[grp - min(grp) < width - 1L] <- FALSE

    need_restart <- c(FALSE, head(can_comp, -1) != can_comp[-1] & can_comp[-1])

    chunk_grp <- cumsum(need_restart) * can_comp
    i_start     <- tapply(istart, chunk_grp, min)[-1]
    i_grp_start <- tapply(seq_along(grp), chunk_grp, min)[-1]
    i_end       <- tapply(seq_along(grp), chunk_grp, max)[-1]

    tmp <- cbind(nobs = out, istart = istart, can_comp = can_comp,
                 need_restart = need_restart, grp = grp, chunk_grp = chunk_grp)

    list(
      tmp = tmp, i_start = i_start, i_grp_start = i_grp_start, i_end = i_end)
  }

  # simple case that we have all observations
  x <- (1:1000 - 1L) %/% 4L + 1L
  wth <- 5L
  min_obs <- 10L

  expect_equal(
    .find_chunks(x, wth, min_obs),
    list(grp_idx_start = 1L, grp_idx_stop = 1000L, has_value_start = 17L))

  set.seed(36385244)
  for(i in 1:25){
    # drop random entries
    x1 <- x[seq_along(x) %in% sample.int(length(x), 500)]
    t1 <-       R_func(x1, wth, min_obs)
    # View(t1$tmp)
    t2 <- .find_chunks(x1, wth, min_obs)

    expect_length(t2$grp_idx_start, length(t2$has_value_start))
    expect_length(t2$grp_idx_stop , length(t2$has_value_start))
    expect_true(all(t2$grp_idx_start   <= t2$has_value_start))
    expect_true(all(t2$has_value_start <= t2$grp_idx_stop))

    expect_equal(t1$i_start, t2$grp_idx_start, check.attributes = FALSE)
    expect_equal(t1$i_grp_start, t2$has_value_start, check.attributes = FALSE)
    expect_equal(t1$i_end, t2$grp_idx_stop, check.attributes = FALSE)
  }
})

