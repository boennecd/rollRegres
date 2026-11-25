test_that(desc = "covariance and standard errors match lm in each window", code = {
  i <- seq_len(length.out = 24L)
  x <- 2 + sin(x = i) + i / 10
  data <- data.frame(x = x, y = 1 + 2 * x + cos(x = 2 * i),
                     row.names = paste0("row", i))
  width <- 8L

  for (rolling in c(TRUE, FALSE)) {
    fit <- roll_regres(formula = y ~ x, data = data, width = width,
                       do_compute = "vcov", do_downdates = rolling)
    coefficients <- coef(object = fit)
    covariance <- vcov(object = fit)
    standard_errors <- vcov(object = fit, se = TRUE)

    expect_identical(object = coefficients, expected = fit$coefs)
    expect_null(object = fit$sigmas)
    expect_identical(object = dim(x = covariance), expected = c(2L, 2L, 24L))
    expect_identical(object = dimnames(x = covariance),
                     expected = dimnames(x = coefficients)[c(2L, 2L, 1L)])
    expect_identical(object = dimnames(x = standard_errors),
                     expected = dimnames(x = coefficients))
    expect_true(object = all(is.na(x = covariance[, , seq_len(length.out = width - 1L)])))
    expect_true(object = all(is.na(x = standard_errors[seq_len(length.out = width - 1L), ])))

    for (end in width:length(x = i)) {
      start <- if (rolling) end - width + 1L else 1L
      reference <- stats::lm(formula = y ~ x, data = data[start:end, ])
      expect_equal(object = covariance[, , end], expected = stats::vcov(object = reference),
                   tolerance = 1e-12)
      expect_equal(object = standard_errors[end, ],
                   expected = summary(object = reference)$coefficients[, "Std. Error"],
                   tolerance = 1e-12)
      expect_equal(object = coefficients[end, ], expected = stats::coef(object = reference),
                   tolerance = 1e-12)
    }
  }
})

test_that(desc = "grouped covariance follows the coefficient rows across gaps", code = {
  group <- rep(x = c(1L, 2L, 3L, 5L, 6L, 7L, 8L), each = 3L)
  i <- seq_along(along.with = group)
  x <- 2 + sin(x = i)
  data <- data.frame(x = x, y = 1 + 2 * x + cos(x = 2 * i))
  fit <- roll_regres.fit(x = stats::model.matrix(object = ~ x, data = data),
                        y = data$y, width = 2L, grp = group,
                        min_obs = 6L, do_compute = c("sigmas", "vcov"))
  covariance <- vcov(object = fit)
  for (g in unique(x = group)) {
    window <- data[group %in% c(g - 1L, g), ]
    rows <- which(x = group == g)
    if (nrow(x = window) < 6L) {
      expect_true(object = all(is.na(x = covariance[, , rows])))
      next
    }
    reference <- stats::lm(formula = y ~ x, data = window)
    for (row in rows)
      expect_equal(object = covariance[, , row], expected = stats::vcov(object = reference),
                   tolerance = 1e-12)
    expect_equal(object = unname(obj = fit$sigmas[rows]),
                 expected = rep(x = stats::sigma(object = reference), times = length(x = rows)),
                 tolerance = 1e-12)
  }
})

test_that(desc = "covariance is optional and retains singleton dimensions", code = {
  data <- data.frame(y = c(1, 2, 4, 5))
  fit <- roll_regres(formula = y ~ 1, data = data, width = 3L)
  expect_equal(object = names(x = fit),
               expected = c("coefs", "sigmas", "r.squareds", "one_step_forecasts"))
  expect_identical(object = coef(object = fit), expected = fit$coefs)
  expect_error(object = vcov(object = fit), regexp = 'do_compute = "vcov"')

  fit <- roll_regres(formula = y ~ 1, data = data, width = 3L, do_compute = "vcov")
  expect_identical(object = dim(x = vcov(object = fit)), expected = c(1L, 1L, 4L))
  expect_identical(object = dim(x = vcov(object = fit, se = TRUE)), expected = c(4L, 1L))
  reference <- stats::lm(formula = y ~ 1, data = data[2:4, , drop = FALSE])
  expect_equal(object = vcov(object = fit)[, , 4L, drop = FALSE][1L],
               expected = stats::vcov(object = reference)[1L], tolerance = 1e-12)
})

test_that(desc = "windows without residual degrees of freedom have undefined covariance", code = {
  data <- data.frame(x = c(1, 3, 2), y = c(2, 4, 3))
  fit <- suppressWarnings(expr = roll_regres(formula = y ~ x, data = data, width = 2L,
                                            min_obs = 2L, do_compute = "vcov"))
  expect_true(object = all(is.nan(x = vcov(object = fit)[, , 2:3])))
  expect_true(object = all(is.nan(x = vcov(object = fit, se = TRUE)[2:3, ])))
})
