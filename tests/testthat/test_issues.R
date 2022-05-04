context("testing reported GitHub issues are fixed")

test_that("issue 10 is solved", {
  ehbSample <- structure(list(
    group_ = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 10L, 10L, 10L, 10L, 13L),
    y = c(7, 0, -1, 1, 3, 13, 9, 4, 32, -15, 28, -23, -40, -15, 57, 7, -24, -4, 10, -6, 7, 10, 16, -2, 16, -10, 4, -6, 5, 1, -2, 2, 4, 10, -28, -11, 5, -8, 17, -9, -4),
    x = c(0, -6, 0, -1, -13, 3, 13, 9, -17, 32, -15, 28, -10, 12, -15, 57, -9, -24, -4, 10, -1, 7, 10, 16, -1, 16, -10, 4, 8, 8, 1, -2, -1, 4, 10, -28, 3, 5, -8, 17, 1)), row.names = c(NA, -41L),
    groups = structure(list(
      group_ = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 13L),
      .rows = structure(
        list(1:4, 5:8, 9:12, 13:16, 17:20, 21:24, 25:28, 29:32, 33:36, 37:40, 41L), ptype = integer(0),
             class = c("vctrs_list_of", "vctrs_vctr", "list"))),
      row.names = c(NA, -11L), class = c("tbl_df", "tbl", "data.frame"), .drop = TRUE),
    class = c("grouped_df", "tbl_df", "tbl", "data.frame"))

  fit <- roll_regres(
    y ~ x,
    ehbSample,
    width = 12,
    grp = ehbSample$group_,
    do_downdates = TRUE)

  true_coef <- lm(y ~ x, ehbSample, group_ %in% 2:13)

  expect_true(all(is.na(fit$coefs[seq_len(NROW(ehbSample) - 1L)])))
  expect_equal(fit$coefs[NROW(ehbSample), ], coef(true_coef))
})
