#' @title Fitting Rolling Linear Models
#' @description Method for fast rolling regression models. I.e., linear models
#' estimated over a moving window of data. The function assumes that
#' @importFrom stats model.frame na.fail
#' @export
roll_regres <- function(
  formula, data, width, contrasts = NULL, do_compute = character()){
  # get model matrix and response
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- quote(na.fail)
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  tt <- terms(mf)
  x <- model.matrix(tt, mf, contrasts = contrasts)
  y <- model.response(mf)

  # run a few checks
  assert_int(width, lower = ncol(x) + 1L, upper = nrow(x))

  # find results
  out <- roll_regres.fit(x = x, y = y, width = width, do_compute = do_compute)

  # add contrasts as an attribuate
  attr(out, "terms") <- attr(x, "contrasts")
  out
}

#' @title Fitter Function for Rolling Linear Models
#' @description Function with a few validation before calling C++ code.
#'
#' @inheritParams roll_regres
#' @param x design matrix of dimension \code{n * p}.
#' @param y numeric vector of observations of length \code{n}.
#'
#' @seealso \code{\link{roll_regres}} for method similar to \code{\link{lm}}.
#'
#' @importFrom checkmate assert_int assert_matrix assert_numeric assert_character
#' @export
roll_regres.fit <- function(x, y, width, do_compute = character()){
  assert_matrix(x, any.missing = FALSE)
  assert_numeric(y, finite = TRUE, any.missing = FALSE, len = nrow(x))
  assert_int(width, lower = ncol(x) + 1L, upper = nrow(x))
  assert_character(do_compute, any.missing = FALSE)
  if(!all(do_compute %in% c("sigmas", "r.squareds")))
    stop(sQuote(do_compute), " contains elements which are not implemented")

  do_compute_sigmas  <- "sigmas"     %in% do_compute
  do_compute_R_sqs   <- "r.squareds" %in% do_compute

  out <- roll_cpp(
    Y = y, X = x, window = width, do_compute_R_sqs = do_compute_R_sqs,
    do_compute_sigmas = do_compute_sigmas)

  # set dimnames
  dimnames(out$coefs) <- dimnames(x)
  if(do_compute_R_sqs){
    out$r.squareds <- drop(out$r.squareds)
    names(out$r.squareds) <- rownames(x)
  }
  if(do_compute_sigmas){
    out$sigmas <- drop(out$sigmas)
    names(out$sigmas) <- rownames(x)
  }

  out
}
