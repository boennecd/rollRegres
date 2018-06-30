#' @title Fitting Rolling and Expanding Linear Models
#' @description Method for fast rolling and expanding regression models. I.e.,
#' linear models estimated over a moving window or expanding window of data.
#' The function assumes that data is ordered.
#'
#' @param formula as \code{formula} in \code{\link{lm}}.
#' @param data an optional \code{data.frame} containing the variables in the model.
#' @param contrasts list passed to \code{\link{model.matrix.default}}s
#' \code{contrasts.arg} argument.
#' @param width integer with the width of the moving window. Only used if
#' \code{do_downdates == TRUE}.
#' @param do_compute character vector with elements \code{"sigmas"},
#' \code{"r.squareds"}, and/or \code{"1_step_forecasts"} for additional output
#' to be computed. See "Details" in \code{\link{roll_regres}}.
#' @param grp integer vector to be used if you e.g., want to run the regression
#' over weekly blocks of data. See "Details" in \code{\link{roll_regres}}.
#' @param do_downdates logical which is \code{TRUE} if you want a rolling
#' window regressions. Otherwise, an expanding window is used.
#'
#' @details
#' \code{do_compute} can contain \code{"sigmas"} if you want the estimated
#' standard deviation of the residuals, \code{"r.squareds"} for the \eqn{R^2}
#' of the models, and \code{"1_step_forecasts"} for the out-of-sample forecast
#' for the next periods value.
#'
#' \code{grp} is a sorted integer vector if you want to make "block" updates.
#' E.g., \code{grp} could be an integer vector with the week number. The
#' \code{width} argument is relative to the \code{grp} argument if the
#' \code{grp} argument is not \code{NULL}. The indices of \code{grp} should
#' match with the other data objects.
#'
#' See \code{vignette("Comparisons", package = "rollRegres")} for further
#' examples.
#'
#' @return
#' List with vector and matrices with the computed output. See the
#' \code{do_compute} argument.
#'
#' @seealso \code{\link{roll_regres.fit}} for method that avoids the call to
#' e.g., \code{\link{model.frame}}.
#'
#' @importFrom stats terms model.matrix model.response
#'
#' @examples
#' # simulate data
#' set.seed(29132867)
#' n <- 50
#' p <- 2
#' X <- cbind(1, matrix(rnorm(p * n), ncol = p))
#' y <- drop(X %*% c(1, -1, 1)) + rnorm(n)
#' df <- data.frame(y, X[, -1])
#'
#' # compute coefs
#' out <- roll_regres(y ~ X1 + X2, df, width = 45L)
#' tail(out$coefs)
#'
#' # compute more output
#' out <- roll_regres(
#'  y ~ X1 + X2, df, width = 45L,
#'  do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
#' lapply(out, tail)
#'
#' @importFrom stats model.frame na.fail
#' @export
roll_regres <- function(
  formula, data, width, contrasts = NULL, do_compute = character(),
  grp = NULL, do_downdates = TRUE){
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

  # find results
  out <- roll_regres.fit(
    x = x, y = y, width = width, do_compute = do_compute, grp = grp)

  # add contrasts as an attribuate
  attr(out, "terms") <- attr(x, "contrasts")
  out
}

#' @title Fitter Function for Rolling and Expanding Linear Models
#' @description Function with a few validations before calling C++ code.
#'
#' @inheritParams roll_regres
#' @param x design matrix of dimension \code{n * p}.
#' @param y numeric vector of observations of length \code{n}.
#'
#' @details
#' First, the \code{dqrdc} routine from LINPACK is used to form the QR
#' decomposition for the first window of data using Householder transformations
#' without pivoting. Then, the LINPACK \code{dchud} and \code{dchdd} routines
#' are used to update and downdate the Cholesky decomposition (the R matrix in
#' the QR decomposition).
#'
#' Notice that unlike \code{lm}, there are no checks of the rank of the matrix.
#'
#' @references
#' Golub, G. H., & Van Loan, C. F. (2013). Matrix computations (4rd ed.).
#' JHU Press. See chapter 5 and section 6.5.
#'
#' @return
#' Same as \code{\link{roll_regres}}.
#'
#' @examples
#' # simulate data
#' set.seed(9623556)
#' n <- 50
#' p <- 2
#' X <- cbind(1, matrix(rnorm(p * n), ncol = p))
#' y <- drop(X %*% c(1, -1, 1)) + rnorm(n)
#'
#' # compute coefs
#' out <- roll_regres.fit(x = X, y = y, width = 45L)
#' tail(out$coefs)
#'
#' @seealso \code{\link{roll_regres}} for method similar to \code{\link{lm}}.
#'
#' @importFrom checkmate assert_int assert_matrix assert_numeric
#' assert_character assert_integer
#' @export
roll_regres.fit <- function(
  x, y, width, do_compute = character(), grp = NULL, do_downdates = TRUE){
  #####
  # checks
  assert_matrix(x, any.missing = FALSE)
  assert_numeric(y, finite = TRUE, any.missing = FALSE, len = nrow(x))
  assert_int(width, lower = ncol(x) + 1L, upper = nrow(x))
  assert_character(do_compute, any.missing = FALSE)
  if(!all(do_compute %in% c("sigmas", "r.squareds", "1_step_forecasts")))
    stop(sQuote(do_compute), " contains elements which are not implemented")

  if(is.null(grp)){
    use_grp <- FALSE
    grp <- 1:nrow(x)

  } else {
    assert_integer(grp, null.ok = FALSE, sorted = TRUE, len = nrow(x))
    if(diff(range(grp)) > nrow(x))
      stop(sQuote("grp"), " has too large differences")
    use_grp <- TRUE

  }

  #####
  # compute
  do_compute_sigmas   <- "sigmas"           %in% do_compute
  do_compute_R_sqs    <- "r.squareds"       %in% do_compute
  do_1_step_forecasts <- "1_step_forecasts" %in% do_compute

  out <- roll_cpp(
    Y = y, X = x, window = width, do_compute_R_sqs = do_compute_R_sqs,
    do_compute_sigmas = do_compute_sigmas,
    do_1_step_forecasts = do_1_step_forecasts, grp = grp, use_grp = use_grp,
    do_downdates = do_downdates)

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
  if(do_1_step_forecasts){
    out$one_step_forecasts <- drop(out$one_step_forecasts)
    names(out$one_step_forecasts) <- rownames(x)
  }

  out
}
