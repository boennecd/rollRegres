#' @title Fitting Rolling Linear Models
#' @description Method for fast rolling regression models. I.e., linear models
#' estimated over a moving window of data. The function assumes that data is
#' orded.
#'
#' @param formula as \code{formula} as in \code{\link{lm}}.
#' @param data an optional \code{data.frame} containing the variables in the model.
#' @param contrasts list passed to \code{\link{model.matrix.default}}.
#' @param width integer with the width of the moving window.
#' @param do_compute character vector with elements \code{"sigmas"},
#' \code{"r.squareds"}, and/or \code{"1_step_forecasts"} for additional output
#' to be computed. See "Details" in \code{\link{roll_regres}}.
#'
#' @details
#' \code{do_compute} can contain \code{"sigmas"} for the estimated standard
#' deviation of the residuals, \code{"r.squareds"} for the \eqn{R^2} of the
#' models, and \code{"1_step_forecasts"} for the out-of-sample forecast for the
#' next periods value.
#'
#' @return
#' List with vector and matrices with the computed output. See the
#' \code{do_compute} argument.
#'
#' @seealso \code{\link{roll_regres.fit}} for method that avoids the call to
#' e.g., \code{\link{model.frame}}.
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
#' @details
#' First, the \code{dqrdc} from routine from LINPACK is used to form the QR
#' decomposition for the first window of data using Householder transformations
#' without pivoting. Then, the LINPACK \code{dchud} and \code{dchdd} routines
#' are used to update and downdate the cholesky decomposition (the R matrix in
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
#' @importFrom checkmate assert_int assert_matrix assert_numeric assert_character
#' @export
roll_regres.fit <- function(x, y, width, do_compute = character()){
  assert_matrix(x, any.missing = FALSE)
  assert_numeric(y, finite = TRUE, any.missing = FALSE, len = nrow(x))
  assert_int(width, lower = ncol(x) + 1L, upper = nrow(x))
  assert_character(do_compute, any.missing = FALSE)
  if(!all(do_compute %in% c("sigmas", "r.squareds", "1_step_forecasts")))
    stop(sQuote(do_compute), " contains elements which are not implemented")

  do_compute_sigmas   <- "sigmas"           %in% do_compute
  do_compute_R_sqs    <- "r.squareds"       %in% do_compute
  do_1_step_forecasts <- "1_step_forecasts" %in% do_compute

  out <- roll_cpp(
    Y = y, X = x, window = width, do_compute_R_sqs = do_compute_R_sqs,
    do_compute_sigmas = do_compute_sigmas,
    do_1_step_forecasts = do_1_step_forecasts)

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
