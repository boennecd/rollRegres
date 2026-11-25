#' @title Fitting Rolling and Expanding Linear Models
#' @description Method for fast rolling and expanding regression models. I.e.,
#' linear models estimated over a moving window or expanding window of data.
#' The function assumes that data is ordered.
#'
#' @param formula as \code{formula} in \code{\link[stats]{lm}}.
#' @param data an optional \code{data.frame} containing the variables in the model.
#' @param contrasts list passed to \code{\link[stats:model.matrix]{model.matrix.default}}s
#' \code{contrasts.arg} argument.
#' @param width integer with the width of the moving window. Only used if
#' \code{do_downdates == TRUE}.
#' @param do_compute character vector of optional outputs. See Details in
#' \code{\link{roll_regres}} for the supported values.
#' @param grp integer vector to be used if you e.g., want to run the regression
#' over weekly blocks of data. See "Details" in \code{\link{roll_regres}}.
#' @param do_downdates logical which is \code{TRUE} if you want a rolling
#' window regressions. Otherwise, an expanding window is used.
#' @param min_obs positive integer with minimum number of observation that are
#' required in a window. Useful if there are gaps in \code{grp} or unequal
#' number of observations for each \code{grp}.
#'
#' @details
#' \code{do_compute} can include any of the following values.
#' \describe{
#'   \item{\code{"sigmas"}}{Residual standard deviations.}
#'   \item{\code{"r.squareds"}}{\eqn{R^2} values.}
#'   \item{\code{"1_step_forecasts"}}{Out-of-sample forecasts for the next period.}
#'   \item{\code{"vcov"}}{Coefficient covariance matrices.}
#' }
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
#' An object of class \code{roll_regres}. This is a list with a coefficient
#' matrix in \code{coefs} and fields for the optional outputs.
#' \code{coef()} returns the coefficient matrix.
#' When \code{do_compute} includes \code{"vcov"}, \code{vcov()} returns the
#' \code{vcov} field, an array with dimensions \code{p * p * n}.
#' Here, \code{p} is the number of coefficients and \code{n} is the number of
#' input rows. Slice \code{[, , i]} corresponds to coefficient row \code{i}.
#' With \code{se = TRUE}, \code{vcov()} returns an \code{n * p} matrix of
#' standard errors with the same dimnames as the coefficient matrix.
#' \code{vcov()} requires a fit computed with \code{do_compute = "vcov"}.
#' Covariance entries are \code{NA} for rows without a fit and \code{NaN} for
#' fitted windows with zero residual degrees of freedom.
#'
#' @seealso \code{\link{roll_regres.fit}} for method that avoids the call to
#' e.g., \code{\link[stats]{model.frame}}.
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
#' # compute coefficients
#' out <- roll_regres(formula = y ~ X1 + X2, data = df, width = 45L)
#' tail(x = coef(object = out))
#'
#' # compute more output
#' out <- roll_regres(
#'  y ~ X1 + X2, df, width = 45L,
#'  do_compute = c("sigmas", "r.squareds", "1_step_forecasts"))
#' lapply(out, tail)
#'
#' # compute coefficient covariance and standard errors
#' out <- roll_regres(formula = y ~ X1 + X2, data = df, width = 45L,
#'                    do_compute = "vcov")
#' vcov(object = out)[, , n]
#' tail(x = vcov(object = out, se = TRUE))
#'
#' @importFrom stats model.frame na.fail
#' @export
roll_regres <- function(
  formula, data, width, contrasts = NULL, do_compute = character(),
  grp = NULL, do_downdates = TRUE, min_obs = NULL){
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
    x = x, y = y, width = width, do_compute = do_compute, grp = grp,
    do_downdates = do_downdates, min_obs = min_obs)

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
#' The fitter uses the LINPACK routine \code{dqrdc} to compute an unpivoted
#' Householder QR decomposition of the first window. It uses \code{dqrsl}
#' to apply the orthogonal transformation to the response. For subsequent
#' windows, it uses \code{dchud} to add rows and \code{dchdd} to remove rows.
#' The fitter assumes that the design matrix has full column rank in each window.
#' It raises an error if a fitted window has fewer rows than coefficients.
#' With equal counts, the residual standard deviation and coefficient
#' covariance matrix are undefined (\code{NaN}).
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
#' @seealso \code{\link{roll_regres}} for method similar to \code{\link[stats]{lm}}.
#'
#' @importFrom checkmate assert_int assert_matrix assert_numeric
#' @importFrom checkmate assert_character assert_integer
#' @importFrom stats complete.cases
#' @export
roll_regres.fit <- function(
  x, y, width, do_compute = character(), grp = NULL, do_downdates = TRUE,
  min_obs = NULL){
  #####
  # checks
  assert_matrix(x, any.missing = FALSE)
  assert_numeric(y, finite = TRUE, any.missing = FALSE, len = nrow(x))
  if(!is.null(grp))
    assert_int(width, upper = nrow(x)) else
      assert_int(width, lower = max(1L, ncol(x = x)), upper = nrow(x))
  assert_character(do_compute, any.missing = FALSE)
  if(!all(do_compute %in% c("sigmas", "r.squareds", "1_step_forecasts", "vcov")))
    stop(sQuote(do_compute), " contains elements which are not implemented")

  if(is.null(grp)){
    use_grp <- FALSE
    grp <- 1:nrow(x)

  } else {
    assert_integer(grp, null.ok = FALSE, sorted = TRUE, len = nrow(x))
    if(any(diff(grp) < 0))
      stop(sQuote("grp"), " is not increasing. Data needs to be sorted")
    use_grp <- TRUE

  }

  if(is.null(min_obs)){
    min_obs <- ncol(x) * 2L

  } else {
    assert_int(min_obs, na.ok = FALSE, lower = 1L)

  }

  #####
  # build up call
  do_compute_sigmas   <- "sigmas"           %in% do_compute
  do_compute_R_sqs    <- "r.squareds"       %in% do_compute
  do_1_step_forecasts <- "1_step_forecasts" %in% do_compute
  cl <- list(
    quote(.roll_regres.fit), width = quote(width), grp = quote(grp),
    use_grp = use_grp, do_downdates = do_downdates,
    do_compute_sigmas   = do_compute_sigmas,
    do_compute_R_sqs    = do_compute_R_sqs,
    do_1_step_forecasts = do_1_step_forecasts,
    do_compute_vcov     = "vcov" %in% do_compute)
  cl <- as.call(cl)

  #####
  # compute and return
  if(!do_downdates){
    cl[c("y", "x")] <- list(quote(y), quote(x))
    return(.set_names(eval(cl, environment()), dimnames(x)))
  }

  # find chunks
  chunks <- .find_chunks(grp, width, min_obs)
  if(length(chunks$grp_idx_start) == 0L){
    warning("No windows with a sufficient number of observations. Returning ",
            sQuote("NULL"))
    return(NULL)
  }

  if(length(chunks$grp_idx_start) == 1L && chunks$grp_idx_start == 1L &&
     chunks$grp_idx_stop == nrow(x)){
    # pass through data from start to end
    cl[c("y", "x")] <- list(quote(y), quote(x))
    return(.set_names(eval(cl, environment()), dimnames(x)))
  }

  if(do_1_step_forecasts)
    warning(sQuote("1_step_forecasts"), " not implemented with large gaps in",
            " ", sQuote("grp"), ". The results do not contain all values")

  # compute for each chunk
  use_min_obs <- logical(length(chunks$grp_idx_start))
  # we may have to also set the first index to false in the following case.
  # `grp` has values 1 3 4 and width is 2. The first group should be 3 4
  # with start index to and end index 4. The first window though only has
  # length 4 - 3 + 1 = 2 but we have seen a 1 before and 4 - 1 >= the `width`
  if(length(use_min_obs) > 1L)
    use_min_obs[-1] <- TRUE
  use_min_obs[1] <- grp[chunks$has_value_start[1]] - grp[1L] >= width
  out <- mapply(function(grp_idx_start, grp_idx_stop, has_value_start,
                         use_min_obs){
    # compute
    idx <- grp_idx_start:grp_idx_stop
    cl[c("y", "x", "grp", "min_obs", "use_min_obs")] <- list(
      quote(y[idx]), quote(x[idx, , drop = FALSE]), quote(grp[idx]),
      min_obs, use_min_obs)
    o <- eval(cl, environment())

    # check and only keep the rows we need to insert into
    idx_out <- has_value_start:grp_idx_stop
    keep <- complete.cases(o$coefs)
    stopifnot(sum(keep) >= length(idx_out),
              all(keep[seq_len(sum(!keep))] == FALSE))
    keep[seq_len(nrow(o$coefs) - length(idx_out))] <- FALSE
    o$coefs <- o$coefs[keep, , drop = FALSE]
    if(!is.null(o$vcov))
      o$vcov <- o$vcov[, , keep, drop = FALSE]
    other <- !names(x = o) %in% c("coefs", "vcov")
    o[other] <- lapply(o[other], "[", keep)

    # return with the index to insert at
    list(out = o, idx = idx_out)
  }, grp_idx_stop = chunks$grp_idx_stop, grp_idx_start = chunks$grp_idx_start,
  has_value_start = chunks$has_value_start, use_min_obs = use_min_obs,
  SIMPLIFY = FALSE)

  # insert values
  n <- nrow(x)
  res <- lapply(out[[1]]$out, function(x){
    if(is.matrix(x)){
      return(matrix(NA_real_, n, ncol(x)))
    } else if(is.array(x)){
      return(array(data = NA_real_, dim = c(dim(x = x)[1:2], n)))
    } else if(is.vector(x) && is.integer(x)){
      return(rep(NA_integer_, n))
    } else if(is.vector(x) && is.numeric(x)){
      return(rep(NA_real_, n))
    } else if(is.null(x))
      return(NULL)

    stop("Invalid type")
  })

  scalar_outputs <- setdiff(x = names(x = res), y = c("coefs", "vcov"))
  for(z in out){
    res$coefs[z$idx, ] <- z$out$coefs
    if(!is.null(res$vcov))
      res$vcov[, , z$idx] <- z$out$vcov

    for(i in scalar_outputs){
      if(is.null(z$out[[i]]))
        next
      res[[i]][z$idx] <- z$out[[i]]
    }
  }

  .set_names(res, dimnames(x))
}

.roll_regres.fit <- function(
  y, x, width, do_compute_sigmas, do_compute_R_sqs, do_1_step_forecasts, grp,
  use_grp, do_downdates, do_compute_vcov, min_obs = 0L, use_min_obs = FALSE){
  roll_cpp(
    Y = y, X = x, window = width, do_compute_R_sqs = do_compute_R_sqs,
    do_compute_sigmas = do_compute_sigmas,
    do_1_step_forecasts = do_1_step_forecasts, grp = grp, use_grp = use_grp,
    do_downdates = do_downdates, min_obs = min_obs, use_min_obs = use_min_obs,
    do_compute_vcov = do_compute_vcov)
}

.set_names <- function(out, .dimnames){
  dimnames(out$coefs) <- .dimnames
  if(!is.null(out$r.squareds))
    names(out$r.squareds) <- .dimnames[[1]]

  if(!is.null(out$sigmas))
    names(out$sigmas) <- .dimnames[[1]]

  if(!is.null(out$one_step_forecasts))
    names(out$one_step_forecasts) <- .dimnames[[1]]

  if(!is.null(out$vcov))
    dimnames(out$vcov) <- .dimnames[c(2L, 2L, 1L)]

  class(out) <- c("roll_regres", "list")
  out
}

#' @rdname roll_regres
#' @param object a fit from \code{roll_regres()} or \code{roll_regres.fit()}.
#' @param ... additional arguments. These methods ignore them.
#' @importFrom stats coef
#' @export
coef.roll_regres <- function(object, ...) {
  object$coefs
}

#' @rdname roll_regres
#' @param se logical. If \code{TRUE}, return standard errors.
#' @importFrom stats vcov
#' @export
vcov.roll_regres <- function(object, se = FALSE, ...) {
  covariance <- object$vcov
  if(is.null(covariance))
    stop('Fit the model with do_compute = "vcov" to compute covariance.', call. = FALSE)
  if(!se)
    return(covariance)

  standard_errors <- matrix(data = NA_real_, nrow = nrow(x = object$coefs),
                            ncol = ncol(x = object$coefs),
                            dimnames = dimnames(x = object$coefs))
  for(j in seq_len(length.out = ncol(x = standard_errors)))
    standard_errors[, j] <- sqrt(x = covariance[j, j, ])
  standard_errors
}
