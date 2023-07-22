#' Exponential series log-likelihood model
#' 
#' Maximum likelihood estimation functions for
#' Exponential series systems from masked data. Functions
#' include the log-likelihood, score, and FIM functions.
#' 
#' Masked component data approximately satisfies the
#' following conditions:
#' 
#' C1: Pr(K in C) = 1
#' C2: Pr(C=c | K=j, T=t) = Pr(C=c | K=j', T=t)
#'     for any j, j' in c.
#' C3: masking probabilities are independent of theta

#'
#' @author Alex Towell
#' @name Exponential series MLE
#' @keywords exponential, distribution, series, statistics, masked data
#' @seealso \code{\link{md_loglike_exp_series_C1_C2_C3}},
#'          \code{\link{md_score_exp_series_C1_C2_C3}},
#'          \code{\link{md_mle_exp_series_C1_C2_C3}}
#'          \code{\link{md_fim_exp_series_C1_C2_C3}}
NULL

#' Maximum likelihood estimator for exponential series
#' in which each component (or competing risk) is
#' parameterized by a single rate parameter.
#'
#' The sample (at least approximately) satisfy
#' competing risks conditions C1, C2, and C3 and
#' is failure time data may be right censored.
#'
#' @param md right-censored masked data
#' @param theta0 initial value for the MLE
#' @param control list of control parameters
#' @param ... pass additional arguments to `control`
#' @return MLE of type `md_mle_exp_series_C1_C2_C3`
#' 
#' @importFrom algebraic.mle mle_numerical sim_anneal mle
#' @importFrom MASS ginv
#' @export
md_mle_exp_series_C1_C2_C3 <- function(md, theta0, control = list(), ...) {
    defaults <- list(
        keep_obs = FALSE,
        maxit = 1000L,
        fnscale = -1,
        sysvar = "t",
        setvar = "x",
        method = "BFGS",
        debug = 0,
        REPORT = 1)

    control <- modifyList(defaults, control)
    optim_control_args <- c("fnscale", "parscale", "ndeps", "maxit", 
                        "abstol", "reltol", "alpha", "beta", 
                        "gamma", "REPORT", "trace", "warn.1d.Nelder-Mead")
    optim_control <- control[names(control) %in% optim_control_args]

    ll <- md_loglike_exp_series_C1_C2_C3(
        md = md,
        setvar = control$setvar,
        sysvar = control$sysvar)

    ll.grad <- md_score_exp_series_C1_C2_C3(
        md = md,
        setvar = control$setvar,
        sysvar = control$sysvar)

    sol <- optim(
        par = theta0,
        fn = ll,
        gr = ll.grad,
        hessian = FALSE,
        method = control$method,
        control = optim_control,
        ...)

    if (sol$convergence != 0) {
        warning("optimization did not converge")
    }

    control$hessian <- md_fim_exp_series_C1_C2_C3(
        md = md, setvar = control$setvar, sysvar = control$sysvar)(sol$par)

    sol <- mle_numerical(
        sol = sol,
        options = control,
        superclass = c("md_mle_exp_series_system_C1_C2_C3",
                       "md_mle_series_system_C1_C2_C3",
                       "md_mle_series_system"))
    if (control$keep_obs) {
        sol$obs <- md
    }
    sol$nobs <- nrow(md)
    sol
}

#' Generates a log-likelihood function for an exponential series system with
#' respect to rate parameter for masked data with candidate sets that satisfy 
#' conditions C1, C2, and C3.
#'
#' @param md masked data
#' @param options list of options
#' - `set.var` prefix of Boolean matrix encoding of candidate sets, defaults
#'   to `x`, e.g., `x1,...,xm`.
#' - `sys.var` system lifetime (optionally right-censored) column name,
#'   defaults to `t`
#' - `delta.var` right-censoring indicator column name. If `NULL`, then
#'   no right-censoring is assumed. If a system lifetime is
#'   right-censored (*not* observed), the right-censoring
#'   indicator is `TRUE`, otherwise it is `FALSE`.
#' @return log-likelihood function
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_loglike_exp_series_C1_C2_C3 <- function(md, options = list(), ...) {

    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    defaults <- list(
        set.var = "x",
        sys.var = "t",
        delta.var = NULL)

    options <- modifyList(defaults, options)
    options <- modifyList(options, list(...))
    stopifnot(
        options$sys.var %in% colnames(md),
        is.null(options$delta.var) || options$delta.var %in% colnames(md))

    sum.t <- sum(md[[options$sys.var]])
    if (!is.null(options$delta.var) && options$delta.var %in% colnames(md)) {
        # only keep the observations that were not right-censored in `md`
        md <- md %>% filter(!.[[options$delta.var]])
    }
    md$C <- md_decode_matrix(md, options$set.var)
    md <- md %>% group_by(.data$C) %>% count()
    m <- ncol(md$C)
    if (m == 0) {
        stop("no candidate sets wih prefix '", options$set.var, "' found")
    }


    # TODO: if m = ncol(md$C) is very large, then this will be slow
    #       having 2^m rows in md. in this case, if 2^m > nrows(md),
    #       then it would be faster to not do any grouping and just
    #       iterate through the original rows of `md`.
    function(theta) {
        if (length(theta) != m) stop("length(theta) != m")
        if (any(theta <= 0)) return(NA)
        f <- -sum.t * sum(theta)
        for (i in seq_len(n)) {
            f <- f + md$n[i] * log(sum(theta[md$C[i, ]]))
        }
        f
    }
}


#' Generates a score function for an exponential series system (or
#' competing risks) with respect to parameter `theta` for masked component
#' failure (or masked competing risks) with candidate sets (risks) that
#' approximately satisfy conditions C1, C2, and C3.
#'
#' @param md right-censored failure-time data with masked competing risks
#' @param options list of options
#' - `set.var` prefix of Boolean matrix encoding of candidate sets, defaults
#'   to `x`, e.g., `x1,...,xm`.
#' - `sys.var` system lifetime (optionally right-censored) column name,
#'   defaults to `t`
#' - `delta.var` right-censoring indicator column name. If `NULL`, then
#'   no right-censoring is assumed. If a system lifetime is
#'   right-censored (*not* observed), the right-censoring
#'   indicator is `TRUE`, otherwise it is `FALSE`.
#' @param ... pass additional arguments to `options`
#' @return score function of type `R^m -> R`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_score_exp_series_C1_C2_C3 <- function(md, options = list(), ...) {
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    defaults <- list(
        set.var = "x",
        sys.var = "t",
        delta.var = NULL)

    options <- modifyList(defaults, options)
    options <- modifyList(options, list(...))
    stopifnot(
        options$sys.var %in% colnames(md),
        is.null(options$delta.var) || options$delta.var %in% colnames(md))

    sum.t <- sum(md[ , options$sys.var])
    if (!is.null(options$delta.var) && options$delta.var %in% colnames(md)) {
        # only keep the observations that were not right-censored in `md`
        md <- md %>% filter(!.[[options$delta.var]])
    }
    md$C <- md_decode_matrix(md, options$set.var)
    md <- md %>% group_by(.data$C) %>% count()
    m <- ncol(md$C)
    if (m == 0) {
        stop("no candidate sets wih prefix '", options$set.var, "' found")
    }

    function(theta) {
        if (length(theta) != m) stop("length(theta) != m")
        if (any(theta <= 0)) return(NA)
        v <- rep(-sum.t, m)
        for (j in seq_len(m)) {
            for (i in seq_len(n)) {
                if (md$C[i, j]) {
                    v[j] <- v[j] + md$n[i] / sum(theta[md$C[i, ]])
                }
            }
        }
        v
    }
}

#' Generates the observed information matrix (FIM) for an exponential series
#' system with respect to parameter `theta` for masked data with candidate
#' sets that approximately satisfy conditions C1, C2, and C3.
#'
#' @param md masked data with candidate sets that meet the
#'           regular candidate model
#' @return observed information matrix of type `R^m -> R^(m x m)`
#' @importFrom dplyr %>% group_by count filter
#' @importFrom md.tools md_decode_matrix
#' @importFrom rlang .data
#' @export
md_fim_exp_series_C1_C2_C3 <- function(md, options = list(), ...) {
    n <- nrow(md)
    if (n == 0) {
        stop("md is empty")
    }

    defaults <- list(
        set.var = "x",
        sys.var = "t",
        delta.var = NULL)

    options <- modifyList(defaults, options)
    options <- modifyList(options, list(...))
    stopifnot(
        options$sys.var %in% colnames(md),
        is.null(options$delta.var) || options$delta.var %in% colnames(md))

    if (!is.null(options$delta.var) && options$delta.var %in% colnames(md))
        md <- md %>% filter(.data$delta == FALSE)
    md$C <- md_decode_matrix(md, options$set.var)
    md <- md %>% group_by(.data$C) %>% count()
    m <- ncol(md$C)
    if (m == 0) {
        stop("no candidate sets wih prefix '", options$set.var, "' found")
    }

    function(theta) {
        if (length(theta) != m) stop("length(theta) != m")
        if (any(theta <= 0)) return(NA)
        I <- matrix(rep(0, m * m), nrow = m)
        for (j in 1:m) {
            for (k in 1:m) {
                for (i in 1:n) {
                    if (md$C[i, j] && md$C[i, k]) {
                        I[j, k] <- I[j, k] + md$n[i] / sum(theta[md$C[i, ]])^2
                    }
                }
            }
        }
        I
    }
}




###################### SIMULATIONS ######################

#' Generate exponential series data with right-censoring and masked
#' component failure satisfying conditions C1, C2, and C3.
#' 
#' @param n number of observations
#' @param rates vector of component failure rates
#' @param tau numeric, right-censoring time, defaults to rep(Inf,n).
#'            If tau[i] = Inf, then the i-th system is not
#'            right-censored.
#' @param p numeric, `p[i]` is the probability that each of the
#'          non-failed components in the i-th series system
#'          appears in the candidate set. Defaults to rep(0,n).
#'          If p[i] = 0, then only the failed component is in the
#'          candidate set (no masking). If p[i] = 1, then all
#'          all components are in the candidate set, in which case
#'          the candidate set conveys no information about the
#'          (latent) components (the lifetime of the system still
#'          conveys information about the components, but if that's
#'          all the sample has, then the MLE is not
#'          identifiable).
#' @export
md_exp_series_system_bernoulli_cand_C1_C2_C3 <- function(
    n,
    rates,
    p = 0,
    tau = Inf) {
    m <- length(rates)
    comp_times <- matrix(nrow = n, ncol = m)
    for (j in 1:m)
        comp_times[,j] <- rexp(n,rates[j])

    data <- md_encode_matrix(comp_times, "t") %>%
        md_series_lifetime_right_censoring(tau)

    if (is.function(p)) {
        data <- md_bernoulli_cand_C1_C2_C3(data, p)
    } else {
        data <- md_bernoulli_cand_C1_C2_C3(data, rep(p, n))
    }
    data %>% md_cand_sampler()
}
