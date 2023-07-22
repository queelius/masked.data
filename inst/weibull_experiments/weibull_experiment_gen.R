library(utils)
library(tidyverse)
library(algebraic.mle)
library(md.tools)
library(stats)

# weibull series system
# parameterization: theta = (shape_1, scale_1, ..., shape_m, scale_m)
# for each component, the failure time is distributed as
#     T_i j ~ Weibull(shape_j, scale_j)
# for i = 1, ..., n and j = 1, ..., m.
weibull_generate_data <- function(n, theta, tau, p) {
    k <- length(theta)
    m <- k / 2
    shapes <- theta[seq(1, k, 2)]
    scales <- theta[seq(2, k, 2)]
    comp_times <- matrix(nrow = n, ncol = m)
    for (j in 1:m)
        comp_times[,j] <- rweibull(n, shapes[j], scales[j])
    comp_times <- md_encode_matrix(comp_times, "t")

    comp_times %>%
        md_series_lifetime_right_censoring(tau) %>%
        md_bernoulli_cand_C1_C2_C3(p) %>%
        md_cand_sampler()
}

weibull_custom_solver <- function(data, theta,
    extra_info = NULL, annealing = TRUE) {
    ll <- md_loglike_weibull_series_C1_C2_C3(data)
    # ll.grad <- md_score_weibull_series_C1_C2_C3(data)
    # fish <- md_fim_weibull_series_C1_C2_C3(data)
    # ll.hess <- function(x) -fish(x)
    theta.hat <- NULL
    start <- NULL

    tryCatch({
        start <- list(par = theta)
        if (annealing) {
            start <- sim_anneal(par = theta, fn = ll, control =
                list(fnscale = -1, maxit = 1000000L, trace = FALSE,
                    t_init = 1, 1e-4, alpha = 0.99,
                    it_per_temp = 100L,
                    sup = function(x) all(x > 0)))
        }
        res_optim <- optim(
            par = start$par,
            fn = ll,
#            gr = ll.grad,
#            method = "L-BFGS-B",
#            lower = 1e-30,
            hessian = TRUE,
            control = list(fnscale = -1, maxit = 1000L))
        theta.hat <- mle_numerical(res_optim)
        summary(theta.hat)
        print(res_optim)

            #options = list(hessian = ll.hess(res_optim$par)))
    }, error = function(e) {
        cat("Sample size", nrow(data), " | Anneal:",
            start$par, " | ", e$message)
        if (!is.null(extra_info)) {
            cat(" | ", extra_info)
        }
        cat("\n")
    })
    theta.hat
}

weibull_experiment_gen <- function(
    R = 999,
    bernoulli_probs,
    taus,
    sample_sizes,
    theta,
    seed,
    use_aneal_start,
    append,
    csv_filename) {

    m <- length(theta)

    if (!append) {
        cnames <- c("R", "p", "tau", "N", paste0("bias",1:m), "mse",
            paste0("se",1:m), paste0("se_asym",1:m), "mse_asym", "mse_asym_hat",
            paste0("coverage",1:m), "sigma.frob")

        # Write column names first
        write.table(t(cnames), file = csv_filename,
            sep = ",", col.names = FALSE,
            row.names = FALSE, append = FALSE, quote = FALSE)
    }

    for (p in bernoulli_probs) {
        cat("Bernoulli probability", p, "\n")

        for (N in sample_sizes) {
            cat("Sample size", N, "\n")

            for (tau in taus) {
                cat("Tau", tau, "\n")
                mles <- matrix(nrow = R, ncol = m)
                CI_lwr <- matrix(nrow = R, ncol = m)
                CI_upr <- matrix(nrow = R, ncol = m)

                j <- 1L
                repeat {
                    data <- weibull_generate_data(N, theta, tau, p)
                    theta.hat <- weibull_custom_solver(
                        data = data,
                        theta = theta,
                        extra_info = paste0("Replicate(", j, ")"),
                        annealing = use_aneal_start)
                    if (is.null(theta.hat)) {
                        next
                    }
                    if (j %% 10 == 0) {
                        cat("Sample size", N, " | Replicate ", j,
                            " | MLE ", point(theta.hat), "\n")
                    }
                    mles[j, ] <- point(theta.hat)

                    CI <- confint(theta.hat)
                    if (any(is.nan(CI)))
                    {
                        print("NaN in CI")
                        print(summary(theta.hat))
                    }
                    CI_lwr[j, ] <- CI[,1]
                    CI_upr[j, ] <- CI[,2]

                    j <- j + 1L
                    if (j > R) {
                        break
                    }
                }

                # compute asymptotics
                theta.mle <- NULL
                while (is.null(theta.mle)) {
                    data <- weibull_generate_data(N, theta, tau, p)
                    theta.mle <- weibull_custom_solver(
                        data = data,
                        theta = theta,
                        extra_info = "asymptotics",
                        annealing = use_aneal_start)
                }
                SE.asym <- se(theta.mle)
                MSE.asym <- mse(theta.hat, theta)
                MSE.asym.hat <- mse(theta.mle)

                # Calculate bias, MSE, and SE for each parameter
                bias <- colMeans(mles) - theta
                MSE <- ((theta - mles)^2)
                sigma <- cov(mles) * ((N - 1) / N)
                SE <- sqrt(diag(sigma))
                fim.theta <- fim(theta.mle)
                sigma.fim <- MASS::ginv(fim.theta)

                frob <- function(x) sqrt(sum(x^2))
                sigma.frob <- frob(sigma - sigma.fim)
                
                # Compute coverage probabilities
                coverage <- colMeans((CI_lwr <= theta) & (theta <= CI_upr))

                datum <- c(R, p, tau, N, bias, MSE, SE, SE.asym,
                    MSE.asym, MSE.asym.hat, coverage, sigma.frob)

                write.table(t(datum), file = csv_filename, sep = ",",
                    col.names = FALSE, row.names = FALSE, append = TRUE)
            }
        }
    }
}

