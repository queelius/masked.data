

true_theta <- c(1, 1.1, 0.9)
p <- .35
R <- 2000
n <- 100
q <- 0.3
tau <- -(1/sum(true_theta))*log(q)

I0 <- matrix(0, length(true_theta), length(true_theta))
for (i in 1:R)
{
    if (i %% 200 == 0) cat(i / R, "\n")
    data <- generate_data(1, true_theta, tau, p)
    v <- md_score_exp_series_C1_C2_C3(data)(true_theta)
    I0 <- I0 + v %*% t(v)
}
(E_I <- I0 / R)

data <- generate_data(R, true_theta, tau, p)
I1 <- md_fim_exp_series_C1_C2_C3(data)(true_theta)
(E_I1 <- I1 / R)

data1 <- generate_data(R*100, true_theta, tau, p)
theta.hat <- mle_numerical(
    optim(
        par = true_theta,
        fn = md_loglike_exp_series_C1_C2_C3(data1),
        control = list(fnscale = -1, maxit = 10000L)))
I2 <- md_fim_exp_series_C1_C2_C3(data1)(point(theta.hat))
(E_I2 <- I2 / R / 100)

I3 <- md_fim_exp_series_C1_C2_C3(data1)(true_theta)
(E_I3 <- I3 / R / 100)
