guo.theta <- c(1.2576, 994.3661, 1.1635, 908.9458, 1.1308, 840.1141)

ll <- md_loglike_weibull_series_C1_C2_C3(
    guo_weibull_series_md,
    deltavar = NULL)

ll.grad <- md_score_weibull_series_C1_C2_C3(
    guo_weibull_series_md,
    deltavar = NULL)

theta.hat <- random_restarts(ll, ll.grad, runif(6, 5, 100),
    lr=1, maxit=1000, tol=1e-3, n_restarts=10,
    lower=5, upper=200, proj=function(x) pmax(x, 1e-3))

res <- optim_dynamic(par=rep(10,6), fn=ll, maxit=10000)

res <- sim_anneal(par=rep(100,6), fn=ll,
                control=list(
                    fnscale=-1,
                    maxit=100000L,
                    it_per_temp=100,
                    t_init=100,
                    t_end=1e-20,
                    alpha=.999,
                    REPORT=1,
                    debug=1,
                    sup=function(x) all(x > 0)))


theta.sann <- optim(par=c(10,200,30,400,50,600), fn=ll, method="SANN",
    control=list(
        fnscale=-1, maxit=1000L, REPORT=100,
        temp=100, tmax=20, parscale=par.scale))

 theta.sann <- c(0.1505641, 37.5694586 ,
    0.1310199, 32.6459371 , 0.2852728, 80.0075632)

par.scale <- c(1, 1/1000, 1, 1/1000, 1, 1/1000)
par.scale <- c(1, 1000, 1, 1000, 1, 1000)


test <- c(10,10,10,10,10,10)


parscale(theta.sann, ll)
