

grad_clip <- function(grad, clip) {
    function(x) {
        g <- grad(x)
        if (norm(g) > clip) {
            g <- g / norm(g) * clip
        }
        return(g)
    }
}


grad_ascent <- function(f, grad, x0, lr, maxit = 1000L, tol = 1e-20, clip = 1) {
    x <- x0

    for (i in 1:maxit) {
        g <- grad(x)
        # apply grad_clip
        x_new <- x + lr * g / norm(g) * clip
        if (all(abs(g) < tol)) {
            break
        }
        x <- x_new
    }
    x
}


random_restarts <- function(
    f, grad, x0,
    lr = 1, maxit = 100L,
    parscale = rep(1, length(x0)),
    tol = 1e-20, clip = 1, proj = function(x) x,
    n_restarts = 100, lower = 0, upper = 200) {

    m <- length(x0)
    best_f <- f(x0)
    best_x <- x0
    for (i in 1:n_restarts) {
        tryCatch ({
            #res <- proj(grad_ascent(f, grad, x0, lr, maxit, tol, clip))
            res <- optim(par=x0, fn=f, control =
                list(parscale=parscale,
                     fnscale=-1))
            #res <- sim_anneal(par=x0, fn=f,
            #    control=list(
            #        fnscale=-1,
            #        maxit=maxit,
            #        it_per_temp=100,
            #        t_init=100,
            #        t_end=1e-6,
            #        alpha=.95,
                    #REPORT=100,
                    #debug = 1,
            #        sup=function(x) all(x > 0),
            #    )
            #)
            x <- res$par
            #fx <- f(x)
            fx <- res$value
            if (fx > best_f) {
                best_f <- fx
                best_x <- x
                cat(sprintf("New best: %s\n", best_f))
            }
        }, error = function(e) {
            cat(sprintf("Error: %s\n", e$message))
        })

        x0 <- runif(m, lower, upper)
    }
    best_x
}


parscale <- function(x, f, h = 1e-1) {
    pd <- function(xx, j) {
        xx[j] <- xx[j] + h
        xx
    }

    scale <- rep(1, length(x))
    for (i in 1:length(x)) {
        scale[i] <- max(min(100, abs(f(x) - f(pd(x, i)))), 1)
    }
    scale
}



# use stats::optim, maxit=10, updating the parscale from the gradient,
# going for another maxit=10, rince and repeat until total iterations is
# 1000

optim_dynamic <- function(par, fn, maxit = 1000)
{
    it <- 0
    res <- NULL
    repeat
    {
        scaling_factors <- parscale(par, fn)
        tryCatch({
            res <- optim(par = par, fn = fn,
                         control = list(maxit = 100,
                                        fnscale = -1,
                                        parscale = 1/scaling_factors))
        }, error = function(e) {
            cat(sprintf("Error: %s\n", e$message))
        })

        par <- res$par
        it <- it + 10
        cat("Iteration: ", it, ", par: ", res$par, ", parscale: ",
            scaling_factors, "\n")
        if (it >= maxit) break
    }

    res
}
