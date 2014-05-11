#int *n, Sample size
#int *m, Number of components
#double *mu, m-vector of current mean estimates
#double *x, data:  vector of length n
#double *h, bandwidth
#double *z, nn*mm vector of normalized posteriors (or indicators in stochastic case), normalized by "column"
#double *f  KDE evaluated at n*m matrix of points, namely, x_i - mu_j for 1<=i<=n and 1<=j<=m

## Implement symmetric kernel density-estimation step for location 
## mixture model, equation (20) in Benaglia et al (2008)
KDEsymloc <-  function(n, m, mu, x, h, z) {
  f <- matrix(nrow=n,ncol=m)
  ## Must loop n^2*m^2 times; each f entry requires n*m calculations
  for(a in 1:n) {
    for(b in 1:m) {
      sum = 0
      u1 = (x[a]-mu[b])
      for(i in 1:n) {
        for(j in 1:m) {
          tmp1=(x[a]-mu[b])-(x[i]-mu[j])
          tmp2=(mu[b]-x[a])-(x[i]-mu[j])
          # Use normal kernel
          sum <- sum + z[i,j] * mean(c(exp(-.5*tmp1**2/h**2),exp(-.5*tmp2**2/h**2)))
        }
      }
      #(a,b) entry of f matrix
      f[a,b] <- sum /(h*n*sqrt(2*pi))
    }
  }
  return(list(f=f))
}

spEMsymloc <- function (x, mu0, bw = bw.nrd0(x), h = bw, eps = 1e-08, maxiter = 100, stochastic = FALSE, verbose = TRUE, use.C=FALSE) {
    bw <- h
    n <- length(x)
    if (length(mu0) > 1) m <- length(mu0) else m <- mu0
    z.hat <- matrix(0, nrow = n, ncol = m)
    fkernel <- matrix(0, nrow = n, ncol = m)
    tt0 <- proc.time()
    lambda <- rep(1/m, m)
    kmeans <- kmeans(x, mu0)
    for (j in 1:m) {
        z.hat[kmeans$cluster == j, j] <- 1
    }
    iter <- 0
    if (stochastic) {
        sumpost <- matrix(0, n, m)
    }
    finished <- FALSE
    lambda <- mu <- matrix(0, maxiter, m)
    while (!finished) {
        iter <- iter + 1
        t0 <- proc.time()
        print(lambda[iter, ] <- colMeans(z.hat))
        mu[iter, ] <- apply(sweep(z.hat, 1, x, "*"), 2, mean)/lambda[iter,]
        cat('z.hat\n')
        print(dim(z.hat))
        if (stochastic) {
            z <- t(apply(z.hat, 1, function(prob) rmultinom(1, 1, prob)))
        if (use.C)
            ans <- .C("KDEsymloc", n = as.integer(n), m = as.integer(m), mu = as.double(mu[iter, ]), x = as.double(x), bw = as.double(bw), z = as.double(z), f = double(n * m), PACKAGE = "mixtools")
        else
            ans <- KDEsymloc(n, m, mu[iter,], x, bw, z)
        }
        else {
          if (use.C)
            ans <- .C("KDEsymloc", n = as.integer(n), m = as.integer(m), mu = as.double(mu[iter, ]), x = as.double(x), bw = as.double(bw), z = as.double(z.hat), f = double(n * m), PACKAGE = "mixtools")
          else
            ans <- KDEsymloc(n, m, mu[iter,], x, bw, z.hat)
        }
        cat('dim f', dim(fkernel <- matrix(ans$f, ncol = m, nrow=n)), '\n')
        lambda.f <- sweep(fkernel, 2, lambda[iter, ], "*")
        z.hat <- lambda.f/rowSums(lambda.f)
        finished <- iter >= maxiter
        if (stochastic) {
            sumpost <- sumpost + z.hat
        }
        else if (iter > 1) {
            change <- c(lambda[iter, ] - lambda[iter - 1, ], mu[iter, ] - mu[iter - 1, ])
            finished <- finished | (max(abs(change)) < eps)
        }
        if (verbose) {
            t1 <- proc.time()
            cat("iteration ", iter, "  lambda ", round(lambda[iter, ], 4), "  mu ", round(mu[iter, ], 4))
            cat(" time", (t1 - t0)[3], "\n")
        }
    }
    if (verbose) {
        tt1 <- proc.time()
        cat("lambda ", round(lambda[iter, ], 4))
        cat(", total time", (tt1 - tt0)[3], "s\n")
    }
    if (stochastic) {
        return(structure(list(data = x, posteriors = sumpost/iter, lambda = lambda, bandwidth = bw, lambdahat = colMeans(lambda), mu = mu, muhat = colMeans(mu), symmetric = TRUE), class = "npEM"))
    }
    else {
        return(structure(list(data = x, posteriors = z.hat, lambda = lambda[1:iter, ], bandwidth = bw, lambdahat = lambda[iter, ], mu = mu[1:iter, ], muhat = mu[iter, ], symmetric = TRUE), class = "npEM"))
    }
}
