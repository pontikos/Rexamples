source('~nikolas/thor/bin/normalise-functions.R')

x0 <- rnorm(10000, mean=1, sd=1)
x1 <- rnorm(1000, mean=5, sd=3)


x0 <- c(rnorm(1000, mean=1, sd=1), rnorm(2000, mean=4, sd=1))
x1 <- c(rnorm(2000, mean=1, sd=1), rnorm(1000, mean=4, sd=1))

x0 <- c( rnorm(1000, mean=1, sd=1), rnorm(2000, mean=5, sd=1), rnorm(2000, mean=10, sd=1))
x1 <- c( rnorm(1000, mean=1, sd=1), rnorm(2000, mean=4, sd=1), rnorm(2000, mean=4, sd=1))


x0 <- c(rnorm(1000, mean=1, sd=1), rnorm(2000, mean=2, sd=1))





#qqn
#quantile normalisation works well when the shape of the distributions
#is the same and only shifted
plot(density(x0), xlim=range(c(x0,x1)))
lines(density(x1), lty=2) 
lines(density(quantile.normalize(x0, x1)), col='red', lty=2)


#features norm
plot(density(x0, kernel='triangular'), xlim=range(c(x0,x1)))
lines(density(x1, kernel='triangular'), lty=2) 
#lines(density(features.normalize(x0, x1)), col='red', lty=2)


#scale
x0 <- scale(x0)
x1 <- scale(x1)
plot(density(x0), xlim=range(c(x0,x1)))
lines(density(x1), lty=2) 
lines(density(features.normalize(x0, x1)), col='red', lty=2)



s <- 1
k <- 5
x <- x0
d <- exp(-dist(x)**2/s)
g <- graph.adjacency(d, diag=FALSE)
g.l <- graph.laplacian(g, normalized=FALSE)
plot(eigen(g.l)$values[1:(2*k)], main=k)
abline(v=k)










