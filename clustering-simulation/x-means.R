library(mvtnorm)
library(car)
library(mclust)


### Variants of kmeans with different weight update functions.
### By playing around with distance functions and weight updates, I can show how kmeans and GMM are related.
### This toy example also shows what can go wrong with kmeans and GMM.
### For example, bad things happen when clusters shrink to size zero.
### Also the covariance matrix cannot be less than p dimensional otherwise it is no longer invertible.
### This that in p dimensions the cluster must contain at least p+1 elements.


### different inverse distance functions
### Inverse distances because they get larger
### when real distance gets smaller.

euclidean.dist <- function(x,y,...) sqrt(sum((x-y)**2))
#gaussian.dist <- function(x,y,sigma, ...) 1/(2*pi*sqrt(det(sigma)))*exp(-.5*mahalanobis(x,y,sigma))
gaussian.dist <- function(x,y,sigma, ...) dmvnorm(x,mean=y,sigma=sigma)
kmeans.dist <- function(x,y, ...) dmvnorm(x,mean=y,sigma=diag(1,2))
exp.dist <- function(x,y, ...) exp(-euclidean.dist(x,y)**2)


### The weight update depends on the distance function.

### soft weight update, allows for weights in [0; 1]
weights.update.soft <- function(x, theta, dist.fun=gaussian.dist, ...) {
  #dim(z <- t(apply(x, 1, function(y) sapply(1:length(theta$mu), function(i) theta$tau[[i]]*dist.fun(y,theta$mu[[i]],theta$sigma[[i]])))))
  dim(z <- t(apply(x, 1, function(y) sapply(1:length(theta$mu), function(i) dist.fun(y,theta$mu[[i]],theta$sigma[[i]])))))
  #all rows must sum up to 1
  z <- z/rowSums(abs(z))
  return(z)
}

### hard weight update, just selects the max, weight is either 0 or 1
weights.update.hard <- function(x, theta, dist.fun=kmeans.dist, ...) {
  z <- weights.update.soft(x, theta, dist.fun, ...)
  z <- t(apply(z, 1, function(z) {
        i <- which.max(z)
        replace(rep(0,length(z)),i,1)
  }))
  return(z)
}

#cost/likelihood/happiness function
happiness <- function(x, z, theta) {
  sum(sapply(1:length(theta$mu),function(k) mean((z[,k]*x-theta$mu[[k]])**2)))
}

# Here I call the M-step before E-step because I start with a z matrix rather than with
# parameter estimates.  I find it simpler to initialize but maybe it's not a good idea?
# Does it make a difference?
# It would be better if I had the choice of passing in either the weights or the parameters
ME.step <- function(x, z, update.weights, ...) {
  # M step
  #update parameters
  #average weight per component
  tau <- colMeans(z)
  #weighted mean
  mu <- lapply(1:ncol(z), function(i) apply(x,2,weighted.mean,z[,i]))
  #weighted covariance
  sigma <- lapply(1:ncol(z), function(i) cov.wt(x,z[,i])$cov)
  print(theta <- list(tau=tau, mu=mu, sigma=sigma))
  #cost/likelihood/happiness function
  print(h <- happiness(x,z,theta))
  # E step
  #update weights
  z <- update.weights(x, theta, ...)
  return(list(theta=theta, z=z, clusters=apply(z,1,which.max)))
}



# Our first toy dataset will be the "wreath" dataset:
# A dataset consisting of 1000 observations drawn from a
# 14-component normal mixture in which the covariances of the
# components have the same size and shape but different orientation.
K <- 2
data(wreath)
x <- scale(wreath)

# First we will introduce the kmeans algorithm.

##hard weights
#initial random assignment
z <- t(replicate(nrow(x),sample(c(numeric(K-1),1))))
res <- ME.step(x, z, weights.update.hard)
plot(x, col=res$clusters, pch=20, main=0)
for (i in 1:ncol(res$z)) {
  points(t(res$theta$mu[[i]]), col=i, pch='X', cex=2)
  lines(ellipse(res$theta$mu[[i]], shape=res$theta$sigma[[i]], radius=1, col=i), col=i)
}

#10 hard ME steps
for (i in 1:10) {
  plot(x, col=res$clusters, pch=20, main=i)
  for (i in 1:ncol(res$z)) {
    points(t(res$theta$mu[[i]]), col=i, pch='X', cex=2)
    lines(ellipse(res$theta$mu[[i]], shape=res$theta$sigma[[i]], radius=1, col=i), col=i)
  }
  res <- ME.step(x, res$z, weights.update.hard)
}


# Next we will look at the GMM which uses the Gaussian distribution as the distance function.


##soft weights
#initial random assignment
#completely flat no update will happen
#z <- t(replicate(nrow(x),rep(1/K,K)))
#dirichlet
require(gtools)
#z <- rdirichlet(nrow(x),alpha=c(.1,.4,.9))
z <- rdirichlet(nrow(x),alpha=rep(1/K,K))
plot(x, col=apply(z, 1, which.max), pch=20, main=0)
res <- ME.step(x, z, weights.update.soft)
for (i in 1:ncol(res$z)) {
  points(t(res$theta$mu[[i]]), col=i, pch='X', cex=2)
  lines(ellipse(res$theta$mu[[i]], shape=res$theta$sigma[[i]], radius=1, col=i), col=i)
}

#10 soft ME steps
for (i in 1:10) {
  plot(x, col=res$clusters, pch=20, main=i)
  for (i in 1:ncol(res$z)) {
    points(t(res$theta$mu[[i]]), col=i, pch='X', cex=2)
    lines(ellipse(res$theta$mu[[i]], shape=res$theta$sigma[[i]], radius=1, col=i), col=i)
  }
  res <- ME.step(x, res$z, weights.update.soft)
}



