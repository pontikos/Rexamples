
d <- data.frame(x=rnorm(100), y=rnorm(100))
d2 <- data.frame(x=rnorm(mean=10,100), y=rnorm(mean=10,100))
d3 <- data.frame(x=rnorm(mean=5,100), y=rnorm(mean=5,100))
d4 <- data.frame(x=rnorm(mean=10,100), y=rnorm(mean=0,100))
d5 <- data.frame(x=rnorm(mean=0,100), y=rnorm(mean=10,100))
d6 <- data.frame(x=rnorm(mean=5,100), y=rnorm(mean=6,100))
d <- rbind(d,d2,d3, d4, d5, d6)


#d <- data.frame(x=rnorm(100))
#d2 <- data.frame(x=rnorm(mean=10,100))
#d <- rbind(d,d2)


K <- 5

ass <- sample(1:K, dim(d)[1], replace=TRUE)
penalty <- c(  sum( do.call('rbind',by(d, ass, function(x) as.list(sum(dist(x))))) ) )
first <- TRUE
#kmeans algorithm
while (first || (penalty[n] != penalty[n-1])) {
    X <- do.call('rbind',by(d, ass, colMeans))
    ass <- apply(sapply(1:K, function(i) sqrt(rowSums((d-X[i,])**2))), 1, which.min)
    penalty <- c( penalty, sum( do.call('rbind',by(d, ass, function(x) as.list(sum(dist(x))))) ) )
plot(d, col=ass)
points(X, col=as.numeric(rownames(X)), pch='X')
print(penalty)
n <- length(penalty)
print(n)
if (n == 500) break
first <- FALSE
}

kmeans.ass <- ass
kmeans.penalty <- penalty


#COV <- function(x) as.list(mean(apply(x, 1, function(x) do.call('*', as.list(x))))-do.call('*', as.list(colMeans(x))))


colMedoids <- function(x) return(as.numeric(x[which.min(rowSums(as.matrix(dist(x)))),]))
ass <- sample(1:K, dim(d)[1], replace=TRUE)
penalty <- c()
#kmedoids algorithm
penalty <- c(  sum( do.call('rbind',by(d, ass, function(x) as.list(sum(dist(x))))) ) )
first <- TRUE
#kmeans algorithm
while (first || (penalty[n] != penalty[n-1])) {
    X <- do.call('rbind',by(d, ass, colMedoids))
    ass <- apply(sapply(1:K, function(i) (rowSums((d-X[i,])**2))), 1, which.min)
    penalty <- c( penalty, sum( do.call('rbind',by(d, ass, function(x) as.list(sum(dist(x))))) ) )
plot(d, col=ass)
points(X, col=as.numeric(rownames(X)), pch='X')
print(penalty)
n <- length(penalty)
print(n)
if (n == 500) break
first <- FALSE
}


kmedoids.ass <- ass
kmedoids.penalty <- penalty


plot(kmeans.penalty, type='l')
lines(kmedoids.penalty, col='red')


plot(d, col=kmedoids.ass)
plot(d, col=kmeans.ass)

d2 <- d2[,c('x','y')]
h <- hclust(dist(d,method='manhattan'),method='ward')
plot( 1:(2*K), sapply(1:(2*K), function(k) sum( do.call('rbind',by(d, cutree(h,k=k), function(x) as.list(sum(dist(x,method='manhattan'))))) )), type='l' )
abline(v=6)
plot(d, col=cutree(h,k=2), pch=20)

library(ggplot2)
m <- ggplot(d, aes(x = x, y = y))
m <- m + geom_point() 
m + geom_density2d()

d2 <- data.frame(y@exprs[,c(1,9)])
d2$x <- d2[,1]
d2$y <- d2[,2]
m <- ggplot(d2, aes(x = x, y = y))
m + geom_density2d()

h <- hclust(dist(d2,method='manhattan'),method='ward')
plot( 1:(2*K), sapply(1:(2*K), function(k) sum( do.call('rbind',by(d2, cutree(h,k=k), function(x) as.list(sum(dist(x,method='manhattan'))))) )), type='l' )
abline(v=6)
plot(d2, col=cutree(h,k=2), pch=20)


### DENSITY based clustering

library(fpc)
dbscan(d,.01,showplot=TRUE)

dist.obj <- dist(d,method='manhattan')
hclust.obj <- hclust(dist.obj,method='ward')
css.obj <- css.hclust(dist.obj,hclust.obj)
elbow.obj <- elbow.batch(css.obj)
print(elbow.obj)
plot(d, col=cutree(h,k=elbow.obj$k), pch=20)



library(feature)
f <- featureSignif(d)
plot(f)
#this is the grid
abline(h=f$fhat$x.grid[[1]], v=f$fhat$x.grid[[2]], lwd=.25)
#plotting data in grid format
image(f$grad)
image(f$curv)
image(f$fhat$est)
image(f$fhat$est>quantile(f$fhat$est,.95))


#fast way to obtain point density instead of grid density
library(KernSmooth)
# compute fast kernel density estimate
b<-bkde2D(as.matrix(d),.1)
# this returns you a grid
# we will use a fast nearest neighbour method
# to find the closest point in the grid
grid <- expand.grid(b$x1, b$x2)
library(RANN)
nn <- nn2(grid,d,k=1)
d.dens <- cbind(d, dens=as.numeric(b$fhat)[nn$nn.idx])
# did this work?
# plot all points
plot(d.dens[,1:2], pch=20)
# plot 95% quantile of high-density points in red
points(d.dens[d.dens[,3]>quantile(d.dens[,3],.95),1:2], pch=20, col='red')



# clustering metrics
# all of these decrease as we increase the number of clusters
# in the limiting case each point is assigned to one cluster 
library(mvtnorm)
set.seed(1234)
n <- 2000
l <- list(
          list(tau=.2, mu=c(0,0), sigma=diag(2)),
          list(tau=.2, mu=c(1,1), sigma=diag(2)),
          list(tau=.2, mu=c(0,2), sigma=diag(2)),
          list(tau=.2, mu=c(2,0), sigma=diag(2)),
          list(tau=.2, mu=c(2,2), sigma=diag(2))
          )
d <- as.matrix( cbind(rmvnorm(n*l[[1]]$tau, mean=l[[1]]$mu, sigma=l[[1]]$sigma), 1) )
for (i in 2:length(l)) d <- rbind(d, cbind(rmvnorm(n*l[[i]]$tau, mean=l[[i]]$mu, sigma=l[[i]]$sigma), i))
d <- data.frame(x=d[,1], y=d[,2], label=factor(d[,3]))

#sum of intra cluster distance
total.intra.dist <- function(d) sum( do.call('rbind',by(d[,1:2], d[,3], function(x) as.list(sum(dist(x))))) )
#sum of total inter cluster distance
#total.inter.dist <- function(d) sum( do.call('rbind',by(d[,1:2], d[,3], function(x) as.list(sum(dist(x))))) )
#sum of intra cluster correlation
total.intra.cor <- function(d) sum( do.call('rbind',by(d[,1:2], d[,3], function(x) as.list((cor(x)**2)[[2]]))) )
#sum of intra cluster covariance
total.intra.cov <- function(d) sum( do.call('rbind',by(d[,1:2], d[,3], function(x) as.list(det(cov(x))))) )
#likelihood

