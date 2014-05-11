library(mvtnorm)
n <- 100
X <- as.data.frame(rmvnorm(N,mean=c(0,0)))
colnames(X) <- c('x','y')
attach(X)
p <- ncol(X)
dof <- n-p
m<-lm(y ~ x, data=X)

y.obs <- X$y
y.exp <- predict(m,list(x=X$x))

RSS.m/SSE 
RSS.m/RSS.0

mean(sqrt((y.obs-y.exp)**2))

cov(y,x)/var(y)

#pearson correlation
cov(x,y)/(sd(x)*sd(y))

#parameters
beta <- cov(y,x)/var(x)
beta.std.err <- sqrt(SSE/98)/sqrt(sum((x-mean(x))**2))
cat('beta', beta, '\n')
cat('beta.std.err', beta.std.err, '\n')
alpha <- mean(y)-beta*mean(x)
alpha.std.err <- beta.std.err*sqrt(sum(x**2)/n)
cat('alpha', alpha, '\n')
cat('apha.std.err', alpha.std.err, '\n')

#Residual standard error: 1.082 on 98 degrees of freedom
cat('Residual standard error:', rse <- sqrt(sum((y.obs-y.exp)**2)/dof), 'on', dof, '\n')

#multiple R^2
multiple.rsquared <- cov(y.obs,y.exp)/cov(y.obs,y.obs)
cat('Multiple R-squared:',round(multiple.rsquared <- cov(y.obs,y.exp)/var(y.obs),5),'\n')

#adjusted R^2

total.SS <- RSS.0 <- sum((y.obs-mean(y.obs))**2)
RSS.m <- explained.SS <- sum((y.exp-mean(y.obs))**2)
SSE <- unexplained.SS <- sum((y.obs-y.exp)**2)
variance.explained <- explained.SS/(p-1)
variance.unexplained <- unexplained.SS/dof 

#F-statistic
cat('F-statistic', F.statistic <- variance.explained/variance.unexplained)
cat(' on', p-1, 'and', n-p, 'DF', '\n')
cat('p-value', 1-pf(F.statistic,p-1,n-p), '\n')

print(summary(m))

predict(m,list(x=mean(x))) == mean(y)


predict(m,list(x=0))

