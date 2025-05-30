library(MASS)
n <- 200
beta0 <- 1

p <- length(beta0)
Sigma <- matrix(0,p,p)
Sigma <- 0.5^(abs(col(Sigma)-row(Sigma)))

X <- mvrnorm(n,rep(0,p),Sigma)
Ystar <- X%*%beta0+rnorm(n)
C <- runif(n,0,2)
Y <- pmin(C,Ystar)
delta <- 1*(Ystar<C)
mean(delta)
data1 =list(X=X,Y=Y,delta=delta)

fit1 <- lm(Y~X+0)
summary(fit1)
AY <- Ystar[delta==0]
AX <- X[delta==0]
BY <- Y[delta==0]
BX <- X[delta==0]
CY <- Y[delta==1]
CX <- X[delta==1]
fit2 <- lm(CY~CX+0)
beta.naive2 <- fit2$coef
beta.naive1 <- fit1$coef
sigma.naive1 <- sqrt(sum(fit1$residuals^2)/(n-p))
sigma.naive2 <- sqrt(sum(fit2$residuals^2)/(n-p))
theta.naive1 <- c(beta.naive1,sigma.naive1^2)
theta.naive2 <- c(beta.naive2,sigma.naive2^2)

fit <- lm(Ystar~X+0)




par(mfrow=c(1,2))
plot(X,Ystar,xlim=c(-2,2),ylim=c(-4,5),xlab="X",ylab="Y")
lines(X,fitted(fit),col=2)


plot(CX,CY,xlim=c(-2,2),ylim=c(-4,5),xlab="X",ylab="Y")
points(AX,AY,pch=2)
points(BX,BY,pch=20)
lines(X,fitted(fit),col=2)
lines(X,fitted(fit1),col=4)
par(mfrow=c(1,1))

