data("iris")
mydata <- iris[51:150,]
mydata[,5] <- as.vector(mydata[,5])
mydata$Species[mydata$Species=="versicolor"] <- 0
mydata$Species[mydata$Species=="virginica"] <- 1
mydata[,5] <- as.numeric(mydata[,5])
X = as.data.frame(cbind(1,mydata[,1:4]))
names(X) <- c("(Intercept)","Sepal.Length",
              "Sepal.Width","Petal.Length","Petal.Width")
X <- as.matrix(X)
Y = mydata[,5]
p = ncol(X)
n = length(Y)
theta_new = rep(0,p)
theta_old = rep(0,p)
Loop = 100
loop = 0
esp = 1
while(loop<Loop&esp>1e-5){
  loop = loop+1
  theta_old = theta_new
  index = X %*% theta_old
  Temp = exp(-index)/(1 + exp(-index))
  Sn = t(X) %*% (Y - 1 + Temp)
  W = diag(as.vector(Temp*(1 - Temp)))
  #Dn = -t(X) %*% W %*% X
  Dn = -t(X) %*% diag(as.vector(exp(index)/(1 + exp(index))^2))  %*% X
  theta_new = theta_old - solve(Dn) %*% Sn
  esp = sqrt(sum((theta_new - theta_old)^2))
}
theta.hat = theta_new
theta.sd = sqrt(diag(solve(-Dn)))
Zvalue = theta.hat/theta.sd
pvalue = 2*(1-pnorm(abs(Zvalue)))
result <- as.data.frame(cbind(round(theta.hat,3),
          round(theta.sd,3),round(Zvalue,3),round(pvalue,3)))
names(result) <- c("beta.hat","beta.sd","Zvalue","pvalue")
Logisticmodel<-glm(formula=Species~.,family=binomial(link=logit),data=mydata)
result
round(summary(Logisticmodel)$coef,3)