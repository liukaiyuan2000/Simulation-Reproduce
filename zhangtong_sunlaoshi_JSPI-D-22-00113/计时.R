rep = 1
t1 = Sys.time()
for(i in c(1:rep)){
  data=generate(n,diaguu,sigma,delta)
  for(m in c(1:10)){
    for(j in c(1:n)){
      cat(j, m, i, '\r')
      x=data$X[j,]
      y=data$Y[j]
      w=data$W[j,]
      # 删一版本
      datax=data$X[-j,]
      datay=data$Y[-j]
      dataw=data$W[-j,]
      H=pca(datax,m,tm)$h
      eigf=pca(datax,m,tm)$eig.functions
      S=H%*%solve(t(H)%*%H)%*%t(H)
      # (4)式
      betahat=solve(t(dataw)%*%(diag(rep(1,(n-1)))-S)%*%dataw-(n-1)*diaguu)%*%t(dataw)%*%(diag(rep(1,(n-1)))-S)%*%datay
      alphahat=solve(t(H)%*%H)%*%t(H)%*%(datay-dataw%*%betahat)
      alphahatf=eigf%*%alphahat
      # CV函数
      cvresult[i,m]=cvresult[i,m]+(y-t(w)%*%betahat-inpro(alphahatf,x,tm))^2-t(betahat)%*%diaguu%*%betahat
    }
  }
  if(i%%100==0){
    print(i)
  }
}
t2 = Sys.time()

t3 = Sys.time()
cvresult = cross_validation(rep, n, cov.u, delta, tm)
t4 = Sys.time()

print(t2-t1)
print(t4-t3)
