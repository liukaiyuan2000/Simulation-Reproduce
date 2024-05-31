##### New test statistics and their wild bootstrap counterparts ----

M.tests = function(X,Y,n,epsilon = 0.0001){
  
  F = ecdf(X)
  G = ecdf(Y)
  Z = c(X,Y)
  E = F(Z)-G(Z)
  KS = sqrt(n)*max(abs(E))
  
  H = ecdf(pmax(X,Y))
  w = sqrt( abs( F(Z)+G(Z)-2*H(Z)-(F(Z)-G(Z))^2 ) )
  wKS = sqrt(n)*max(abs(ifelse(w == 0,0,E/w)))
  
  Z.t = sort(Z)[ceiling(epsilon*2*n):floor((1-epsilon)*2*n)]
  E.t = F(Z.t)-G(Z.t)
  w.t = sqrt( abs(F(Z.t)+G(Z.t)-2*H(Z.t)-(F(Z.t)-G(Z.t))^2 ) )
  WKS = sqrt(n)*max(abs(ifelse(w.t == 0,0,E.t/w.t)))
  
  M = max(KS,WKS)
  
  all = c(KS,wKS,WKS,M)
  return(all)
  
}

M.tests.boot = function(X,Y,n,epsilon = 0.0001,vec.Xi, par.c){
  
  Z = c(X,Y)
  Z.s = sort(Z)
  E.boot = matrix(0,1,2*n)
  for(i in 1:(2*n)){
    E.boot[1,i] = sum( vec.Xi*( ifelse(X <= Z.s[i],1,0) - ifelse(Y <= Z.s[i],1,0) ) )/sqrt(n)
  }
  KS.boot = max(abs(E.boot))
  
  w.boot = matrix(0,1,2*n)
  for(i in 1:(2*n)){
    w.boot[1,i] = sum( abs(vec.Xi/sqrt(2/pi))*( ifelse(X <= Z.s[i],1,0) + ifelse(Y <= Z.s[i],1,0) - 2*ifelse(pmax(X,Y) <= Z.s[i],1,0) ) )/n -
      (sum( abs(vec.Xi/sqrt(2/pi))*(ifelse(X <= Z.s[i],1,0) - ifelse(Y <= Z.s[i],1,0)) ))^2/n^2 
  }
  wKS.boot = max(abs(ifelse(w.boot == 0,0,E.boot/sqrt(abs(par.c*w.boot)))))
  
  Z.t = sort(Z)[ceiling(epsilon*2*n):floor((1-epsilon)*2*n)]
  dl = length(Z.t)
  E.t.boot = matrix(0,1,dl)
  w.t.boot = matrix(0,1,dl)
  for(i in 1:dl){
    E.t.boot[1,i] = sum( vec.Xi*( ifelse(X <= Z.t[i],1,0) - ifelse(Y <= Z.t[i],1,0) ) )/sqrt(n)
    w.t.boot[1,i] = sum( abs(vec.Xi/sqrt(2/pi))*( ifelse(X <= Z.t[i],1,0) + ifelse(Y <= Z.t[i],1,0) - 2*ifelse(pmax(X,Y) <= Z.t[i],1,0) ) )/n -
      (sum( abs(vec.Xi/sqrt(2/pi))*(ifelse(X <= Z.t[i],1,0) - ifelse(Y <= Z.t[i],1,0)) ))^2/n^2 
  }
  
  WKS.boot = max(abs(ifelse(w.t.boot == 0,0,E.t.boot/sqrt(abs(par.c*w.t.boot)))))
  
  M.boot = max(KS.boot,WKS.boot)
  
  all.boot = c(KS.boot,wKS.boot,WKS.boot,M.boot)
  return(all.boot)
  
}

###################################################