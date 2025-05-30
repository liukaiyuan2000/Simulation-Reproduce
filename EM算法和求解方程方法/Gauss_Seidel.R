setwd("F:/advancedmath")
rm(list = ls(all = TRUE))
library(MASS)
A1<-matrix(c(1,2,3,4),nrow=2,byrow=T)
A2<-matrix(c(2,2,3,4,7,7,-2,4,5),nrow=3,byrow=T)
A3<-matrix(c(6,2,1,-1,2,4,1,0,1,1,4,-1,-1,0,-1,3),nrow=4,byrow=T)
A4<-matrix(c(6,2,1,-1,1,2,4,1,0,1,1,1,4,-1,1,-1,0,-1,3,1,1,1,1,1,1),nrow=5,byrow=T)
b1<-c(1,2)
b2<-c(3,1,-7)
b3<-c(6,-1,5,-5)
b4<-c(7,0,6,-4,1)
##xz为正确的解
Gauss_Seidel<-function(A,B,IXH,x,wc){
  H<-nrow(A)
  C<-ncol(A)
  if(H!=C)
  {
    return(paste("Row of matrix is not equal to Column"))
  }
  if(length(B)!=H)
  {
    return(c("length of Right side of the equation is not equal to matrix's Column"))
  }#判断方程两边是否平齐
  if(det(A)==0)
  {
    return(paste("A is singular"))
  }#判断是否奇异
  D<-matrix(nrow = H,ncol = C)
  L<-matrix(nrow = H,ncol = C)
  U<-matrix(nrow = H,ncol = C)
  for(i in 1:H)
  {
    for(j in 1:C)
    {
      if(i==j) { D[i,j]=A[i,j]}
      else{
        D[i,j]=0
      }
      if(i>j)
      {
        L[i,j]<--A[i,j]
      }
      else{L[i,j]=0}
      if(i<j)
      {
        U[i,j]<--A[i,j]
      }
      else{U[i,j]=0}
    }
  
  }
  G<-solve(D-L)%*%U#计算(D-L)^(-1)*U
  pbj<-vector()
  for(i in 1:H)
  {
    if(is.numeric(eigen(G)$values[i]))
    {
      pbj[i]<-eigen(G)$values[i]
    }
    else
    {
      pbj[i]<-Mod(eigen(G)$values[i])
    }
    
  }##如果特征值为复数时计算谱半径
  if(max(pbj)>=1)
  {
    return(paste("It will not converge"))
  }##判断谱半径是否小于1
  for(k in 1:IXH)
 { xp<-x
    for(i in 1:H)
 {a<-B[i]
 if(i==1)
 {
   for(j in (i+1):H)
   {
     a<-a-A[i,j]*x[j]
   }
   x[i]<-a/A[i,i]
 }
 else if(i==H)
 {
   for(j in 1:(i-1) )
   {
     a<-a-A[i,j]*x[j]
   }
   x[i]<-a/A[i,i]
 }
 else
 {
   for(j in 1:(i-1) )
   {
     a<-a-A[i,j]*x[j]
   }
   for(j in (i+1):H)
   {
     a<-a-A[i,j]*x[j]
   }
   x[i]<- a/A[i,i]
 }
 }
    w<-0
    tt<-vector()
    for(i in 1:H)
    {
    tt[i]<-abs(x[i]-xp[i])
    }
    w<-max(tt)
    if(w<=wc)
    {  result <- list('x'=x,'tcxh'=k,'wc'=w)
      return(result)
    }
  }
  result <- list('x'=x,'tcxh'=k,'wc'=w)
  return(result)
}
jj<-5000
x<-c(0,0)
wc=1e-6
Gauss_Seidel(A1,b1,jj,x,wc)
jj<-5000
x<-c(0,0,0)
wc=1e-6
Gauss_Seidel(A2,b2,jj,x,wc)
jj<-5000
x<-c(0,0,0,0)
wc=1e-6
Gauss_Seidel(A3,b3,jj,x,wc)
jj<-5000
x<-c(0,0,0,0,0)
wc=1e-6
Gauss_Seidel(A4,b4,jj,x,wc)

