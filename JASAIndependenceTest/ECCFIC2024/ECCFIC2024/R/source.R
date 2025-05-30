########################################################################################
### Expected Conditional Characteristic Function-based Measures for Testing Independence
### (C. Ke and X. Yin, 2019)
### R code author: Chenlu Ke (chenlu.ke@uky.edu)
########################################################################################

#---------------------------------------------------------------------------------------------
#' Slicing method
#' @name H2d
#' @description compute the ECCFIC between x and y, where y is categorical
#' 
#' @param x a numeric vector/matrix
#' @param y a numeric vector
#' @param kernel a kernel to use for x, 'gaussian' or 'distance'
#' @param d exponent on distance in (0,2], if a 'distance' kernel is used
#' @param sigma bandwidth, if a 'gaussian' kernel is used, 
#'              default is heuristic median pairwise distances of x
#' @usage H2d(x, y, kernel='gaussian', d=1, sigma='default')
#' @examples
#' y <- rep(1:4, 25)
#' x <- y + rnorm(100)
#' result <- H2d(x,y)

H2d <- function(x, y, kernel='gaussian', d=1, sigma='default'){
  
  x = as.matrix(x)
  n = nrow(x)
  
  #n=1, trivial
  if(n==1){
    return(0)
  }
  
  #compute heuristic sigma
  if(kernel=='gaussian' & sigma=='default'){
    sigma = sqrt(0.5*median(dist(x)^2))
  }
  
  #compute H2d(x,x)
  if(identical(x,as.matrix(y))){
    if(kernel=='gaussian'){
      K = dnorm(as.matrix(dist(x,diag=T,upper=T)), mean=0, sd=sigma)
      return((dnorm(0,0,sigma)-mean(K))*sqrt(2*pi)*sigma)
    }else{
      if(kernel=='distance'){
        K = 0.5*(as.matrix(dist(x,diag=T,upper=T))^d)
        return(mean(K))
      }
    }
  }
  
  #compute ECCFIC
  if(kernel=='gaussian'){
    K = dnorm(as.matrix(dist(x,diag=T,upper=T)), mean=0, sd=sigma)
    L = ifelse(as.matrix(dist(y,diag = T,upper = T))==0,1,0)
    L = n*L/rowSums(L)-1
    return(sum(K*L)/(n^2)*sqrt(2*pi)*sigma)
  }else{
    if(kernel=='distance'){
      K = -0.5*(as.matrix(dist(x,diag=T,upper=T))^d)
      L = ifelse(as.matrix(dist(y,diag = T,upper = T))==0,1,0)
      L = n*L/rowSums(L)-1
      return(sum(K*L)/(n^2))
    }
  }
}


#' Permutation test with slicing estimator
#' @name ptest_H2d
#' @param B number of replicates
#' @examples
#' y <- rep(1:4, 25)
#' x <- y + rnorm(100)
#' result <- ptest_H2d(x,y,B=199)

ptest_H2d <- function(x,y,kernel='gaussian',d=1,sigma='default',B){
  H = H2d(x,y,kernel,d,sigma)
  p = (sum(replicate(B,H2d(x,sample(y),kernel,d,sigma))>=H)+1)/(B+1)
  return(list(ECCFIC=H,pvalue=p))
}


#---------------------------------------------------------------------------------------------
#' Kernel regression method
#' @name H2c
#' @description compute the ECCFIC between continuous x and y, 
#'              need to load the library kernlab
#' 
#' @param x a numeric vector/matrix
#' @param y a numeric vector/matrix
#' @param kernel a kernel to use for x, 'gaussian' or 'distance'
#' @param d exponent on distance in (0,2], if a 'distance' kernel is used for x
#' @param sigma bandwidth, if a 'gaussian' kernel is used for x, 
#'              default is heuristic median pairwise distances of x
#' @param bw bandwidth of the gaussian smoothing kernel applied on y,
#'           bandwidths suggested by Silverman (1986) are used unless otherwise specified.
#' @param wt use weight function or not, if TRUE, f^2 is used by default
#' @usage H2c(x, y, kernel='gaussian', d=1, sigma='default', bw='default', wt=F)
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100)
#' result <- H2c(x,y)

H2c <- function(x, y, kernel='gaussian', d=1, sigma='default', bw='default', wt=F){
  
  x = as.matrix(x)
  n = nrow(x)
  
  #n=1, trivial
  if (n==1){
    return(0)
  }
  
  #compute heuristic sigma for x
  if(sigma=='default'){
    sigma = sqrt(0.5*median(dist(x)^2))
  }
  

  y = as.matrix(y)
  p = ncol(y)
  
  #compute default bandwidth
  if(bw=='default'){
    if(p==1){
      bw = 1.06*sd(y)*n^(-1/5)
    }else{
      bw = (4/n/(p+2))^(1/(p+4))*sum(diag(cov(y)))/p
    }
  }
  
  H = diag(n) - matrix(1, n, n)/n
  if(kernel=='gaussian'){
    #Incomplete Cholesky decomposition of HKH
    tempx = inchol(as.matrix(x), kernel = 'rbfdot', kpar = list(sigma = 1/2/(sigma^2)))
    iCHx = H%*%tempx@.Data
  }else{
    if(kernel=='distance'){
      #Cholesky decomposition of HKH
      normx = matrix(apply(x*x,1,sum)^(d/2),n,n)
      tempx = chol(0.5*(normx+t(normx)-as.matrix(dist(x,diag=T,upper=T))^d))
      iCHx = H%*%tempx
    }
  }
  a = svd(iCHx)
  Ux = a$u
  Sx2 = (a$d)^2
  
  if(wt==T){
    #HGGH
    tempy = dnorm(as.matrix(dist(y,diag=T,upper=T)), mean=0, sd=bw)
    iCHy = H%*%tempy
    b = svd(iCHy)
    Uy = b$u
    Sy2 = (b$d)^2
    #ECCFIC = trace(KHGGH)/n^3
    return(sum(diag(Ux%*%diag(Sx2)%*%t(Ux)%*%Uy%*%diag(Sy2)%*%t(Uy)))/n^3)
  }else{
    #no weight
    tempy = dnorm(as.matrix(dist(y,diag=T,upper=T)), mean=0, sd=bw)
    tempy = t(tempy/rowSums(tempy))
    iCHy = H%*%tempy
    b = svd(iCHy)
    Uy = b$u
    Sy2 = (b$d)^2
    return(sum(diag(Ux%*%diag(Sx2)%*%t(Ux)%*%Uy%*%diag(Sy2)%*%t(Uy)))/n)
  }
}

#' Permutation test with kernel regression estimator
#' @name ptest_H2c
#' @param B number of replicates
#' @examples
#' x <- rnorm(100)
#' y <- x + rnorm(100)
#' result <- ptest_H2c(x,y,B=199)

ptest_H2c <- function(x,y,kernel='gaussian',d=1,sigma='default',bw='default',wt=F,B){
  n = nrow(as.matrix(x))
  H = H2c(x,y,kernel,d,sigma,bw,wt)
  y = as.matrix(y)
  p = (sum(replicate(B,H2c(x,y[sample(n),],kernel,d,sigma,bw,wt))>=H)+1)/(B+1)
  return(list(ECCFIC=H,pvalue=p))
}
