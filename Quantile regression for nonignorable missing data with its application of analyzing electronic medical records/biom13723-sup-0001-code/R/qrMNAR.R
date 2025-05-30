## usethis namespace: start
#' @useDynLib qrMNAR, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @export
## ---------------- randomly draw m replicates for each missing response Y from estimated f(Y|X) ---------------
get_ymis.c <- function( n2, X.unobs, tau_seq.ms, beta0_tau.hat, m){
  # 分位数间隔---tau_{k+1}-tau_k
  delta_tau <- tau_seq.ms[2]-tau_seq.ms[1]
  d <- ncol(X.unobs)
  # 抽取随机数
  u_mt <- matrix( runif( n2*m, min = min(tau_seq.ms), max = max(tau_seq.ms)), nrow = n2, ncol = m )
  # 确定分位数的索引
  int.u_mt <- floor( u_mt/delta_tau)
  ymis_mt1 <- matrix( 0, nrow = n2, ncol = m)
  ymis_mt1 <- beta_tau_mt( n2, m, d, delta_tau, as.matrix( beta0_tau.hat), tau_seq.ms, as.matrix( u_mt),
                           as.matrix( int.u_mt), as.matrix( X.unobs), ymis_mt1)
  return(ymis_mt1)
}

#' @export
## ------------------- get the number of levels for categorical variables in X1, X2 ------------------------------
Nlev_catVar <- function(sim_data, p1, q1, p2){
  nlev.catV <- vector(mode="numeric",length=0)  # 初始化一个空向量，用于存储分类变量的水平数
  if(q1==0){  # 如果没有分类变量
    sim_data.cat <- sim_data[, c( 1:p1)]  # 只提取前 p1 列数值型变量
  }else{  # 如果有分类变量
    sim_data.cat <- sim_data[, c( 1:p1, (p2+p1+1):(p2+p1+q1))]  # 提取前 p1 列数值型变量和分类变量
  }
  if( q1+p1 == 1){  # 如果只有一个分类变量
    nlev.catV <- length(levels(factor(sim_data.cat)))  # 计算该变量的水平数
  }else{  # 如果有多个变量
    nlev.catV <- sapply( sim_data.cat, function(x) length(levels(factor(x))))  # 对每个分类变量计算水平数
  }
  return(nlev.catV)  # 返回水平数的向量
}

#### ============================ main function =========================
#' Estimates quantile regression coefficients with nonignorable missing data
#'
#' This function is used to get the proposed estimates of quantile regression coefficients for nonignorable missing data
#'
#' @param sim_data imported data frame of the form sim_data = data.frame(X,Y).
#'  the last column of sim_data is the respones Y which subject to missing (NA if missing).
#'  the other columns of sim_data are the covariates X = (X1, X2)= (X1.cat, X1.cont, X2.cat, X2.cont) in order, where
#'  X1 refers to covariates which are related to the response and also included in the misssing mechanism model;
#'  X2 refers to the auxiliary variable, which means given (X1, Y), X2 is not included in the misssing mechanism model;
#'  X1.cat and X2.cat are the categorical variables in X1 and X2, respectively;
#'  X1.cont and X2.cont are the continuous variables in X1 and X2, respectively.
#' @param p1 the number of categorical variables contained in X1
#' @param p2 the number of continuous variables contained in X1; p1>0 since X1.cont Include intercept term
#' @param q1 the number of categorical variables contained in X2
#' @param q2 the number of continuous variables contained in X2; q1+q2>0 since X2 has a dimension of at least 1
#' @param thresh the threshold of this algorithm for convergence
#' @param max_iter.ms The maximum iterations of this algorithm;
#'  When the algorithm reaches the maximum number of iterations or threshold, the algorithm stops iterating
#' @param tau_seq.out a vector; the quantiles at which the regression coefficients to be estimated for the output
#'  The result will output the coefficient estimates of the quantile regression of these quantile points
#' @param tau_seq.ms a vector; the quantiles at which the regression coefficients to be estimated for f(y|x)
#'  when estimating theta, similar as Wei. et al(2009) to estimate the density f(y|x)
#' @param m the number of samples drawn from f(y|x) for the missing response Y
#'
#' @return the output is a list which include seven elements;
#'  the first element of the output is a vector, including the number of iterations and the maximum
#'  difference of the estimates between the previous and last iterations when the algorithm terminates;
#'  the second element of the output is a vector, the estimates of missing mechanism coefficients;
#'  the third element of the output is a matrix, the estimates of quantile regression coefficients at quantile levels "tau_seq.out"
#'  the fourth element of the output is a matrix, the estimates of missing mechanism coefficients for each iteration
#'  the fifth element of the output is a matrix, the estimates of quantile regression coefficients for each iteration at quantile levels "tau_seq.out".
#'  the sixth element of the output is a list, the summary of missing mechanism model,
#'  more details of missing mechanism model can be found here, including p-values which can be used to determine which variables are significant
#'  the seventh element of the output is a list, the summary of quantile regression model at quantile levels "tau_seq.out".
#'  more details of missing mechanism model can be found here, including p-values which can be used to determine which variables are significant

#' @export
#'
#' @examples
#' ## generate sim_data for example
#' nsize <- 1000
#' beta <- c( 1,-2,2,0.5 ); gam <- c( 0.5, 0.5, 0.2, 0); theta <- c(-2, 0.5,0.5,0.5)
#' x1 <- rbinom( nsize, 1, 0.5); x2 <- rnorm( nsize, mean = 2, sd = 0.5); x3 <- runif( nsize, min = 0, max = 2)
#' xx <- cbind( 1, x1, x2, x3); err <- xx%*%gam*rnorm( nsize, mean = 0, sd = 0.5); yy <- xx %*% beta + err
#' xy <- cbind( 1, x1, x2, yy); prob.m <- 1/( 1+ exp( - xy %*% theta ))
#' delta.ind <- sapply( prob.m, function(p){rbinom(n=1, size=1, prob=p)}); yy[delta.ind==0] <- NA
#' sim_data <- data.frame( x1, x2, x3, yy)
#' p1<- 1; p2<- 1; q1<- 0; q2 <- 1
#'
#' ## set the initial value related to the algorithm
#' thresh = 1e-04; max_iter.ms <- 50;  m <- 100
#' tau_seq.ms  <- seq(0.05, 0.95, by = 0.05); tau_seq.out <- seq(0.25, 0.75, by = 0.25)
#'
#' ## using sim_data to estimate quantile coefficients
#' out <- msIpwQr(sim_data, p1, p2, q1, q2, thresh, max_iter.ms, tau_seq.out, tau_seq.ms, m )
#' beta_tau_true <- beta + gam %*% t(qnorm( tau_seq.out, mean = 0, sd = 0.5)) # true value of beta_tau
#' out[[3]] # quantile coefficients estimates
#' out[[7]] # summary of quantile regression model
msIpwQr <- function(sim_data, p1, p2, q1, q2,  thresh, max_iter.ms, tau_seq.out , tau_seq.ms, m ){
  #convert cat. variables in X to dummy variable; including intercept
  colnames(sim_data)[ncol(sim_data)] <- "Y"
  delta.ind <- ifelse(is.na(sim_data$Y) == 1, 0, 1) # delta.ind =1 if response Y is obs'd; o.w., delta.ind = 0
  X.name <- colnames(sim_data)[-ncol(sim_data)]
  # 生成设计矩阵
  X <- data.frame(model.matrix(as.formula(paste("~", paste (X.name ,collapse = "+"))),data = sim_data))
  X <- X[, which(apply(X[delta.ind==1, ],2,sum)>0)]

  sim_data1 <- data.frame( X, Y = sim_data[, ncol(sim_data)])
  sim_data1.obs <- sim_data1[delta.ind == 1, ]
  sim_data1.unobs <- sim_data1[delta.ind == 0, ]
  X.obs <- X[delta.ind==1,];  X.unobs <- X[delta.ind==0,]
  nsize <- nrow(sim_data); n1 <- nrow(X.obs) ; n2 <- nsize - n1
  l <- length(tau_seq.ms)
  nlev.catV  <- vector(mode="numeric",length=0)
  if( p1+q1 >=1){
    nlev.catV <- Nlev_catVar(sim_data, p1, q1, p2) # get the num of levels of category variables in covariates X1, X2
  }
  if( p1>=1){
    X1.obs <- X.obs[ , 1:( sum(nlev.catV[1:p1])-p1 + p2 +1)] # including intercept
    X1.unobs<- X.unobs[ , 1:( sum(nlev.catV[1:p1])-p1 + p2 +1)]
  }else{
    X1.obs <- X.obs[, 1:(p2 +1)]
    X1.unobs<- X.unobs[ ,1:(p2 +1)]
  }
  ## obtain initial values for parameters beta and thata
  beta0_tau.hat <- rq(Y ~ -1+ ., tau = tau_seq.ms , data =sim_data1.obs)$coefficients
  theta0_hat <- rep( 0, ncol(X1.obs) + 1 )
  ## to keep updated values for parameters
  beta1_tau.hat <- matrix(10, nrow = nrow(beta0_tau.hat), ncol = ncol(beta0_tau.hat ))
  theta1_hat <- rep(10, length(theta0_hat))
  beta0_tau.hat.out <-  matrix( 0, nrow = nrow(beta0_tau.hat), ncol = length(tau_seq.out))
  beta1_tau.hat.out <-  matrix( 10, nrow = nrow(beta0_tau.hat), ncol = length(tau_seq.out))
  beta_hat.iters <- rq(Y ~ -1+ ., tau = tau_seq.out , data =sim_data1.obs)$coefficients
  theta_hat.iters <- matrix(0, ncol = length(theta0_hat) )
  ymis_mt <- matrix( 0, nrow = n2, ncol =  m)
  sw1 <- matrix( 0, nrow = n2, ncol = m)
  sw2 <- rep( 0, n2)
  sw <- matrix(0, nrow = n2, ncol = m)
  iter <- 0
  max_diff <- 10
  while( iter <=  max_iter.ms & ( max_diff >= thresh) ){
    # using mean score function to update theta with beta0_tau.hat and obtain m replicates of missing Y at current parameter values
    theta1_hat <- theta0_hat
    # 抽取缺失部分的样本
    ymis_mt <- get_ymis.c( n2, X.unobs, tau_seq.ms, beta0_tau.hat, m)
    exp1.unobs <-  matrix( as.matrix( X1.unobs) %*% theta0_hat[c(1:ncol(X1.unobs))], nrow = n2, ncol = m, byrow = F)
    sw1 <- 1/(1+exp(exp1.unobs + theta0_hat[ncol(X1.unobs)+1]*ymis_mt))
    # 这里是remark 2中的权重
    sw <- sw1/matrix( apply(sw1, 1, sum),  nrow = n2, ncol = m, byrow = F)
    df <- rbind( X1.obs, X1.unobs[rep(1:n2, m), ])
    df$Y <- c( sim_data1.obs$Y,  as.vector(ymis_mt))
    df$delta <- c(rep(1,n1), rep(0,n2*m))
    miss.out <- try( glm( delta ~ -1 + ., family=binomial, data = df, weights = c( rep( 1, n1), as.vector(sw))), TRUE)
    fit <- try( miss.out$coefficient, TRUE)
    if((is.numeric(fit) == 0) | (NaN %in% fit) | (Inf %in% fit)  | (-Inf %in% fit)  | (NA %in% fit)){
      print(fit)
      break
    }else if(max(abs(fit))>= 1000){
      print(fit)
      break
    }else{
      theta0_hat <- fit
    }
    ## using IPW+QR method to update beta_tau
    wi <- (1+exp(- as.matrix( cbind(X1.obs, sim_data1.obs$Y)) %*% theta0_hat))
    beta1_tau.hat <- beta0_tau.hat
    beta0_tau.hat <- rq(Y ~ -1 + . , data = sim_data1.obs, tau = tau_seq.ms, weights = wi)$coefficients

    beta1_tau.hat.out <- beta0_tau.hat.out
    rq.out <- rq(Y ~ -1 + . , data = sim_data1.obs, tau = tau_seq.out, weights = wi)
    beta0_tau.hat.out <- rq.out$coefficients
    beta_hat.iters <- rbind( beta_hat.iters, beta0_tau.hat.out)
    theta_hat.iters <- rbind( theta_hat.iters, theta0_hat)

    iter <- iter+1
    max_diff <- mean(abs(beta0_tau.hat - beta1_tau.hat))
    cat("iter = ", iter, sep="", "\n")
    cat("diff.beta = ", max_diff , '\n')
  }
  out <- list( n.iter_diff = c(iter-1, max_diff),
               theta_hat = c(theta0_hat), beta_tau.msIpwQr = data.frame(beta0_tau.hat.out),
               theta_hat.iters = theta_hat.iters, beta_hat.iters = beta_hat.iters,
               summary.miss = summary(miss.out), summary.rq =  summary(rq.out,  se="boot"))
  print("msIpwQr Ends!")
  return(out)
}

.onUnload <- function (libpath) { library.dynam.unload("qrMNAR", libpath)}

