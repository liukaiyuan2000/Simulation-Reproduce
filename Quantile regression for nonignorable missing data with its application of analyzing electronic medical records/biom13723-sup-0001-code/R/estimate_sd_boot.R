## usethis namespace: start
#' @useDynLib qrMNAR, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' @export
#### ---------------------- Function for weighted bootstrap ---------------------------------------
msIpwQr.wiboot <- function(sim_data, p1, p2, q1, q2,  thresh, max_iter.ms, tau_seq.out , tau_seq.ms, m ){
  #convert cat. variables in X to dummy variable; including intercept
  colnames(sim_data)[ncol(sim_data)] <- "Y"
  delta.ind <- ifelse(is.na(sim_data$Y) == 1, 0, 1) # delta.ind =1 if response Y is obs'd; o.w., delta.ind = 0
  X.name <- colnames(sim_data)[-ncol(sim_data)]
  # 构造设计矩阵
  X <- data.frame(model.matrix(as.formula(paste("~", paste (X.name ,collapse = "+"))),data = sim_data)) # create indicate var. for cat. var.
  X <- X[, which(apply(X[delta.ind==1, ],2,sum)>0)] #ensure the matrix of X is not Singular matrix
  # Y <- Y[which(apply(X[delta.ind==1, ],2,sum)>0)]

  sim_data1 <- data.frame( X, Y = sim_data[, ncol(sim_data)])
  sim_data1.obs <- sim_data1[delta.ind == 1, ]
  sim_data1.unobs <- sim_data1[delta.ind == 0, ]
  X.obs <- X[delta.ind==1,];  X.unobs <- X[delta.ind==0,]
  nsize <- nrow(sim_data); n1 <- nrow(X.obs) ; n2 <- nsize - n1

  # generate bootstrapping weights
  wi_boot <- rexp( nsize)
  wi_boot.obs <-  wi_boot[ delta.ind==1]
  wi_boot.unobs <-  wi_boot[ delta.ind==0]

  l <- length( tau_seq.ms)
  if(p1+q1 >=1){
    nlev.catV <- Nlev_catVar(sim_data, p1, q1, p2) # get the num of levels of category variables in covariates X1, X2
  }
  if(p1>=1){
    X1.obs <- X.obs[, 1:( sum(nlev.catV[1:p1]) - p1 + p2 + 1)] # including intercept
    X1.unobs<- X.unobs[ ,1:( sum(nlev.catV[1:p1]) - p1 + p2 + 1)]
  }else{
    X1.obs <- X.obs[, 1:(p2 +1)]
    X1.unobs<- X.unobs[ ,1:(p2 +1)]
  }
  ## obtain initial values for parameters beta and thata
  beta0_tau.hat <- rq(Y ~ -1+ ., tau = tau_seq.ms , data =sim_data1.obs)$coefficients
  theta0_hat <- rep( 0, ncol(X1.obs) + 1)
  ## to keep updated values for parameters
  beta1_tau.hat <- matrix( 10, nrow = nrow(beta0_tau.hat), ncol = ncol(beta0_tau.hat ))
  theta1_hat <- rep( 10, length(theta0_hat))
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
    ymis_mt <- get_ymis.c( n2, X.unobs, tau_seq.ms, beta0_tau.hat, m)
    exp1.unobs <- matrix( as.matrix(X1.unobs)%*% theta0_hat[c(1:ncol(X1.unobs))], nrow = n2, ncol = m, byrow = F)
    sw1 <- 1/(1+exp(exp1.unobs + theta0_hat[ncol(X1.unobs)+1]*ymis_mt))
    sw <- sw1/matrix( apply(sw1, 1, sum),  nrow = n2, ncol = m, byrow = F)
    df <- rbind( X1.obs, X1.unobs[rep(1:n2, m), ])
    df$Y <- c( sim_data1.obs$Y,  as.vector(ymis_mt))
    df$delta <- c( rep(1,n1), rep(0,n2*m))
    # 这里的glm中乘了bootstrap的权重
    fit <- try(glm( delta ~ -1 + ., family=binomial, data = df, weights = c( wi_boot.obs,
           as.vector(sw)*rep( wi_boot.unobs, m )))$coefficient, TRUE)
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
    # 逆概率权重
    wi <-  (1+exp(- as.matrix( cbind(X1.obs, sim_data1.obs$Y)) %*%  theta0_hat))*wi_boot.obs
    beta1_tau.hat <- beta0_tau.hat
    beta0_tau.hat <- rq(Y ~ -1 + . , data = sim_data1.obs, tau =  tau_seq.ms, weights = wi)$coefficients
    beta1_tau.hat.out <- beta0_tau.hat.out
    beta0_tau.hat.out <- rq(Y ~ -1 + . , data = sim_data1.obs, tau = tau_seq.out, weights = wi)$coefficients
    beta_hat.iters <- rbind( beta_hat.iters, beta0_tau.hat.out)
    theta_hat.iters <- rbind( theta_hat.iters, theta0_hat)

    iter <- iter+1
    max_diff <- mean(abs(beta0_tau.hat - beta1_tau.hat))
    cat("iter = ", iter, sep="", "\n")
    cat("diff.beta = ", max_diff , '\n')
  }
  out <- list(n.iter_diff = c( iter-1, max_diff), theta_hat = c(theta0_hat), beta_tau.msIpwQr = data.frame(beta0_tau.hat.out),
              theta_hat.iters = theta_hat.iters, beta_hat.iters = beta_hat.iters )
  return(out)
}

####==================== using bootstrap method to estimate variance by parallel computing ==============================
#' Estimate standard deviation of the proposed estimates via bootstrapping
#'
#' This function is used to estimate standard deviation of the proposed estimates via weighted bootstrapping by parallel computing
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
#' @param B bootstrap times to estimate estimate standard deviation of the proposed coefficient estimate
#' @param num_cores number of CPU cores used for parallel computing for bootstrapping
#'
#' @return  the output is a list; the first element of the list is a vector, including the average number of iteration of
#' bootstrapping B times and the average maximum difference of the estimates between the previous and last iterations
#' when the algorithm terminates;
#' the second element of the list is a vector, the estimated standard deviation of missing mechanism coefficients
#' estimation via bootstrapping B times;
#' the third element of the list is a matrix, the estimated standard deviation of quantile regression coefficients
#' estimation at quantile "tau_seq.out" via bootstrapping B times.
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
#' thresh = 1e-04; max_iter.ms <- 50;  m <- 100; B = 100; num_cores <- 5
#' tau_seq.ms  <- seq(0.05, 0.95, by = 0.05); tau_seq.out <- seq(0.25, 0.75, by = 0.25)
#'
#' ## to obtain standard errors of quantile coefficients estimates via bootstrapping by parallel computing
#' out_sd <- estimate_sd.boot( sim_data, p1, p2, q1, q2, thresh, max_iter.ms, tau_seq.out, tau_seq.ms, m, B, num_cores)
#' out_sd[[3]] # standard errors of quantile coefficients estimates via bootstrapping
estimate_sd.boot <- function(sim_data, p1, p2, q1, q2, thresh, max_iter.ms, tau_seq.out, tau_seq.ms, m, B, num_cores){
  cl <- snow::makeSOCKcluster(num_cores)
  doSNOW::registerDoSNOW(cl)
  # snow::clusterEvalQ(cl, source(file.path(getSrcDirectory(function(x) {x}), "qrMNAR1.R")))
  varlist <- c( 'sim_data', 'p1', 'p2', 'q1', 'q2', 'thresh', 'max_iter.ms', 'tau_seq.out',
                 'tau_seq.ms', 'm', 'get_ymis.c', 'Nlev_catVar', 'msIpwQr.wiboot')
  snow::clusterExport( cl, varlist)
  print("Show the Progress of Bootstrapping:")
  pb <- txtProgressBar( min = 0, max = B, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list( progress = progress)
  result <- foreach::foreach( i=1:B, .options.snow= opts,.combine = rbind, .packages = c( "quantreg"), .noexport = c("beta_tau_mt")) %dopar% {
    msIpwQr.wiboot( sim_data, p1, p2, q1, q2, thresh, max_iter.ms, tau_seq.out, tau_seq.ms, m )[1:3]
  }
  close(pb)
  snow::stopCluster(cl)

  rows <- nrow(result[[2*B+1]])
  len_theta <- length(result[[B+1]])
  d.iter <- length(result[[1]])
  iter_diff_mt <- matrix(0, nrow = B, ncol =  d.iter)
  theta_mt.basic <- matrix(0, nrow = B, ncol = len_theta)
  beta_tau_mt.basic <- matrix(0, nrow = B*rows, ncol = length(tau_seq.out) )
  for(i in 1:B){
    iter_diff_mt[i, ] <-  as.matrix(result[[i]])
    theta_mt.basic[i,] <- c(result[[B+i]])
    beta_tau_mt.basic[((i-1)*rows+1):(i*rows),] <- as.matrix(result[[2*B+i]])
  }

  ese_beta_tau.basic <- matrix(0, nrow = rows, ncol = length(tau_seq.out))
  for(i in 1:rows){
    ese_beta_tau.basic[ i, ] <- apply( beta_tau_mt.basic[ seq(i, B*rows, by = rows), ], 2, sd)
  }
  ese_theta.basic <- apply(theta_mt.basic, 2, sd)
  mean_iter_diff <- apply( iter_diff_mt, 2, mean)

  return( list( mean_iter_diff, ese_theta.basic, ese_beta_tau.basic))
}

.onUnload <- function (libpath) { library.dynam.unload("qrMNAR", libpath)}
