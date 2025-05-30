library(purrr)
DGP <- function(R, m, rho, case.DGP, case.var, d = 0, time.num = 100){
  n = R * m
  time.points = (1:time.num - 0.5) / time.num
  #Generating weight matrix
  I_R <- diag(R)
  I_m <- diag(m)#weight
  l_m <- rep(1, m)#vector
  B_m <- (l_m %*% t(l_m) - I_m) / (m - 1)
  Wn <- I_R %x% B_m
  ginv.mat <- solve(diag(n) - rho*Wn)

  lambda.j = (pi*(1:100 - 0.5))^(-2)
  ## jt
  phi.jt  = sqrt(2)*sin(pi * (1:100 - 0.5) %*% t(time.points))
  switch(
    case.DGP,
    {
      ## DGP1
      Zn = matrix(rnorm(n*2), n, 2)
      gamma.t <- sqrt(2)*sin(pi * time.points / 2) + 3*sqrt(2)*sin(3 * pi * time.points / 2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt     <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      I1 = rowSums(Zn) + colMeans(t(Xnt)*gamma.t)
    },
    {
      ## DGP2
      Zn = matrix(runif(n*2), n, 2)
      gamma.t <- sqrt(2)*sin(pi * time.points / 2) + 3*sqrt(2)*sin(3 * pi * time.points / 2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt    <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      I1 = rowSums(Zn) + colMeans(t(Xnt)*gamma.t)
    },
    {
      ## DGP3
      Zn = matrix(rnorm(n*5), n, 5)
      gamma.t <- sqrt(2)*sin(pi * time.points / 2) + 3*sqrt(2)*sin(3 * pi * time.points / 2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt     <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      temp1 = colMeans(t(Xnt)*gamma.t)
      temp2 = rowSums(Zn)
      I1 = 2*temp2 + temp1 + 0.8 * (cos(0.3*pi*temp1)) + sin(temp2)
    },
    {
      ## DGP4
      Zn = matrix(rnorm(n*2), n, 2)
      gamma.t <- sqrt(2)*sin(pi * time.points / 2) + 3*sqrt(2)*sin(3 * pi * time.points / 2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt     <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      temp1 = colMeans(t(Xnt)*gamma.t)
      I1 = Zn[, 1] + 2*Zn[, 2] + temp1 + 0.35 * (rowSums(Zn) + exp(temp1))
    },
    {
      ## DGP5
      Zn = matrix(rnorm(n*2), n, 2)
      gamma.t <- sqrt(2)*sin(pi * time.points / 2) + 3*sqrt(2)*sin(3 * pi * time.points / 2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt     <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      temp1 = colMeans(t(Xnt)*gamma.t)
      I1 = rowSums(Zn) + temp1 + d*exp(0.1*rowSums(Zn^2) + temp1)
    },
    {
      ## DGP6
      Zn = matrix(runif(n*3), n, 3)
      gamma.t <- sin(pi*time.points/2) + 0.5*sin(3*pi*time.points/2) + 0.25*sin(5*pi*time.points/2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt     <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      temp1 = colMeans(t(Xnt)*gamma.t)
      I1 = rowSums(Zn) + Zn[, 1] + temp1 + d*exp(2*rowSums(Zn^2) + temp1)
    },
    {
      ## DGP7
      Zn = matrix(rnorm(n*2), n, 2)
      gamma.t <- sin(pi*time.points/2) + 0.5*sin(3*pi*time.points/2) + 0.25*sin(5*pi*time.points/2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt     <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      temp1 = colMeans(t(Xnt^2)*gamma.t)
      I1 = 2*rowSums(Zn) + Zn[, 2] + 2*temp1
    },
    {
      ## DGP8
      Zn = matrix(runif(n*3), n, 3)
      gamma.t <- sin(pi*time.points/2) + 0.5*sin(3*pi*time.points/2) + 0.25*sin(5*pi*time.points/2)
      eta.nj  <- map(1:n, \(y) map_dbl(lambda.j, \(x) rnorm(1, 0, sd = sqrt(x))))
      Xnt     <- do.call(rbind, map(eta.nj, \(x) colSums(x * phi.jt)))

      temp1 = colMeans(t(Xnt^2)*gamma.t)
      I1 = 3*rowSums(Zn) + 0.2*temp1 + d*exp(Zn[,1]^2 + Zn[,2]^2 + 2*temp1)
    }
  )
  eps <- switch(
    case.var,
    rnorm(n),
    sqrt(2)*(1 + Zn[, 1])*rnorm(n)
  )
  Yn = ginv.mat %*% (I1 + eps)

  return(
    list(
      Wn = Wn, Zn = Zn, Xnt = Xnt, Yn = Yn
    )
  )
}


phi.hat <- function(Xnt)
{
  t.num = ncol(Xnt)
  n = nrow(Xnt)
  # Calculate phi_hat
  K_hat_st <- t(Xnt) %*% Xnt / (n - 1)
  # Spectral decomposition
  eigen_decomp <- eigen(K_hat_st)
  # Extract eigenvalue
  eigenvalues <- eigen_decomp$values
  # Extract feature vector
  eigenvectors <- eigen_decomp$vectors
  # Calculate the sum of the eigenvalues
  total_variance <- sum(eigenvalues)
  # Sort the eigenvalues by size
  sorted_indices <- order(eigenvalues, decreasing = TRUE)
  sorted_eigenvalues <- eigenvalues[sorted_indices]
  sorted_eigenvectors <- eigenvectors[, sorted_indices]
  # Calculate the cumulative proportion of eigenvalues
  variance_ratio <- cumsum(sorted_eigenvalues) / total_variance
  # Select the eigenvalues with a cumulative proportion of not less than 80%
  selected_indices <- min(which(variance_ratio >= 0.95))
  # sorted_eigenvalues[selected_indices]
  lambda_hat <-  sorted_eigenvalues[1:selected_indices]
  # sorted_eigenvectors[, selected_indices]
  phi <- sorted_eigenvectors[, 1:selected_indices]
  r <- 0.95
  cumulative_variance <- cumsum(sorted_eigenvalues^2) / sum(sorted_eigenvalues^2)
  j_n <- min(which(cumulative_variance >= r))
  eta <- rep(0, j_n)
  h_2_tilde <- rep(0, t.num)
  for(i in 1:j_n)
  {
    eta[i] <- rnorm(1, 0, sqrt(lambda_hat[j_n]))
    h_2_tilde <- h_2_tilde+eta[i]*sorted_eigenvectors[, i]
  }
  #Calculate Phi
  Phi <- Xnt %*% phi / t.num

  return(list(Phi = Phi, h_2_tilde = h_2_tilde))
}


para.hat <- function(Zn, Wn, Yn, Phi)
{
  n = length(Yn)
  # Calculate Qn
  In <- diag(n)
  Qn <- cbind(Wn %*% Yn, Zn)
  P <- Phi %*% ginv(t(Phi) %*% Phi) %*% t(Phi)
  In.minus.P <- In - P
  ######2sls###
  #step1:
  Qnp <- Qn %*% ginv(t(Qn) %*% Qn) %*% t(Qn)
  temp1 <- rbind(t(Wn%*%Yn),t(Zn),t(Phi))
  temp2 <- temp1 %*% Qnp
  para_tilde <- ginv(temp2 %*% t(temp1)) %*% temp2 %*% Yn
  H_tilde <- cbind(Wn %*% ginv(In - para_tilde[1] * Wn) %*% cbind(Phi * para_tilde[4], Zn), Zn)

  #step2:
  M_tilde <- H_tilde %*% ginv(t(H_tilde) %*% H_tilde) %*% t(H_tilde)
  temp3 <- t(Qn) %*% In.minus.P %*% M_tilde %*% In.minus.P
  temp4 <- ginv(t(Phi) %*% Phi) %*% t(Phi)
  theta_bar <- ginv(temp3 %*% Qn) %*% temp3 %*% Yn
  alpha_bar <- temp4 %*% (Yn - Qn %*% theta_bar)
  H <- cbind(Wn %*% ginv(In - theta_bar[1] * Wn) %*% (Phi %*% alpha_bar + Zn %*% theta_bar[-1]), Zn)
  M <- H %*% ginv(t(H) %*% H) %*% t(H)

  temp5 <- t(Qn) %*% In.minus.P %*% M %*% In.minus.P
  theta_hat <- ginv(temp5 %*% Qn) %*% temp5 %*% Yn

  alpha_hat <- temp4 %*% (Yn - Qn %*% theta_hat)

  return(list(theta_hat = drop(theta_hat), alpha_hat = drop(alpha_hat)))
}

TSsingle <- function(Zn, Wn, Yn, Phi, Xnt, theta.hat, alpha.hat, h_2_tilde)
{
  time.num = ncol(Xnt)
  u = (1:time.num - 0.5) / time.num
  n = length(Yn)
  varepsilon.hat <- Yn - theta.hat[1] * Wn %*% Yn - Zn %*% theta.hat[-1] - Phi %*% alpha.hat

  h_mod <- sum(h_2_tilde^2) + sum(theta.hat[-1]^2)

  h <- theta.hat[-1] / sqrt(h_mod)

  h_2 <- h_2_tilde / sqrt(h_mod)

  Uh <- Zn %*% h + Xnt %*% h_2 / time.num

  CRn.hat <- rep(0, time.num)

  for(i in 1:time.num)
  {
    sum <- 0

    for(j in 1:n)
    {
      sum <- sum + varepsilon.hat[j] * ifelse(Uh[j] <= u[i], 1, 0)

    }
    CRn.hat[i] <- sum

  }
  CT_n <- mean(CRn.hat^2)

  return(
    list(
      varepsilon_hat = drop(varepsilon.hat), CT_n = CT_n
    )
  )
}

bootstrap.test <- function(Zn, Wn, Phi, Xnt, theta.hat, alpha.hat, h_2_tilde, varepsilon.hat)
{
  n = nrow(Xnt)
  time.num = ncol(Xnt)
  u.star = (1:time.num - 0.5) / time.num
  #step1: obtain varepsilon.hat and CTn
  #step2:
  varepsilon <- rnorm(n)
  varepsilon.star <- varepsilon.hat * varepsilon

  Y.star <- ginv(diag(n) - theta.hat[1] * Wn) %*% (Zn %*% theta.hat[-1] - Phi %*% alpha.hat + varepsilon.star)

  #step3:
  pavalue <- para.hat(Zn, Wn, Y.star, Phi)
  theta.hat.star <- pavalue$theta_hat
  alpha.hat.star <- pavalue$alpha_hat
  varepsilon.hat.star <- Y.star - theta.hat.star[1] * Wn %*% Y.star - Zn %*% theta.hat.star[-1] - Phi %*% alpha.hat.star

  #step4:
  h.mod <- sum(h_2_tilde^2) + sum((theta.hat[-1])^2)

  h <- theta.hat[-1] / sqrt(h.mod)

  h.2 <- h_2_tilde / sqrt(h.mod)

  Uh <- Zn %*% h + Xnt %*% h.2 / time.num

  CRn.hat.star <- rep(0,time.num)

  for(i in 1:time.num)
  {
    sum <- 0

    for(j in 1:n)
    {
      sum <- sum + varepsilon.hat.star[j] * ifelse(Uh[j] <= u.star[i], 1, 0)
    }
    CRn.hat.star[i] <- sum

  }
  CTn.star <- mean(CRn.hat.star^2)

  return(CTn.star)
}















