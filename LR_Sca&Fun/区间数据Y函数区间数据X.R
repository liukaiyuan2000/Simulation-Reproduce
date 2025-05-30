#install.packages(c("fda", "MASS"))
library(fda)
library(MASS)

# Function to estimate coefficient matrix B

# response            : a matrix containing response curves (center)
# response_l          : a matrix containing lower limits of the response curves
# response_u          : a matrix containing upper limits of the response curves
# response_r          : a matrix containing half-ranges of the response curves
# predictor           : a list whose elements consist of predictors (center)
# predictor_l         : a list whose elements consist of lower limits of the predictors
# predictor_u         : a list whose elements consist of upper limits of the predictors
# predictor_r         : a list whose elements consist of half-ranges of the predictors
# predictor_test      : a list whose elements consist of predictors (center) for the test sample
# predictor_test_l    : a list whose elements consist of lower limits of the predictors for the test sample
# predictor_test_u    : a list whose elements consist of lower upper of the predictors for the test sample
# predictor_test_r    : a list whose elements consist of half-ranges of the predictors for the test sample
# nbf_vec_predictors  : a vector of number of basis function corresponding to predictors
# num_t               : Number of discrete time points
# np                  : Number of predictors
# n_train             : Number of functions in the training samples
iv_fof = function(response, response_l, response_u, response_r, 
                  predictor, predictor_l, predictor_u, predictor_r,
                  predictor_test, predictor_test_l, predictor_test_u, predictor_test_r,
                  nbf_vec_predictors, num_t, np, n_train)
{
  # Discrete time points
  dtp = 1:num_t
  # B-spline for predictors
  B_spline_basis_x = vector("list",)
  B_spline_basis_funs_x = vector("list",)
  
  for(i in 1:np){
    B_spline_basis_x[[i]] = 
      create.bspline.basis(c(1, num_t), nbasis = nbf_vec_predictors[i])
    B_spline_basis_funs_x[[i]] = eval.basis(dtp, B_spline_basis_x[[i]])
  }
  
  #Inner products
  Inner_prod = vector("list",)
  
  for(i in 1:np)
    Inner_prod[[i]] = inprod(B_spline_basis_x[[i]], B_spline_basis_x[[i]])
  
  # Weight argument
  w_arg = matrix(dtp, nrow = n_train, ncol = num_t, byrow = T)
  
  # Weight matrices of the predictors
  W_x = vector("list",)
  W_xl = vector("list",)
  W_xu = vector("list",)
  W_xr = vector("list",)
  
  W_x_test = vector("list",)
  W_xl_test = vector("list",)
  W_xu_test = vector("list",)
  W_xr_test = vector("list",)
  
  for(i in 1:np){
    W_x[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xl[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_l[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xu[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_u[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xr[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_r[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    
    W_x_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xl_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test_l[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xu_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test_u[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
    W_xr_test[[i]] = t(smooth.basis(argvals=t(w_arg), y=t(predictor_test_r[[i]]), fdParobj=B_spline_basis_x[[i]])$fd$coefs)
  }
  
  # Mean of variables
  mean_y = mean(response)
  mean_yl = mean(response_l)
  mean_yu = mean(response_u)
  mean_yr = mean(response_r)
  
  mean_x = vector("list",)
  mean_xl = vector("list",)
  mean_xu = vector("list",)
  mean_xr = vector("list",)
  
  mean_x_test = vector("list",)
  mean_xl_test = vector("list",)
  mean_xu_test = vector("list",)
  mean_xr_test = vector("list",)
  
  for(i in 1:np){
    mean_x[[i]] = apply(W_x[[i]], 2, mean)
    mean_xl[[i]] = apply(W_xl[[i]], 2, mean)
    mean_xu[[i]] = apply(W_xu[[i]], 2, mean)
    mean_xr[[i]] = apply(W_xr[[i]], 2, mean)
    
    mean_x_test[[i]] = apply(W_x_test[[i]], 2, mean)
    mean_xl_test[[i]] = apply(W_xl_test[[i]], 2, mean)
    mean_xu_test[[i]] = apply(W_xu_test[[i]], 2, mean)
    mean_xr_test[[i]] = apply(W_xr_test[[i]], 2, mean)
  }
  
  # Centered variables
  W_y <- response - mean_y
  W_yl <- response_l - mean_yl
  W_yu <- response_u - mean_yu
  W_yr <- response_r - mean_yr
  
  for(i in 1:np){
    for(j in 1:length(response)){
      W_x[[i]][j,] = W_x[[i]][j,] - mean_x[[i]]
      W_xl[[i]][j,] = W_xl[[i]][j,] - mean_xl[[i]]
      W_xu[[i]][j,] = W_xu[[i]][j,] - mean_xu[[i]]
      W_xr[[i]][j,] = W_xr[[i]][j,] - mean_xr[[i]]
    }
  }
  
  for(i in 1:np){
    for(j in 1:nrow(predictor_test[[1]])){
      W_x_test[[i]][j,] = W_x_test[[i]][j,] - mean_x_test[[i]]
      W_xl_test[[i]][j,] = W_xl_test[[i]][j,] - mean_xl_test[[i]]
      W_xu_test[[i]][j,] = W_xu_test[[i]][j,] - mean_xu_test[[i]]
      W_xr_test[[i]][j,] = W_xr_test[[i]][j,] - mean_xr_test[[i]]
    }
  }
  
  
  # Matrices for the regressions
  Reg_mat = vector("list",)
  Reg_mat_l = vector("list",)
  Reg_mat_u = vector("list",)
  Reg_mat_r = vector("list",)
  
  Reg_mat_test = vector("list",)
  Reg_mat_test_l = vector("list",)
  Reg_mat_test_u = vector("list",)
  Reg_mat_test_r = vector("list",)
  
  for(i in 1:np){
    Reg_mat[[i]] = W_x[[i]] %*% Inner_prod[[i]]
    Reg_mat_l[[i]] = W_xl[[i]] %*% Inner_prod[[i]]
    Reg_mat_u[[i]] = W_xu[[i]] %*% Inner_prod[[i]]
    Reg_mat_r[[i]] = W_xr[[i]] %*% Inner_prod[[i]]
    
    Reg_mat_test[[i]] = W_x_test[[i]] %*% Inner_prod[[i]]
    Reg_mat_test_l[[i]] = W_xl_test[[i]] %*% Inner_prod[[i]]
    Reg_mat_test_u[[i]] = W_xu_test[[i]] %*% Inner_prod[[i]]
    Reg_mat_test_r[[i]] = W_xr_test[[i]] %*% Inner_prod[[i]]
  }
  
  Reg_mat = cbind(1, do.call(cbind, Reg_mat))
  Reg_mat_l = cbind(1, do.call(cbind, Reg_mat_l))
  Reg_mat_u = cbind(1, do.call(cbind, Reg_mat_u))
  Reg_mat_r = cbind(1, do.call(cbind, Reg_mat_r))
  Reg_mat_bcrm = cbind(1, cbind(Reg_mat[, -1], Reg_mat_r[, -1]))
  
  Reg_mat_test = cbind(1, do.call(cbind, Reg_mat_test))
  Reg_mat_test_l = cbind(1, do.call(cbind, Reg_mat_test_l))
  Reg_mat_test_u = cbind(1, do.call(cbind, Reg_mat_test_u))
  Reg_mat_test_r = cbind(1, do.call(cbind, Reg_mat_test_r))
  Reg_mat_test_bcrm = cbind(1, cbind(Reg_mat_test[, -1], Reg_mat_test_r[, -1]))
  
  
  # Model estimation
  coeff_c = ginv(t(Reg_mat)%*%Reg_mat) %*% t(Reg_mat)%*%W_y
  coeff_l = ginv(t(Reg_mat_l)%*%Reg_mat_l) %*% t(Reg_mat)%*%W_yl
  coeff_u = ginv(t(Reg_mat_u)%*%Reg_mat_u) %*% t(Reg_mat_u)%*%W_yu
  coeff_r = ginv(t(Reg_mat_r)%*%Reg_mat_r) %*% t(Reg_mat_r)%*%W_yr
  coeff_bcrm_c = ginv(t(Reg_mat_bcrm)%*%Reg_mat_bcrm) %*% t(Reg_mat_bcrm)%*%W_y
  coeff_bcrm_r = ginv(t(Reg_mat_bcrm)%*%Reg_mat_bcrm) %*% t(Reg_mat_bcrm)%*%W_yr
  
  pred_l = as.vector(Reg_mat_test_l %*% coeff_l) + mean_yl
  pred_u = as.vector(Reg_mat_test_u %*% coeff_u) + mean_yu
  
  pred_cm_l_in = as.vector(Reg_mat_test_l %*% coeff_c) + mean_yl
  pred_cm_u_in = as.vector(Reg_mat_test_u %*% coeff_c) + mean_yu
  
  pred_bcrm_c = as.vector(Reg_mat_test_bcrm %*% coeff_bcrm_c) + mean_y
  pred_bcrm_r = as.vector(Reg_mat_test_bcrm %*% coeff_bcrm_r) + mean_yr
  pred_bcrm_l = pred_bcrm_c - pred_bcrm_r
  pred_bcrm_u = pred_bcrm_c + pred_bcrm_r
  
  pred_cm_l = apply(cbind(pred_cm_l_in, pred_cm_u_in), 1, min)
  pred_cm_u = apply(cbind(pred_cm_l_in, pred_cm_u_in), 1, max)
  
  return(list("flm_l" = pred_l, "flm_u" = pred_u, 
              "cm_l" = pred_cm_l, "cm_u" = pred_cm_u, 
              "bcrm_l" = pred_bcrm_l, "bcrm_u" = pred_bcrm_u))
}



