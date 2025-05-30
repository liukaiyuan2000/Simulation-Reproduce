rm(list = ls())
setwd("D:/桌面/QR_Sca&Fun")
source("qr_fof_lqSolve.R")
source("lr_fof.R")
library(lubridate)
## data set preprocessing ----
orgin <- read.csv('D:/桌面/alldata_skab/data.csv', header = T)
Hour <- hour(orgin$datetime)
orgin <- cbind(orgin, Hour)

temp <- list()
time_dif <- 1:31
a = 0
for(i in time_dif){
  a <- a + 1
  temp[[a]] <- orgin[orgin$Day == i, ]
}

N = sum(!duplicated(orgin$market)) # 905 markets
Idx <- orgin$market[!duplicated(orgin$market)]
temp2 <- matrix(0, N, length(time_dif))
min_price <- data.frame(temp2, row.names = Idx)
max_price <- data.frame(temp2, row.names = Idx)
min_modal <- data.frame(temp2, row.names = Idx)
max_modal <- data.frame(temp2, row.names = Idx)


for(j in Idx){
  for(i in time_dif){
    min_price[j, i] <- min(temp[[i]][temp[[i]]$market == j, ]$min_price)
    max_price[j, i] <- max(temp[[i]][temp[[i]]$market == j, ]$max_price)
    min_modal[j, i] <- min(temp[[i]][temp[[i]]$market == j, ]$modal_price)
    max_modal[j, i] <- max(temp[[i]][temp[[i]]$market == j, ]$modal_price)
  }
}

Idx2 <- (apply(min_price < Inf, 1, sum) == 31)
min_price <- min_price[Idx2, ]
max_price <- max_price[Idx2, ]
min_modal <- min_modal[Idx2, ]
max_modal <- max_modal[Idx2, ]
Idx3 <- (apply((max_price - min_price) == 0, 1, sum) == 0)
min_price <- min_price[Idx3, ]
max_price <- max_price[Idx3, ]
min_modal <- min_modal[Idx3, ]
max_modal <- max_modal[Idx3, ]
#####

## fundemental set ######
# Number of basis functions for predictors
nbf_vec_predictors = rep(20, 2)
# Number of simulations
nsim = 50
# Sample sizes of training and test sets
N = nrow(min_price)
n_train = floor(0.8 * N)
n_test = N - n_train
# Quantile para
tau = seq(0.1, 0.9, by=0.1)
# data.frame rownames
rowname <- c()
for(i in 1:length(tau)){
  rowname[i] <- paste('tau = ', tau[i], sep = "")
}
rowname <- c(rowname, 'LR')
#####

## generate data set ----
#　response - Purchase price
y_l <- apply(min_modal, 1, min)
y_u <- apply(max_modal, 1, max)
y_c <- (y_u + y_l) / 2
y_r <- (y_u - y_l) / 2
# predictor - Selling price
x_l <- min_price
x_u <- max_price
x_c <- (x_u + x_l) / 2
x_r <- (x_u - x_l) / 2
#####

###################### Simulations start #############################
mae_qr_l <- matrix(0, nsim, length(tau))
mae_qr_u <- matrix(0, nsim, length(tau))
mae_lr_l <- matrix(0, nsim, length(tau))
mae_lr_u <- matrix(0, nsim, length(tau))

for(j in 1:length(tau)){
  for(sim in 1:nsim){
    cat('tau = ', tau[j], '; 第 ', sim, '/', nsim,' 次模拟\r', sep = "")
    test_samp = sample(1:(n_train + n_test), n_train, replace = F)
    
    ## generate data
    response_test_l = y_l[-test_samp]
    response_test_u = y_u[-test_samp]
    response = y_c[test_samp]
    response_l = y_l[test_samp]
    response_u = y_u[test_samp]
    response_r = y_r[test_samp]
    
    predictor_test = list(x_c[-test_samp, ])
    predictor_test_l = list(x_l[-test_samp, ])
    predictor_test_u = list(x_u[-test_samp, ])
    predictor_test_r = list(x_r[-test_samp, ])
    predictor = list(x_c[test_samp, ])
    predictor_l = list(x_l[test_samp, ])
    predictor_u = list(x_u[test_samp, ])
    predictor_r = list(x_r[test_samp, ])
    
    ## estimation
    qr_run = qr_fof(
      response, response_l, response_u, response_r, 
      predictor, predictor_l, predictor_u, predictor_r, 
      predictor_test, predictor_test_l, predictor_test_u, predictor_test_r,
      nbf_vec_predictors, num_t = 31, np = 1, n_train, tau = tau[j]
    )
    lr_run = lr_fof(
      response, response_l, response_u, response_r, 
      predictor, predictor_l, predictor_u, predictor_r, 
      predictor_test, predictor_test_l, predictor_test_u, predictor_test_r,
      nbf_vec_predictors, num_t = 31, np = 1, n_train
    )
    
    ## mae of quantile regression and linear regression
    mae_qr_l[sim, j] <- mean(abs((qr_run$bcrm_l - response_test_l)))
    mae_qr_u[sim, j] <- mean(abs((qr_run$bcrm_u - response_test_u)))
    mae_lr_l[sim, j] <- mean(abs((lr_run$bcrm_l - response_test_l)))
    mae_lr_u[sim, j] <- mean(abs((lr_run$bcrm_u - response_test_u)))
  }
  cat('\ndone~\n')
}
############################# End ###########################

## result output ----
result_l <- data.frame(Mean.l = round(c(apply(mae_qr_l, 2, mean), 
                                        mean(mae_lr_l)), 4), 
                       Sd.l = round(c(apply(mae_qr_l, 2, sd), 
                                      sd(mae_lr_l)), 4), 
                       row.names = rowname)
result_u <- data.frame(Mean.u = round(c(apply(mae_qr_u, 2, mean), 
                                        mean(mae_lr_u)), 4), 
                       Sd.u = round(c(apply(mae_qr_u, 2, sd), 
                                      sd(mae_lr_u)), 4), 
                       row.names = rowname)
cbind(result_l, result_u)
## Latex code
# stargazer(cbind(result_l, result_u), summary = F)
#####

















