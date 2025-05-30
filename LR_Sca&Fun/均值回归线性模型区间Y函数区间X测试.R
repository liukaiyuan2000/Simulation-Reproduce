rm(list=ls())
setwd("D:/桌面/FLM_interval_valued_data-master")
# source("auxiliary_functions.R")

# Number of basis functions for predictors
nbf_vec_predictors = rep(10,2)
# Number of basis functions for response
nbf_response = 10
# Number of Monte Carlo simulations
B = 100
# Number of simulations
nsim = 100
# Number of functions in the training and test samples
n_train = 40
n_test = 8
# Number of predictors
np = 2

mse_flm_l = numeric()
mse_flm_u = numeric()
mse_cm_l = numeric()
mse_cm_u = numeric()
mse_bcrm_l = numeric()
mse_bcrm_u = numeric()


# Load datasets
max_t = load("max_temperature.RData")
min_t = load("min_temperature.RData")
mean_t = load("mean_temperature.RData")
hr_t = load("hr_temperature.RData")

max_h = load("max_humidity.RData")
min_h = load("min_humidity.RData")
mean_h = load("mean_humidity.RData")
hr_h = load("hr_humidity.RData")

max_e = load("max_evaporation.RData")
min_e = load("min_evaporation.RData")
mean_e = load("mean_evaporation.RData")
hr_e = load("hr_evaporation.RData")

# List
predictor_orig = list(mean_temperature, mean_humidity)
predictor_l_orig = list(min_temperature, min_humidity)
predictor_u_orig = list(max_temperature, max_humidity)
predictor_r_orig = list(hr_temperature, hr_humidity)

y <- apply(mean_evaporation, 1, mean)
y_l <- apply(min_evaporation, 1, mean)
y_u <- apply(max_evaporation, 1, mean)
y_r <- apply(hr_evaporation, 1, mean)

##### Simulation ----
for(sim in 1:nsim){
  cat('第 ', sim, '/', nsim,' 次模拟\r', sep = "")
  test_samp = sample(1:(n_train+n_test), n_train, replace = FALSE)
  
  response = as.vector(y[test_samp])
  response_l = as.vector(y_l[test_samp])
  response_u = as.vector(y_u[test_samp])
  response_r = as.vector(y_r[test_samp])
  response_test_l = as.vector(y_l[-(test_samp)])
  response_test_u = as.vector(y_u[-(test_samp)])
  
  predictor = list()
  predictor_l = list()
  predictor_u = list()
  predictor_r = list()
  predictor_test = list()
  predictor_test_l = list()
  predictor_test_u = list()
  predictor_test_r = list()
  
  for(i in 1:np){
    predictor[[i]] = predictor_orig[[i]][test_samp,]
    predictor_l[[i]] = predictor_l_orig[[i]][test_samp,]
    predictor_u[[i]] = predictor_u_orig[[i]][test_samp,]
    predictor_r[[i]] = predictor_r_orig[[i]][test_samp,]
    predictor_test[[i]] = predictor_orig[[i]][-(test_samp),]
    predictor_test_l[[i]] = predictor_l_orig[[i]][-(test_samp),]
    predictor_test_u[[i]] = predictor_u_orig[[i]][-(test_samp),]
    predictor_test_r[[i]] = predictor_r_orig[[i]][-(test_samp),]
  }
  
  m_run = iv_fof(
    response, response_l, response_u, response_r, 
    predictor, predictor_l, predictor_u, predictor_r, 
    predictor_test, predictor_test_l, predictor_test_u, predictor_test_r, 
    nbf_vec_predictors, num_t = 12, np, n_train)
  
  mse_flm_l[sim] = mean((m_run$flm_l - response_test_l)^2)
  mse_flm_u[sim] = mean((m_run$flm_u - response_test_u)^2)
  mse_cm_l[sim] = mean((m_run$cm_l - response_test_l)^2)
  mse_cm_u[sim] = mean((m_run$cm_u - response_test_u)^2)
  mse_bcrm_l[sim] = mean((m_run$bcrm_l - response_test_l)^2)
  mse_bcrm_u[sim] = mean((m_run$bcrm_u - response_test_u)^2)
}
#####


##### Plot ----
mse_l <- as.data.frame(cbind(mse_flm_l, mse_cm_l, mse_bcrm_l))
mse_u <- as.data.frame(cbind(mse_flm_u, mse_cm_u, mse_bcrm_u))
names(mse_l) <- c('FLM', 'CM', 'BCRM')
names(mse_u) <- c('FLM', 'CM', 'BCRM')
# AMSE
par(mfrow = c(1,2))
boxplot(mse_l, ylab = 'AMSE')
boxplot(mse_u, ylab = 'AMSE')
par(mfrow = c(1,1))


# Temperature
plot(0, 0, xlab = 's', ylab = '摄氏度', main = 'Temperature', 
     xlim = c(1, 12), ylim = c(0, 50), type = 'n')
for(i in 1:48){
  lines(1:12, max_temperature[i, ])
  lines(1:12, min_temperature[i, ], col = 2)
}
# Humidity
plot(0, 0, xlab = 's', ylab = '%', main = 'Humidity', 
     xlim = c(1, 12), ylim = c(0, 100), type = 'n')
for(i in 1:48){
  lines(1:12, max_humidity[i, ])
  lines(1:12, min_humidity[i, ], col = 2)
}
# Evaporation
plot(0, 0, xlab = 't', ylab = 'mm', main = 'Evaporation', 
     xlim = c(1, 12), ylim = c(0, 80), type = 'n')
for(i in 1:48){
  lines(1:12, max_evaporation[i, ])
  lines(1:12, min_evaporation[i, ], col = 2)
}

#####






