setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("scatterplot3d")
library(np) # non parametric library
library(rgl)
library(plot.matrix)
library(MASS)

source("library_smooth.R")

# resolution: how many pixels per:
K1 <- 50 # moneyness variable
K2 <- K1 # days_to_expiration_variable

# bounds for the variables
moneyness_min <- -0.5 #0.5
moneyness_max <-  0.5 #1.5
days_to_expiration_min <- 14
days_to_expiration_max <- 365


################
## load data ###
################
surfaces_all <- read.csv("calls_clean_all_more_columns.csv")

# log transform
surfaces_all$IV <- log(surfaces_all$IV) # log transformation here !!! XXX !!!
surfaces_all$moneyness <- log(surfaces_all$moneyness)

# create a subset
surfaces_all <- subset(surfaces_all, days_to_expiration <= days_to_expiration_max)
surfaces_count <- length(unique(surfaces_all$surface_id))

#############################
## round to a common grid ###
#############################

# data matrix
X <- array(NA, dim=c(surfaces_count, K1, K2))
X_counts <- array(0, dim=c(surfaces_count, K1, K2))

surface_date <- numeric(surfaces_count)
surface_DataDate <- numeric(surfaces_count)
surface_symbol <- character(surfaces_count)

# this is for the mean smoother
X_for_mean <- matrix(NA, ncol=2, nrow=nrow(surfaces_all))
Y_for_mean <- rep(NA, nrow(surfaces_all))

# go through all observations
for (obs_i in 1:nrow(surfaces_all)){
  
  # extract data from the current observation
  surface_id <- surfaces_all$surface_id[obs_i]
  moneyness <- surfaces_all$moneyness[obs_i]
  days_to_expiration <- surfaces_all$days_to_expiration[obs_i]
  IV <- surfaces_all$IV[obs_i]
  surface_DataDate[surface_id] <- surfaces_all$DataDate[obs_i]
  surface_date[surface_id] <- difftime( as.Date(surfaces_all$DataDate[obs_i], format="%m/%d/%Y"), as.Date("01/01/2003", format="%m/%d/%Y") )
  surface_symbol[surface_id] <- surfaces_all$UnderlyingSymbol[obs_i]
  
  if ((moneyness > moneyness_min) && (moneyness < moneyness_max) && (days_to_expiration > days_to_expiration_min) && (days_to_expiration < days_to_expiration_max)){
    
    # round moneyness
    moneyness_01 <- (moneyness - moneyness_min)/(moneyness_max - moneyness_min) # scale to [0,1]
    moneyness_grid <- round(moneyness_01 * K1 + 0.5)
    if (moneyness_grid == 0){ moneyness_grid <- 1 }
    if (moneyness_grid > K1){ moneyness_grid <- K1 }
    
    # round days_to_expiration
    days_to_expiration_01 <- (days_to_expiration - days_to_expiration_min)/(days_to_expiration_max - days_to_expiration_min) # scale to [0,1]
    days_to_expiration_grid <- round(days_to_expiration_01 * K2 + 0.5)
    if (days_to_expiration_grid == 0){ days_to_expirationgrid <- 1 }
    if (days_to_expiration_grid > K2){ days_to_expiration_grid <- K2 }
    
    # write in the data
    count_now <- X_counts[ surface_id, moneyness_grid, days_to_expiration_grid]
    if (length(count_now) > 1){print(obs_i)}
    if (count_now == 0){
      X[ surface_id, moneyness_grid, days_to_expiration_grid] <- IV
      X_counts[ surface_id, moneyness_grid, days_to_expiration_grid] <- count_now + 1
    } else {
      X[ surface_id, moneyness_grid, days_to_expiration_grid] <- (X[ surface_id, moneyness_grid, days_to_expiration_grid]*count_now + IV)/(count_now+1)
      X_counts[ surface_id, moneyness_grid, days_to_expiration_grid] <- count_now + 1
    }
    
    # for the mean smoother
    X_for_mean[obs_i,] <- c( moneyness_grid, days_to_expiration_grid)
    Y_for_mean[obs_i] <- IV
 
  } 
}

##################
### Estimation ###
##################

# mean
xy_grid <- as.matrix(expand.grid(1:K1, 1:K2))
surface_mean_res <- npreg(bws = c(K1/10,K2/10), xdat=X_for_mean, ydat=Y_for_mean, exdat = xy_grid, regtype = "ll")
surface_mean <- matrix( surface_mean_res$mean, nrow=K1, ncol=K2)

X_centred <- array(NA, dim = dim(X))
for (surface_id in 1:surfaces_count){
  X_centred[surface_id,,] <- X[surface_id,,] - surface_mean
}

# covariance structure and noise level
make_matrix_posdefinite <- function(m){
  res <- eigen( (m+t(m))/2 )
  res$values[res$values < 0] <- 0
  return(res$vectors %*% diag(res$values) %*% t(res$vectors))
}

Res <- estimate_sparse(X_centred, WN_flag=T)

# make them pos def
Res_A_posdef <- make_matrix_posdefinite(Res$A)
Res_B_posdef <- make_matrix_posdefinite(Res$B)

# eigenfunctions
eigen_A <- eigen(Res_A_posdef)
eigen_B <- eigen(Res_B_posdef)

# for prediction function (currently slow implementation using full covariance tensor)
prior_mu <- c(surface_mean)
prior_var <- kronecker(Res_B_posdef, Res_A_posdef)

# function for prediction
predict_surface_var <- function(surface_condition, prior_mu, prior_var, sigma2){
  
  # vectorize observed data
  obs_data_which <- !is.na(c(surface_condition))
  obs_data <- c(surface_condition)[obs_data_which]
  obs_data_dim <- length(obs_data)
  
  # create the censor matrix: mapping from the latent vectorized surface into the vector of observations
  K1K2 <- dim(surface_condition)
  K1 <- K1K2[1]
  K2 <- K1K2[2]
  censor_matrix <- matrix( 0, nrow=obs_data_dim, ncol=K1*K2 )
  obs_i <- 1
  for ( i in 1:(K1*K2) ){
    if (obs_data_which[i]){
      censor_matrix[obs_i, i] <- 1
      obs_i <- obs_i + 1
    }
  }
  
  ## prediction formula
  posterior <- list()
  condition_var_inv <- solve( censor_matrix %*% prior_var %*% t(censor_matrix) + sigma2 * diag(obs_data_dim) )
  posterior$mu_flattened <- prior_mu + prior_var %*% t(censor_matrix) %*% condition_var_inv %*% ( obs_data - censor_matrix %*% prior_mu ) 
  posterior$mu <- matrix(posterior$mu_flattened, nrow=K1, ncol=K2)
  posterior$var <- prior_var - prior_var %*% t(censor_matrix) %*% condition_var_inv %*% censor_matrix %*% prior_var
  posterior$diag_var <- diag(posterior$var)
  posterior$cor <- posterior$var / outer(sqrt(posterior$diag_var),sqrt(posterior$diag_var))
  
  
  ## pointwise confidence bands
  posterior$pointwise_band_lower_flattened <- posterior$mu - 1.96 * sqrt(posterior$diag_var)
  posterior$pointwise_band_upper_flattened <- posterior$mu + 1.96 * sqrt(posterior$diag_var)
  posterior$pointwise_band_lower <- matrix(posterior$pointwise_band_lower_flattened, nrow=K1, ncol=K2)
  posterior$pointwise_band_upper <- matrix(posterior$pointwise_band_upper_flattened, nrow=K1, ncol=K2)
  
  ## simultaneous confidence bands
  # simulate the quantile z
  sample <- mvrnorm( n=1000, mu=numeric(K1*K2), Sigma=posterior$cor )
  maxima <- apply(abs(sample), 1, max)
  z <- quantile(maxima, 0.95)
  
  
  # construct the bands
  posterior$simult_band_lower_flattened <- posterior$mu - z * sqrt(posterior$diag_var)
  posterior$simult_band_upper_flattened <- posterior$mu + z * sqrt(posterior$diag_var)
  posterior$simult_band_lower <- matrix(posterior$simult_band_lower_flattened, nrow=K1, ncol=K2)
  posterior$simult_band_upper <- matrix(posterior$simult_band_upper_flattened, nrow=K1, ncol=K2)
  
  return(posterior)
}

# predict Dell
id <- 37
surface_full <- X[id,,]

## predict the surface
posterior <- predict_surface_var(surface_full, prior_mu, prior_var, Res$sigma2)

### save data for plotting later
posterior_new <- list(Mean = posterior$mu, Upper = posterior$simult_band_upper, Lower=posterior$simult_band_lower)
TRout <- list(surface_mean=surface_mean, Dell=X[37,,], Qualcomm=X[182,,], Res=Res,
                  EIG1=eigen_A, EIG2=eigen_B, posterior=posterior_new)
save(TRout, file="TRout.RData")

