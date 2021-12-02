
# function for prediction
predict_surface <- function(surface_condition, prior_mu, prior_var, sigma2){
  
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
  
  # prediction formula
  posterior_mu <- prior_mu + prior_var %*% t(censor_matrix) %*% solve( censor_matrix %*% prior_var %*% t(censor_matrix) +
                                                                         sigma2 * diag(obs_data_dim) ) %*% ( obs_data - censor_matrix %*% prior_mu ) 
  
  return(matrix(posterior_mu, nrow=K1, ncol=K2))
  
  
}


###########################################################################################################
## parallel setting
# 
# library(doParallel)
# if (!exists("global_setting_cores")){
#   global_setting_cores <- 4
#   registerDoParallel(global_setting_cores)
# }

###########################################################################################################

# library("scatterplot3d")
library(np) # non parametric library
library(foreach)
# library(rgl)
# library(plot.matrix)



# resolution: how many pixels per:
K1 <- 20 # moneyness variable
K2 <- K1 # days_to_expiration_variable

# bounds for the variables
moneyness_min <- -0.5
moneyness_max <-  0.5
days_to_expiration_min <- 14
days_to_expiration_max <- 365


###############################
## load data

# setwd("C:/C-epfl results/Sparse separable covariances/v11 - k-fold comparison - separable vs presmooth/")
surfaces_all <- read.csv("calls_clean_all.csv")

# log transform
surfaces_all$IV <- log(surfaces_all$IV) # log transformation here !!! XXX !!!
surfaces_all$moneyness <- log(surfaces_all$moneyness) # log transformation here !!! XXX !!!

# create a subset
surfaces_all <- subset(surfaces_all, days_to_expiration <= days_to_expiration_max)
surfaces_count <- length(unique(surfaces_all$surface_id))

##########################################
## K-Fold cross validation
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
k_fold_k <- 10  # here K=10 !!!XXX!!!
k_folds <- chunk2(1:surfaces_count, k_fold_k)





# out <- foreach(k_fold_i = 1:k_fold_k, .packages=c('np','rgl') ) %dopar% {
out <- foreach(k_fold_i = 1:k_fold_k, .packages=c('np') ) %do% {
# out <- for(k_fold_i in 1:k_fold_k){
  print(k_fold_i)
  
  
  ## split train and test
  surfaces_test_ids <- k_folds[[k_fold_i]]
  surfaces_train_ids <- (1:surfaces_count)[is.na(pmatch((1:surfaces_count),surfaces_test_ids))]
  
  # prepare out structure
  out_ratios_abcde <- mm<-matrix(list(), 1, 5)
  out_rmse_4D_abcde <- mm<-matrix(list(), 1, 5)
  out_rmse_our_abcde <- mm<-matrix(list(), 1, 5)
  out_rmse_smoother_abcde <- mm<-matrix(list(), 1, 5)
  for (abcde in 1:5){
    out_ratios_abcde[[abcde]] <- rep(NA, length(surfaces_test_ids))
    out_rmse_4D_abcde[[abcde]] <- rep(NA, length(surfaces_count))
    out_rmse_our_abcde[[abcde]] <- rep(NA, length(surfaces_count))
    out_rmse_smoother_abcde[[abcde]] <- rep(NA, length(surfaces_count))
  }
  
  ################################################################################
  ## round to a common grid
  
  # data matrix
  X <- array(NA, dim=c(surfaces_count, K1, K2))
  X_test <- array(NA, dim=c(surfaces_count, K1, K2))
  X_counts <- array(0, dim=c(surfaces_count, K1, K2))
  X_test_counts <- array(0, dim=c(surfaces_count, K1, K2))
  
  # this is for mean smoother
  X_for_mean <- matrix(NA, ncol=2, nrow=nrow(surfaces_all))
  Y_for_mean <- rep(NA, nrow(surfaces_all))
  
  # go through all observations
  for (obs_i in 1:nrow(surfaces_all)){
    
    # extract data from the current observation
    surface_id <- surfaces_all$surface_id[obs_i]
    moneyness <- surfaces_all$moneyness[obs_i]
    days_to_expiration <- surfaces_all$days_to_expiration[obs_i]
    IV <- surfaces_all$IV[obs_i]
    
    # round moneyness
    moneyness_01 <- (moneyness - moneyness_min)/(moneyness_max - moneyness_min) # scale to [0,1]
    if ((moneyness_01 >= 0) & (moneyness_01 <= 1)){
      
      moneyness_grid <- round(moneyness_01 * K1 + 0.5)
      if (moneyness_grid == 0){ moneyness_grid <- 1 }
      if (moneyness_grid > K1){ moneyness_grid <- K1 }
      
      # round days_to_expiration
      days_to_expiration_01 <- (days_to_expiration - days_to_expiration_min)/(days_to_expiration_max - days_to_expiration_min) # scale to [0,1]
      days_to_expiration_grid <- round(days_to_expiration_01 * K2 + 0.5)
      if (days_to_expiration_grid == 0){ days_to_expirationgrid <- 1 }
      if (days_to_expiration_grid > K2){ days_to_expiration_grid <- K2 }
      
      if (surface_id %in% surfaces_train_ids){
        # write in the data
        count_now <- X_counts[ surface_id, moneyness_grid, days_to_expiration_grid]
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
        
        
      } else { # test partition
        count_now_test <- X_test_counts[ surface_id, moneyness_grid, days_to_expiration_grid]
        if (count_now_test == 0){
          X_test[ surface_id, moneyness_grid, days_to_expiration_grid] <- IV
          X_test_counts[ surface_id, moneyness_grid, days_to_expiration_grid] <- count_now_test + 1
        } else {
          X_test[ surface_id, moneyness_grid, days_to_expiration_grid] <- (X_test[ surface_id, moneyness_grid, days_to_expiration_grid]*count_now_test + IV)/(count_now_test+1)
          X_test_counts[ surface_id, moneyness_grid, days_to_expiration_grid] <- count_now_test + 1
        }
        
      }
    }
  }
  
  
  ############################################################################################################################
  ## mean estimation 
  
  
  surface_mean_naive <- apply(X, c(2,3), mean, na.rm = TRUE)
  obs_counts <- apply(X_counts, c(2,3), sum, na.rm = TRUE)
  
  ## mean estimator - smoother
  # surface_mean <- npregbw(xdat = X_for_mean, ydat = Y_for_mean, regtype = "ll")
  xy_grid <- as.matrix(expand.grid(1:K1, 1:K2))
  surface_mean_res <- npreg(bws = c(K1/10,K2/10), xdat=X_for_mean, ydat=Y_for_mean, exdat = xy_grid, regtype = "ll")
  surface_mean <- matrix( surface_mean_res$mean, nrow=K1, ncol=K2)
  
  # surface_mean <- surface_mean_naive
  
  # center the surfacase in order to use Masak's code
  X_centred <- array(NA, dim = dim(X))
  for (surface_id in 1:surfaces_count){
    X_centred[surface_id,,] <- X[surface_id,,] - surface_mean
  }
  
  
  ############################################################################################################################
  ## function to adjust matrix to be pos-definite
  
  make_matrix_posdefinite <- function(m){
    res <- eigen( (m+t(m))/2 )
    res$values[res$values < 0] <- 0
    return(res$vectors %*% diag(res$values) %*% t(res$vectors))
  }
  
  
  
  
  ############################################################################################################################
  ## separable smoother
  # Res <- estimate_sparse(X_centred, WN_flag=T,bwidths=c(6,10))
  # Res <- estimate_sparse(X_centred, WN_flag=T,bwidths=c(8,8))
  
  source("library_smooth.R")
  start_time <- Sys.time()
  # Res <- estimate_sparse(X_centred, WN_flag=T,bwidths=c(K1/10*5,K2/10*5))
  Res <- estimate_sparse(X_centred, WN_flag=T)
  end_time <- Sys.time()
  timing_sparse <- end_time - start_time
  
  # make them pos def
  Res_A_posdef <- make_matrix_posdefinite(Res$A)
  Res_B_posdef <- make_matrix_posdefinite(Res$B)
  
  
  
  ############################################################################################################################
  ## 4D smoother
  
  start_time <- Sys.time()
  Res_4D <- smooth4D(X_centred,bw1=Res$A_bandwidth,bw2=Res$B_bandwidth,WN_flag=T)
  end_time <- Sys.time()
  timing_4D <- end_time - start_time
  # Res$C has dimension .... dim=c(K1,K2,K1,K2)
  
  
  
  
  #######################################################################################################################
  ## prepare prior distribution
  prior_mu <- c(surface_mean)
  prior_var <- matrix(NA, nrow=K1*K2, ncol=K1*K2)
  prior_var_4D <- matrix(NA, nrow=K1*K2, ncol=K1*K2)
  for (i in 1:K1){
    for (j in 1:K2){
      for (k in 1:K1){
        for (l in 1:K2){
          prior_var[i+(j-1)*K1, k+(l-1)*K1] <- Res_A_posdef[i,k] * Res_B_posdef[j,l]
          prior_var_4D[i+(j-1)*K1, k+(l-1)*K1] <- Res_4D$C[i,j,k,l]
        }
      }
    }
  }
  
  prior_var_4D <- make_matrix_posdefinite(prior_var_4D)
  
  
  ## print timing
  print(timing_4D)
  print(timing_sparse)
  
  
  
  #######################################################################################################################
  ## prediction - other strategies
  
  
  
  for (abcde in 2:5){
    
    
    # prior distribution
    prior_mu <- c(surface_mean)
    prior_var <- kronecker(Res_B_posdef, Res_A_posdef)
    
    
    # prepare all ids
    ids_all <- surfaces_test_ids
    
    
    # ids_all <- c(37)
    mse_our_all <- rmse_our_all <- rep(NA, length(ids_all))
    mse_4D_all <- rmse_4D_all <- rep(NA, length(ids_all))
    mse_smoother_all <- rmse_smoother_all <- rep(NA, length(ids_all))
    
    for (id_i in 1:length(ids_all)){
      id <- ids_all[id_i]
      # the surface to predict
      # plot(X_test[id,,], xlab="time to expiration", ylab="moneyness")
      surface_full <- X_test[id,,]
      
      
      # (b) keep out-the-money (i.e. those with moneyness > 1)
      if (abcde == 2){
        surface_condition <- surface_test <- surface_full
        surface_condition[1:ceiling(K1/2),] <- NA # delete observation with moneyness <= 1
      }
      
      # (c) keep in-the-money (i.e. those with moneyness < 1)
      if (abcde == 3){
        surface_condition <- surface_test <- surface_full
        surface_condition[floor(K1/2):K1,] <- NA # delete observation with moneyness >= 1
      }
      
      # (d) keep only short expirations
      if (abcde == 4){
        surface_condition <- surface_test <- surface_full
        surface_condition[, floor(K2/2):K2] <- NA # delete observation with expiration in the second half of the interval
      }
      
      # (e) keep long expirations
      if (abcde == 5){
        surface_condition <- surface_test <- surface_full
        surface_condition[, 1:ceiling(K1/2)] <- NA # delete observation with expiration in the first half of the interval
      }
      
      
      ## the test partition
      surface_test[!is.na(surface_condition)] <- NA # delete obs. that are present in the condition
      # plot(surface_condition)
      # plot(surface_test)
      
      #############################################################################################################
      ## smooth first approach
      
      surface_condition_exp <- exp(surface_condition)
      # plot(surface_condition_exp)
      surface_condition_num_obs <- sum(!is.na(surface_condition_exp))
      
      # interpolate
      X_npreg <- matrix(NA, nrow=surface_condition_num_obs, ncol=2)
      Y_npreg <- rep(NA, surface_condition_num_obs)
      obs_i <- 1
      for (i in 1:K1){
        for (j in 1:K2){
          if (!is.na(surface_condition_exp[i,j])){
            X_npreg[obs_i,1] <- i
            X_npreg[obs_i,2] <- j
            Y_npreg[obs_i] <- surface_condition_exp[i,j]
            obs_i <- obs_i + 1
          }
        }
      }
      xy_grid <- as.matrix(expand.grid(1:K1, 1:K2))
      
      # continue only if there are at least 2 different expiration dates, and if at least 2 strikes
      if ((length(unique(X_npreg[,1])) >= 2) && (length(unique(X_npreg[,2])) >= 2)){
        
        
        kre <- npreg(bws=c(2,3), xdat=X_npreg, ydat=Y_npreg, exdat=xy_grid, regtype = "ll")
        surface_smoothed <- matrix(kre$mean, nrow=K1, ncol=K2)
        
        # RMSE denominator
        relativity <- sqrt(var(c(surface_full), na.rm = TRUE))
        
        # evaluate MSE for this surface
        test_mask <- !is.na(surface_test)
        if (is.na(mse_smoother_all[id_i])){mse_smoother_all[id_i] <- 0}
        if (is.na(rmse_smoother_all[id_i])){rmse_smoother_all[id_i] <- 0}
        mse_smoother_all[id_i] <- mse_smoother_all[id_i] + sqrt(mean((exp(surface_test[test_mask]) - surface_smoothed[test_mask])^2))
        rmse_smoother_all[id_i] <- rmse_smoother_all[id_i] +
          sqrt(mean((exp(surface_test[test_mask]) - surface_smoothed[test_mask])^2)) / relativity
        
        
        #############################################################################################################
        ## predict the surface using our method
        surface_posterior <- predict_surface(surface_condition, prior_mu, prior_var, Res$sigma2)
        surface_posterior_4D <- predict_surface(surface_condition, prior_mu, prior_var_4D, Res_4D$sigma2)
        
        # separable error
        test_mask <- !is.na(surface_test)
        if (is.na(mse_our_all[id_i])){mse_our_all[id_i] <- 0}
        if (is.na(rmse_our_all[id_i])){rmse_our_all[id_i] <- 0}
        mse_our_all[id_i] <- mse_our_all[id_i] + sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior[test_mask]))^2 ))
        rmse_our_all[id_i] <- rmse_our_all[id_i] +
          sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior[test_mask]))^2 )) / relativity
        
        # 4D error
        test_mask <- !is.na(surface_test)
        if (is.na(mse_4D_all[id_i])){mse_4D_all[id_i] <- 0}
        if (is.na(rmse_4D_all[id_i])){rmse_4D_all[id_i] <- 0}
        mse_4D_all[id_i] <- mse_4D_all[id_i] + sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior_4D[test_mask]))^2 ))
        rmse_4D_all[id_i] <- rmse_4D_all[id_i] +
          sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior_4D[test_mask]))^2 )) / relativity
        
        
      }
      
    }
    
    out_rmse_our_abcde[[abcde]]<- rmse_our_all
    out_rmse_4D_abcde[[abcde]] <- rmse_4D_all
    out_rmse_smoother_abcde[[abcde]] <- rmse_smoother_all
  }
  
  #############################################################################################################
  ## leave one chain out
  
  
  # prepare all ids
  ids_all <- surfaces_test_ids
  
  
  # ids_all <- c(37)
  mse_our_all <- rmse_our_all <- rep(NA, length(ids_all))
  mse_4D_all <- rmse_4D_all <- rep(NA, length(ids_all))
  mse_smoother_all <- rmse_smoother_all <- rep(NA, length(ids_all))
  
  for (id_i in 1:length(ids_all)){
    id <- ids_all[id_i]
    # the surface to predict
    # plot(X_test[id,,], xlab="time to expiration", ylab="moneyness")
    surface_full <- X_test[id,,]
    
    ## leave one string out
    surface_condition <- surface_test <- surface_full
    strings_all <- (1:K2)[apply(apply(surface_condition, c(1,2), is.na ), 2, sum) < K2] # identify the strings
    strings_all_n <- length(strings_all)
    
    for (string_to_leave_i in 1:strings_all_n){
      
      string_to_leave <- strings_all[string_to_leave_i] # pick one string to leave out
      surface_condition[,string_to_leave] <- NA # remove that string from the condition
      
      
      ## the test partition
      surface_test[!is.na(surface_condition)] <- NA # delete obs. that are present in the condition
      # plot(surface_condition)
      # plot(surface_test)
      
      #############################################################################################################
      ## smooth first approach
      
      surface_condition_exp <- exp(surface_condition)
      # plot(surface_condition_exp)
      surface_condition_num_obs <- sum(!is.na(surface_condition_exp))
      
      # interpolate
      X_npreg <- matrix(NA, nrow=surface_condition_num_obs, ncol=2)
      Y_npreg <- rep(NA, surface_condition_num_obs)
      obs_i <- 1
      for (i in 1:K1){
        for (j in 1:K2){
          if (!is.na(surface_condition_exp[i,j])){
            X_npreg[obs_i,1] <- i
            X_npreg[obs_i,2] <- j
            Y_npreg[obs_i] <- surface_condition_exp[i,j]
            obs_i <- obs_i + 1
          }
        }
      }
      xy_grid <- as.matrix(expand.grid(1:K1, 1:K2))
      
      # continue only if there are at least 2 different expiration dates, and if at least 2 strikes
      if ((length(unique(X_npreg[,1])) >= 2) && (length(unique(X_npreg[,2])) >= 2)){
        
        
        kre <- npreg(bws=c(2,3), xdat=X_npreg, ydat=Y_npreg, exdat=xy_grid, regtype = "ll")
        surface_smoothed <- matrix(kre$mean, nrow=K1, ncol=K2)
        
        # RMSE denominator
        relativity <- sqrt(var(c(surface_full), na.rm = TRUE))
        
        # evaluate MSE for this surface
        test_mask <- !is.na(surface_test)
        if (is.na(mse_smoother_all[id_i])){mse_smoother_all[id_i] <- 0}
        if (is.na(rmse_smoother_all[id_i])){rmse_smoother_all[id_i] <- 0}
        mse_smoother_all[id_i] <- mse_smoother_all[id_i] + sqrt(mean((exp(surface_test[test_mask]) - surface_smoothed[test_mask])^2))/strings_all_n
        rmse_smoother_all[id_i] <- rmse_smoother_all[id_i] +
          sqrt(mean((exp(surface_test[test_mask]) - surface_smoothed[test_mask])^2))/strings_all_n / relativity
        
        
        #############################################################################################################
        ## predict the surface using our method
        surface_posterior <- predict_surface(surface_condition, prior_mu, prior_var, Res$sigma2)
        surface_posterior_4D <- predict_surface(surface_condition, prior_mu, prior_var_4D, Res_4D$sigma2)
        
        # separable error
        test_mask <- !is.na(surface_test)
        if (is.na(mse_our_all[id_i])){mse_our_all[id_i] <- 0}
        if (is.na(rmse_our_all[id_i])){rmse_our_all[id_i] <- 0}
        mse_our_all[id_i] <- mse_our_all[id_i] + sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior[test_mask]))^2 ))/strings_all_n
        rmse_our_all[id_i] <- rmse_our_all[id_i] +
          sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior[test_mask]))^2 ))/strings_all_n / relativity
        
        # 4D error
        test_mask <- !is.na(surface_test)
        if (is.na(mse_4D_all[id_i])){mse_4D_all[id_i] <- 0}
        if (is.na(rmse_4D_all[id_i])){rmse_4D_all[id_i] <- 0}
        mse_4D_all[id_i] <- mse_4D_all[id_i] + sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior_4D[test_mask]))^2 ))/strings_all_n
        rmse_4D_all[id_i] <- rmse_4D_all[id_i] +
          sqrt(mean( (exp(surface_test[test_mask]) - exp(surface_posterior_4D[test_mask]))^2 ))/strings_all_n / relativity
        
        
        
      }
    }
    
  }
  
  
  ## plot
  out_rmse_our_abcde[[1]] <- rmse_our_all
  out_rmse_4D_abcde[[1]] <- rmse_4D_all
  out_rmse_smoother_abcde[[1]] <- rmse_smoother_all


  out <- matrix(list(),1,1+3*5)
  out[[1]] <- surfaces_test_ids
  for (abcde in 1:5){
    out[[1+3*(abcde-1)+1]]   <- out_rmse_our_abcde[[abcde]]
    out[[1+3*(abcde-1)+2]] <- out_rmse_4D_abcde[[abcde]]
    out[[1+3*(abcde-1)+3]] <- out_rmse_smoother_abcde[[abcde]]
  }
  
  out
}


#########################################################################################
# convert out into what structure I want

ratios_abcde <- matrix(list(), 1, 5)
rmse_4D_abcde <- matrix(list(), 1, 5)
rmse_our_abcde <- matrix(list(), 1, 5)
rmse_smoother_abcde <- matrix(list(), 1, 5)
for (abcde in 1:5){
  ratios_abcde[[abcde]] <- rep(NA, surfaces_count)
  rmse_our_abcde[[abcde]] <- rep(NA, surfaces_count)
  rmse_4D_abcde[[abcde]] <- rep(NA, surfaces_count)
  rmse_smoother_abcde[[abcde]] <- rep(NA, surfaces_count)
}

for(k_fold_i in 1:k_fold_k){
  
  # out[[k_fold_i]] has the following elements
  # - list of test_id's
  # - (a) rmse_our
  # - (a) rmse_4D
  # - (b) etc.
  
  surfaces_test_ids <- out[[k_fold_i]][[1]]
  
  for (abcde in 1:5){
    rmse_our_abcde[[abcde]][surfaces_test_ids] <- out[[k_fold_i]][[1+3*(abcde-1)+1 ]]
    rmse_4D_abcde[[abcde]][surfaces_test_ids] <- out[[k_fold_i]][[ 1+3*(abcde-1)+2 ]]
    rmse_smoother_abcde[[abcde]][surfaces_test_ids] <- out[[k_fold_i]][[ 1+3*(abcde-1)+3 ]]
    ratios_abcde[[abcde]][surfaces_test_ids] <- out[[k_fold_i]][[ 1+3*(abcde-1)+1 ]] / out[[k_fold_i]][[ 1+3*(abcde-1)+2 ]] # our / 4D
  }
}


#########################################################################################
# abcde_names <- c("(a) Leave one chain out",
#                  "(b) Predict in-the-money options",
#                  "(c) Predict out-of-the-money options",
#                  "(d) Predict short maturities",
#                  "(e) Predict long maturities")
# par(mfrow=c(1,5))
# for (abcde in 1:5){
#   
#   # boxplot( 1/ratios_abcde[[abcde]], ylim=c(0,3) ) # 1/ratio > 1 ... our is better 
#   # 
#   # # separable vs 4D
#   boxplot( sqrt(rmse_our_abcde[[abcde]]) / sqrt(rmse_4D_abcde[[abcde]]), ylim=c(0,3) ) # <1 ... separable is better 
#   
#   # separable vs pre-smooth
#   # boxplot( sqrt(rmse_our_abcde[[abcde]]) / sqrt(rmse_smoother_abcde[[abcde]]), ylim=c(0,3) ) # <1 ... separable is better 
#   
#   title(abcde_names[abcde])
#   abline(h=1)
#   
# }
#  
# 
# 
# # mean on the level of RMSE
# for (abcde in 1:5){
#   print(mean(sqrt(rmse_our_abcde[[abcde]]), na.rm = T) / mean(sqrt(rmse_4D_abcde[[abcde]]), na.rm = T)) # <1 ... separable is better
# }
# 
# # mean, separable vs presmooth
# for (abcde in 1:5){
#   print(mean(sqrt(rmse_our_abcde[[abcde]]), na.rm = T) / mean(sqrt(rmse_smoother_abcde[[abcde]]), na.rm = T)) # <1 ... separable is better
# }
# 
# 
# # median on the level of RMSE
# for (abcde in 1:5){
#   print(median(sqrt(rmse_our_abcde[[abcde]]), na.rm = T) / median(sqrt(rmse_4D_abcde[[abcde]]), na.rm = T)) # <1 ... separable is better
# }
# 
# # median, separable vs presmooth
# for (abcde in 1:5){
#   print(median(sqrt(rmse_our_abcde[[abcde]]), na.rm = T) / median(sqrt(rmse_smoother_abcde[[abcde]]), na.rm = T)) # <1 ... separable is better
# }

save.image("grid20.Rdata")


