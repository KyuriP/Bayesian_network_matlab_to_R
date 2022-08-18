
PROC_YEAST <- function (data_1, data_2, steps, step_iterations, k_max, k_transition, lambda_coup, lambda_snr, nue_var, MATRIX){
  
  ###### global DATA_ALL    ### DATA_ALL needs to be global!
  DATA_1 <- SHIFT_DATA(data_1)
  DATA_2 <- SHIFT_DATA(data_2)
  
  DATA_ALL <- list()
  for (i in 1:length(DATA_1)){               
    DATA_ALL[[i]] <- cbind(DATA_1[[i]], DATA_2[[i]])
  }
  
  ##### preliminary check: it works ######
  # D1 <- SHIFT_DATA(data_on)
  # D2 <- SHIFT_DATA(data_off)
  # 
  # DATA_ALL <- list()
  # for (i in 1:length(D1)){               
  #   DATA_ALL[[i]] <- cbind(D1[[i]], D2[[i]])
  # }
  

  ### global Prior  #### Prior needs to be global!
  
  Prior <- SET_HYPERPARAMETERS(DATA_ALL) 
  
  n_nodes <- length(DATA_ALL)    
  
  DAG <- matrix(0, n_nodes, n_nodes)
  
  lambda_coup_vec <- lambda_coup * matrix(1, n_nodes, 1) 
  lambda_snr_vec  <- lambda_snr  * matrix(1, n_nodes, 1) 
  
  ###########################################################################
  VECTORS <- list()
  for (i_node in 1:n_nodes){
    VECTORS[[i_node]] <- rbind(0, -1 * matrix(1, n_nodes, 1))
  }
  
  ###########################################################################
  
  [Run] = INITIALISE(DATA_ALL, DAG, steps, nue_var, lambda_snr_vec, lambda_coup_vec, MATRIX, VECTORS);
  [Run] = START(DATA_ALL, steps, step_iterations, k_max, Run, k_transition, nue_var);
  
  return;
}

###########################################################################
###########################################################################
SHIFT_DATA <- function(data){
  n <- dim(data)[1]
  t <- dim(data)[2]
  n_plus <- n + 1 # Covariance matrices will be of size (n+1)x(n+1)
  
  DATA_SHIFTED <- list()
  for (node_i in 1:n){
    # For each variable X_i consider
    # data_1[X_1(t-1),...,X_(n)(t-1),...,X_1(t-slice),...,X_(n)(t-slice), X_(i)(t)] 
    obs <- 1:(t-1)  
    data_new <- matrix(0, n_plus, (t-1))   
    data_new[1:n, obs] <- data[1:n, 1:(t-1)]       
    data_new[n+1, obs] <- data[node_i, 2:t]   
    DATA_SHIFTED[[node_i]] <- data_new     
  }
  return(DATA_SHIFTED)
}

###########################################################################
###########################################################################
SET_HYPERPARAMETERS <- function(DATA_ALL){
  n_plus <- dim(DATA_ALL[[1]])[1]
  n_obs <- dim(DATA_ALL[[1]])[2]
  n_nodes <- n_plus - 1
  Prior <- matrix(0, (n_nodes + 1), 1)
  
  return(Prior)
}

