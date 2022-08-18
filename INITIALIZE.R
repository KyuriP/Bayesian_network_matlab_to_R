
INITIALISE <- function(DATA_ALL, DAG, steps, nue_var, lambda_snr_vec, lambda_coup_vec, MATRIX, VECTORS){
  
  n_plus <- dim(DATA_ALL[[1]])[1]
  m <- dim(DATA_ALL[[1]])[2]
  
  n_nodes <- length(DATA_ALL)
  
  
  DATA <- list()
  for (node_i in 1:n_nodes){
    data <- DATA_ALL[[node_i]]
    nested_DATA <- list()
    for (component in 1:max(MATRIX[node_i, ])){
      nested_DATA[[component]] <- data[ , which(MATRIX[node_i, ] == component)]
    } 
    DATA[[node_i]] <- nested_DATA
  }
  
  log_score <- COMPUTE_LOG_SCORE(DATA, DAG, MATRIX, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
  
  Run <- list()
  for (i in 1:(steps+1)){
    Run$dag[[i]]             <- 0
    Run$matrix[[i]]          <- 0   ## Check the data structure for each again!
    Run$Log_Scores[[i]]      <- 0
    Run$lambda_snr_vec[[i]]  <- 0
    Run$lambda_coup_vec[[i]] <- 0
    Run$VECTORS[[i]]         <- 0
  }
  
  # Initialisation:
  Run$dag[[1]] <- DAG
  Run$matrix[[1]] <- MATRIX
  Run$Log_Scores[[1]] <- log_score
  Run$steps <- 1                     
  Run$lambda_snr_vec[[1]]  <- lambda_snr_vec
  Run$lambda_coup_vec[[1]] <-  lambda_coup_vec
  Run$VECTORS[[1]] <- VECTORS
  
  return(Run)
}

##########################################################
##########################################################
COMPUTE_LOG_SCORE <- function(DATA, DAG, MATRIX, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS) {
  
  #global Prior; Prior needs to be a global variable!
  
  log_prob_breaks <- 0
  
  n_nodes <- dim(MATRIX)[1]
  m <- dim(MATRIX)[2]
  
  for (i_node in 1:n_nodes){
    k <- length(DATA[[i_node]]) 
    log_prob_k <- log(dpois(k,1))
    k_cps <- k-1
    
    breakpoints <- which((MATRIX[i_node, 2:ncol(MATRIX)] - MATRIX[i_node, 1:(ncol(MATRIX)-1)])!=0) # does this work?
    if (length(breakpoints) == 0){
      log_prob_break <- 0
    } else {
      breakpoints <- matrix( c(0, breakpoints, m), nrow=1)
      log_prob_break <- log(prod(1:(2*k_cps+1))) - log(prod(((m-1)-(2*k_cps+1)+1):(m-1)))
      for (i in 2:length(breakpoints)){
        log_prob_break <- log_prob_break + log(breakpoints[i]-breakpoints[i-1] - 1)
      }
    }
    log_prob_breaks <- log_prob_breaks + log_prob_break + log_prob_k
  }
  
  #################################################################################
  
  log_prob_graph = 0
  
  for (node in 1:n_nodes){
    log_prob_graph = log_prob_graph + Prior[length(which(DAG[,node] != 0))+1]
  }
  
  #################################################################################
  
  log_prob_data = 0
  
  for (i_node in 1:n_nodes){
    k_i = length(DATA[[i_node]])
    
    parents = which(DAG[ ,i_node] != 0)
    
    lambda_coup = lambda_coup_vec[i_node, 1]
    lambda_snr  = lambda_snr_vec[i_node, 1]
    
    sum_log_det_Sigma_tilde = 0
    sum_Delta2              = 0
    
    vector_i = VECTORS[[i_node]]
    
    ind1 = which(vector_i==1)
    ind0 = which(vector_i==0)
    
    LAMBDA_VEC       = vector_i
    LAMBDA_VEC[ind0] = lambda_snr
    LAMBDA_VEC[ind1] = lambda_coup
    
    LAMBDA_MAT = diag(drop(LAMBDA_VEC))
    LAMBDA_MAT = LAMBDA_MAT[c(1,parents+1),c(1,parents+1)] 
    
    ### FOR THE FIRST SEGMENT:
    
    LAMBDA_VEC_first       = vector_i
    LAMBDA_VEC_first[ind0] = lambda_snr
    LAMBDA_VEC_first[ind1] = lambda_snr
    
    LAMBDA_MAT_first = diag(drop(LAMBDA_VEC_first))
    LAMBDA_MAT_first = LAMBDA_MAT_first[c(1,parents+1),c(1,parents+1)]  
    
    for (component in 1:k_i){
      data = DATA[[i_node]][component][[1]]  # to make data as a matrix [[1]]
      
      n_plus <- dim(data)[1]
      n_obs <- dim(data)[2]
      
      if (n_obs == 0){
        # do nothing
      } else {
        X = rbind(matrix(1, 1, n_obs), data[parents,])
        y = as.matrix(data[nrow(data),])  # transpose  no need?
        
        if (component == 1){
          mue_prior = matrix(0, length(parents)+1, 1)
          LAMBDA    = LAMBDA_MAT_first
        } else {
          if (length(parents) > 0) {
            mue_prior = vector_i[c(1, parents+1), 1] * mue_apost
          } else {
            mue_prior = vector_i[1,1] * mue_apost
          }
          LAMBDA = LAMBDA_MAT
        }
 
        m_tilde  = t(X) %*% mue_prior    
        
        Sigma_tilde = diag(n_obs) + (t(X) * LAMBDA) %*% X  ### not sure if here %*% or * is correct.
        
        # pred * obs
        inv_Sigma_tilde = diag(n_obs) - t(X) %*% solve(solve(LAMBDA) + X %*% t(X)) %*% X
        
        # (1 x obs) * (obs x obs) * (obs x 1) 
        sum_Delta2 = sum_Delta2 + (t(y - m_tilde) %*% inv_Sigma_tilde %*% (y - m_tilde))
        
        sum_log_det_Sigma_tilde = sum_log_det_Sigma_tilde + log(det(Sigma_tilde))
        
        Sigma_inv = solve(LAMBDA) + X %*% t(X)  # pred * pred 
        
        mue_apost  = solve(Sigma_inv) %*% (solve(LAMBDA) %*% mue_prior + X %*% y) # pred x 1
      }
    }
    sum_1 = log(gamma((m + nue_var)/2)) - log(gamma(nue_var/2))
    
    sum_2 = (nue_var/2) * log(nue_var) - (m/2)*log(pi) - 0.5 * sum_log_det_Sigma_tilde
    
    sum_3 = -(m + nue_var)/2 * log(nue_var + sum_Delta2)
    
    log_score_i = sum_1 + sum_2 + sum_3
    
    log_prob_data = log_prob_data + log_score_i
  }
  
  #################################################################################
  # global alpha_snr;                 # They all need to be global variables! 
  # global beta_snr; 
  # global alpha_coup;
  # global beta_coup;
  
  
  log_prob_lambda = 0
  
  for (i_node in 1:n_nodes){
    log_prob_lambda_snr_i  = log(dgamma(1/lambda_snr_vec[i_node, 1], alpha_snr, scale = (1/beta_snr)))
    log_prob_lambda_coup_i = log(dgamma(1/lambda_coup_vec[i_node, 1], alpha_coup, scale = (1/beta_coup)))
    log_prob_lambda        = log_prob_lambda + log_prob_lambda_snr_i + log_prob_lambda_coup_i
    
  }
  
  
  #################################################################################
  
  log_prob_VECTOR = (sum(sum(DAG)) + n_nodes) * log(0.5)
  
  #################################################################################
  
  log_score = log_prob_breaks + log_prob_graph + log_prob_data + log_prob_lambda + log_prob_VECTOR
  
  return (log_score)
}

## need to check whether the matrix multiplication works. this is just vague draft. compare with matlab by stopping each step !!!
