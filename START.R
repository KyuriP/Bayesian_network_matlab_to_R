
START <- function(DATA_ALL, steps, step_iterations, k_max, Run, k_transition, nue_var){
  
  fan_in = 3
  
  DAG             = Run$dag[[1]]         # the current DAG
  MATRIX          = Run$matrix[[1]]      # the current allocation matrix
  log_score       = Run$Log_Scores[[1]]  # the current log score
  lambda_snr_vec  = Run$lambda_snr_vec[[1]]
  lambda_coup_vec = Run$lambda_coup_vec[[1]]
  VECTORS         = Run$VECTORS[[1]]
  
  cat("##############################################################", "\n", 
      "An MCMC simulation for the EWC NH-DBN model has been started","\n",
      "##############################################################", "\n", 
      "Log-score of the initial graph:", round(log_score,5))
  
  ### Start of the MCMC simulation
  Counter = 2
  
  for (i in 1:steps){
    mcmc_output = MCMC_INNER(step_iterations, DAG, fan_in, log_score, MATRIX, DATA_ALL, k_max, nue_var, lambda_snr_vec, lambda_coup_vec, k_transition, VECTORS)
  
    DAG = mcmc_output[[1]]
    log_score = mcmc_output[[2]]
    MATRIX = mcmc_output[[3]]
    lambda_snr_vec = mcmc_output[[4]]
    lambda_coup_vec = mcmc_output[[5]]
    nue_var = mcmc_output[[6]]
    VECTORS = mcmc_output[[7]]
    
    Run$dag[[Counter]]              = DAG
    Run$matrix[[Counter]]           = MATRIX
    Run$Log_Scores[[Counter]]       = log_score
    Run$steps                       = i+1           # this updates the steps, check?
    Run$lambda_snr_vec[[Counter]]   = lambda_snr_vec
    Run$lambda_coup_vec[[Counter]]  = lambda_coup_vec
    Run$VECTORS[[Counter]]          = VECTORS
    
    Counter = Counter+1
  }
  return(Run)
}

#################################################################################
#################################################################################

MCMC_INNER <- function(step_iterations, DAG, fan_in, log_score, MATRIX, DATA_ALL, k_max, nue_var, lambda_snr_vec, lambda_coup_vec, k_transition, VECTORS){
  
  n_plus <- dim(DATA_ALL[[1]])[1]
  m      <- dim(DATA_ALL[[1]])[2]
  
  n_parents_keep <- dim(DAG)[1]
  n_nodes        <- dim(DAG)[2] 
  
  for (t in 1:step_iterations){
    for (i_node in 1:n_nodes){
      
      x = runif(1)
      
      ##### We need to initialize an empty list by the "number of nodes" and "H_max" later to fill in with the for-loop ... see line 127
      ##### # k_i cannot be larger than H_max (maximal number of data segments per gene) right??
      
      DATA     = replicate(n_nodes, vector(mode = "list", length = H_max), simplify = FALSE)
      DATA_NEW = replicate(n_nodes, vector(mode = "list", length = H_max), simplify = FALSE)   
      
      p_0 = 0.5 # structure MCMC
      p_1 = 0.7 # Birth Move
      p_2 = 0.9 # Death Move
      # else    # Reallocation Move
      
      if(x <= p_0){              # PERFORM A STRUCTURE-MCMC MOVE
        y = runif(1)      # x = rand y = rand(1) --> what is the difference?
        
        if(y < 0.5){
          
          child_node    = i_node
          old_parents   = as.matrix(DAG[ ,child_node]) # column vector
          new_parents = old_parents
    
          parents_card_old = sum(old_parents)
          
          parents_indicis_old = which(old_parents != 0)
          
          if(parents_card_old < fan_in){
            
            neighbours_old = n_parents_keep - 1
            
            indicis = sample(n_parents_keep)
            x_ind = indicis[1]
            
              if (x_ind == child_node){
                x_ind = indicis[2]
              }
            
            parent_index = x_ind # delete or add this parent node
            
            new_parents[parent_index, 1] = 1 - new_parents[parent_index, 1]
            
          }else{   # elseif (parent_card_old==fan_in)
            
            x_ind = sample(fan_in)
            x_ind = x_ind[1]
            
            parent_index = parents_indicis_old[x_ind] # delete this parent node
            
            new_parents[parent_index,1] = 0
            
            neighbours_old = parents_card_old     # = fan_in
          }
          
          parents_card_new = sum(new_parents)
          
          if (parents_card_new < fan_in){
            neighbours_new = n_parents_keep - 1
          } else {
            neighbours_new = parents_card_new # = fan_in
          }
          
          DAG_NEW = DAG
          DAG_NEW[ ,child_node] = new_parents
          
          
          data = DATA_ALL[[child_node]]
          k_i = max(MATRIX[child_node, ]) 
          for (component in 1:k_i){
            DATA[[child_node]][[component]] = data[ ,which(MATRIX[child_node, ]==component)]
          }      
          
          ###################################################
          VECTORS_NEW  = VECTORS  
          vector_i_new = VECTORS_NEW[[child_node]]
          
          if (new_parents[parent_index, 1] == 1){    # proposal is to add a parent
            vector_i_new[parent_index + 1 ,1] = as.numeric(runif(1) < 0.5)
            log_HR_supplement = log(2)
            
          } else{    # proposal is to delete a parent
            vector_i_new[parent_index + 1 ,1] = -1
            log_HR_supplement = -log(2)
          } 
          
          VECTORS_NEW[[child_node]] = vector_i_new
          ###################################################
          
          local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG_NEW, MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW)
          local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
          
          A = exp(local_score_new - local_score_old + log(neighbours_old) - log(neighbours_new) + log_HR_supplement)
          
          u = runif(1)
          
          if (u < A){      # accept the move:
            DAG       = DAG_NEW 
            VECTORS   = VECTORS_NEW
            
            log_score = log_score + local_score_new - local_score_old
          } 
          
          rm(DATA)
          
        } else {     # Perform an exchange move
          child_node    = i_node
          
          old_parents = which(DAG[ ,child_node] != 0) 
          
          parents_card_old = length(old_parents)
          
          if (parents_card_old==0){
            # do nothing
  
          } else { # perform an exchange move
          
          DAG_NEW = DAG
          
          indicis        = sample(parents_card_old)
          index          = indicis[1]
          
          parent_old_index = old_parents[index] # delete this parent node
          
          candidate_parents = which(DAG[ ,child_node] == 0)
          
          candidates_card = length(candidate_parents)
          
          indicis = sample(candidates_card)
          index   = indicis[1]
          
          parent_new_index = candidate_parents[index]
          
          if (parent_new_index == child_node){
            index = indicis[2]
          parent_new_index = candidate_parents[index]
          }
          
          DAG_NEW[parent_old_index, child_node] = 0
          DAG_NEW[parent_new_index, child_node] = 1
          
          data = DATA_ALL[[child_node]]
          k_i = max(MATRIX[child_node, ])
          
          for (component in 1:k_i){
            DATA[[child_node]][[component]] = data[ ,which(MATRIX[child_node, ]==component)]
          }   
          
          ###################################################
          VECTORS_NEW = VECTORS
          vector_i_new = VECTORS_NEW[[child_node]]
          
          vector_i_new[parent_old_index+1,1] = -1
          vector_i_new[parent_new_index+1,1] = as.numeric(runif(1) < 0.5)
          
          VECTORS_NEW[[child_node]] = vector_i_new
          
          log_HR_supplement = 0
          
          ###################################################
          
          local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG_NEW, MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW)
          local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
          
          log_hastings = 0
          
          A = exp(local_score_new - local_score_old + log_hastings + log_HR_supplement)
          
          u = runif(1)
          
          if (u < A){      # accept the move:
            DAG       = DAG_NEW
            VECTORS   = VECTORS_NEW
            log_score = log_score + local_score_new - local_score_old
          }
          rm(DATA)
          
          }
        } #########  if(y < 0.5) end
      
        
      } else if (x<=p_1){   # PERFORM A BREAKPOINT BIRTH MOVE
        
        if(max(MATRIX[i_node, ]) < k_max){
          
          NEW_CANDIDATES = matrix(1, 1, m-1)
          
          break_vec     = which(MATRIX[i_node, 2:ncol(MATRIX)] - MATRIX[i_node, 1:(ncol(MATRIX)-1)] != 0)
          break_vec_add = matrix(c(0, break_vec, m), nrow=1)
          

          for (i in break_vec_add){  
            NEW_CANDIDATES[1, max(c(i-(k_transition-1), 1)) : min(c(i+(k_transition-1)),m-1)] = 0
          }
          
          NEW_CANDIDATES = which(NEW_CANDIDATES==1)
          
          n_candidates = length(NEW_CANDIDATES)
          
          if (n_candidates > 0){
            indicis      = sample(n_candidates)
            index        = indicis[1]
            
            index        = NEW_CANDIDATES[index]
            
            MATRIX_NEW = MATRIX
            MATRIX_NEW[i_node, (index+1):ncol(MATRIX_NEW)] = MATRIX_NEW[i_node,(index+1):ncol(MATRIX_NEW)] + 1
            
            rm(DATA)
            rm(DATA_NEW)
            
            data = DATA_ALL[[i_node]]
            
            ###### I am confused... in Matlab your remove the variable but it automatically initialize when just calling it?...
            
            for(component in 1:max(MATRIX[i_node, ])){
              DATA[[i_node]][[component]] = data[ ,which(MATRIX[i_node, ] == component)]
            }   
            
            for(component in 1:max(MATRIX_NEW[i_node, ])){
              DATA_NEW[[i_node]][component] = data[ ,which(MATRIX_NEW[i_node, ] == component)]
            }   
            
            ###################################################
            
            vector_i_new                  = VECTORS[[i_node]]
            
            parents   = which(DAG[,i_node] != 0)
            n_parents = length(parents)
            
            vector_i_new[c(1, parents+1),1] = matrix(as.numeric(runif(n_parents+1) < 0.5), n_parents+1, 1)
            
            log_HR_supplement  = 0
            
            VECTORS_NEW         = VECTORS
            VECTORS_NEW{i_node} = vector_i_new 
            ###################################################
            
            local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_NEW, DAG, MATRIX_NEW, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW)
            local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX,     i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
            
            n_breakpoints_new = max(MATRIX_NEW[i_node, ])
            
            log_hastings = log_HR_supplement + log(n_candidates) - log(n_breakpoints_new)
            
            A = exp(local_score_new - local_score_old + log_hastings)
            
            u = runif(1)
            
            if (u < A){     # accept the move:
              MATRIX    = MATRIX_NEW
              VECTORS   = VECTORS_NEW
              log_score = log_score + local_score_new - local_score_old
            } 
          }
        }
        rm(DATA)
        rm(DATA_NEW)
      }      ##### if(x <= p_0) --> continue .. first else if(x<=p_1) ends here
      
      else if (x<=p_2){     # PERFORM A BREAKPOINT DEATH MOVE
    
        CANDIDATES = which(MATRIX[i_node, 2:ncol(MATRIX)] - MATRIX[i_node, 1:ncol(MATRIX)-1] != 0)
        
        if(length(CANDIDATES)==0){  # then there is no breakpoint which can be removed
          
        } else {
          n_candidates = length(CANDIDATES)
          indicis = sample(n_candidates)
          index   = indicis[1]
          
          candidate = CANDIDATES[index]
          
          MATRIX_NEW                           = MATRIX
          MATRIX_NEW[i_node, (candidate+1):ncol(MATRIX_NEW)] = MATRIX_NEW[i_node, (candidate+1):ncol(MATRIX_NEW)] - 1
          
          rm(DATA)
          rm(DATA_NEW)
          
          data = DATA_ALL[[i_node]]
          
          for (component in 1:max(MATRIX[i_node, ])){     ###### probably need to DEFINE DATA and DATA NEW AGAIN HERE BECUZ THEY ARE RMEMOVED ABOVE... ######
            DATA[[i_node]][component] = data[ , which(MATRIX[i_node, ] == component)]
          } 
          
          for (component in 1:max(MATRIX_NEW[i_node, ])){
            DATA_NEW[[i_node]][component] = data[ ,which(MATRIX_NEW[i_node, ]==component)]
          } 

          #######################################################
          vector_i_new                  = VECTORS[[i_node]]
          
          parents   = which(DAG[ ,i_node] != 0)
          n_parents = length(parents)
          
          vector_i_new[c(1, parents+1),1] = matrix(as.numeric(runif(n_parents+1) < 0.5), n_parents+1, 1)
          
          log_HR_supplement  = 0
          
          #######################################################
         
          VECTORS_NEW = VECTORS
          VECTORS_NEW[[i_node]] = vector_i_new
          
          
          local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_NEW, DAG, MATRIX_NEW, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW)
          local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX,     i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
          
          BIRTH_CANDIDATES = matrix(1, 1, m-1)
          break_vec_new = find(MATRIX_NEW[i_node, 2:ncol(MATRIX_NEW)] - MATRIX_NEW[i_node, 1:ncol(MATRIX_NEW)-1])
          
          break_vec_new_add = matrix(c(0, break_vec_new, m), nrow = 1)
          
          for (i in break_vec_new_add){
            BIRTH_CANDIDATES[1, max(c(i-(k_transition-1), 1)) : min(c(i+(k_transition-1)),m-1)] = 0
          }

          
          BIRTH_CANDIDATES   = which(BIRTH_CANDIDATES==1)
          n_birth_candidates = length(BIRTH_CANDIDATES)
          
          log_hastings =  log_HR_supplement + log(n_candidates) - log(n_birth_candidates)
          
          A = exp(local_score_new - local_score_old + log_hastings)
          
          u = unif(1)
          
          if (u < A){# accept the move:
            MATRIX      = MATRIX_NEW 
            VECTORS     = VECTORS_NEW
            log_score   = log_score + local_score_new - local_score_old
          } 
          
          rm(DATA)
          rm(DATA_NEW)
          
        }
        
      }  ##### if(x <= p_0) --> continue .. second else if(x<=p_2) ends here 
      else {     #PERFORM A BREAKPOINT REALLOCATION MOVE
        
        CANDIDATES = which(MATRIX[i_node, 2:ncol(MATRIX)] - MATRIX[i_node, 1:ncol(MATRIX)-1] != 0)
        n_candidates = length(CANDIDATES)
        
        if(n_candidates > 0){    # then there are breakpoints which can be re-allocated
          
          indicis   = sample(n_candidates)
          index_old = indicis[1]
          candidate = CANDIDATES[index_old]
          
          CANDIDATES_ADD = matrix(c(0,CANDIDATES,m), nrow=1)
          
          index_candidate = which(CANDIDATES_ADD == candidate)
          
          k_minus = CANDIDATES_ADD[index_candidate-1]
          k_i     = CANDIDATES_ADD[index_candidate]
          k_plus  = CANDIDATES_ADD[index_candidate+1] 
          
          k_i_new_possible = (k_minus + k_transition):(k_plus-k_transition)
          n_candidates_new = length(k_i_new_possible)
          
          indicis_new = sample[n_candidates_new]
          index_new   = indicis_new[1]
          
          candidate_new = k_i_new_possible[index_new]
          
          MATRIX_NEW                                = MATRIX
          MATRIX_NEW[i_node, (candidate+1):ncol(MATRIX_NEW)]     = MATRIX_NEW[i_node, (candidate+1):ncol(MATRIX_NEW)] - 1
          MATRIX_NEW[i_node, (candidate_new+1):ncol(MATRIX_NEW)] = MATRIX_NEW[i_node, (candidate_new+1):ncol(MATRIX_NEW)] + 1
          
          rm(DATA)
          rm(DATA_NEW)
          
          data = DATA_ALL[[i_node]]
          
          ############ AGAIN PROABLY DATA AND DATA_NEW SHOULD BE REDIFINED SINCE THEY ARE REMOVED ABOVE
          for (component in 1:max(MATRIX[i_node, ])){
            DATA[[i_node]][component] = data[ ,which(MATRIX[i_node, ]==component)]
          }   
        
          for (component in 1:max(MATRIX_NEW[i_node, ])){
            DATA_NEW[[i_node]][component] = data[ ,which(MATRIX_NEW[i_node, ]==component)]
          } 
          
          
          #######################################################
          
          vector_i_new                  = VECTORS[[i_node]]
          
          parents   = which(DAG[ ,i_node] != 0)
          n_parents = length(parents)
          
          vector_i_new[c(1, parents+1),1] = matrix(as.numeric(runif(n_parents+1) < 0.5), n_parents+1, 1)
          
          log_HR_supplement  = 0
          
          #######################################################
          
          VECTORS_NEW         = VECTORS
          VECTORS_NEW[[i_node]] = vector_i_new
          
          #######################################################
          
          local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_NEW, DAG, MATRIX_NEW, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW)
          local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX,i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
          
          A = exp(local_score_new - local_score_old + log_HR_supplement)
          
          u = runif(1)
          
          if (u < A){# accept the move:
            MATRIX      = MATRIX_NEW 
            VECTORS     = VECTORS_NEW
            log_score   = log_score + local_score_new - local_score_old
          } 
        }
      }
      rm(DATA)
      rm(DATA_NEW)
    }  ## i_node inner forloop ends   % MOVE-TYPE-LOOP
    
    # Update hyperparameters:
      
    rm(DATA)
    ## HERE THEY AGAIN DEFINE AN EMPTY MATRIX (in matlab)
    DATA = replicate(n_nodes, vector(mode = "list", length = H_max), simplify = FALSE)
    
    for (node_i in 1:n_nodes){
      data = DATA_ALL[[node_i]] 
      for (component in 1:max(MATRIX[node_i, ])){
        DATA[[node_i]][[component]] = data[ ,which(MATRIX[node_i,]==component)]
      } 
    }

    update_output = UPDATE(DATA, DAG, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
    lambda_snr_vec = update_output[[1]]
    lambda_coup_vec = update_output[[2]]
    VECTORS = update_output[[3]]
    
    # ADDITIONAL VECTOR_MOVES
    
    for (i_node in 1:n_nodes){
      
      VECTORS_NEW = VECTORS
      vector_i    = VECTORS[[i_node]]
      
      parents   = which(DAG[ ,i_node] != 0)
      n_parents = length(parents)
      
      indicis = sample(n_parents+1)
      index   = indicis[1]
      
      coefficients = matrix(c(1,parents+1), ncol = 1)
      coefficient  = coefficients[index]
      
      vector_i_new = vector_i
      vector_i_new[coefficient, 1] = 1-vector_i_new[coefficient, 1]
      
      VECTORS_NEW[[i_node]] = vector_i_new
      
      ##### check this whether the function assigns the output correctly!! 
      local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW)
      local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
      
      log_hastings = 0
      
      A = exp(local_score_new - local_score_old + log_hastings)
      u = runif(1)
      
      if (u<A) {# accept the move:
        VECTORS   = VECTORS_NEW
        log_score = log_score + local_score_new - local_score_old     
      }
    }
    
    log_score = COMPUTE_LOG_SCORE(DATA, DAG, MATRIX, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)
  } # t forloop ends
  
  output <- list(DAG, log_score, MATRIX, lambda_snr_vec, lambda_coup_vec, nue_var, VECTORS)
  return (output)
} # MCMC function ends



###########################################################################
###########################################################################
COMPUTE_LOCAL_LOG_SCORE <- function(DATA, DAG, MATRIX, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS){
  
  # global Prior     ## make Prior a global variable
  
  n_nodes <- dim(MATRIX)[1]
  m <- dim(MATRIX)[2]
  
  k = length(DATA[[i_node]])    
  
  log_prob_k <- log(dpois(k,1))
  
  k_cps = k-1
  
  breakpoints <- which((MATRIX[i_node, 2:ncol(MATRIX)] - MATRIX[i_node, 1:(ncol(MATRIX)-1)])!=0) 
  
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
  
  ###########################################################################
  
  log_prob_graph =  Prior[length(which(DAG[ ,i_node] != 0))+1]
  
  ###########################################################################
  
  # Compute the local score for i_node:
  
  k_i = length(!is.null(DATA[[i_node]]))  ##### try to only count the non-NULL (this is because we initialize the empty list above by the number of nodes and H_max)
  
  parents   = which(DAG[,i_node] != 0)
  n_parents = length(parents)
  
  sum_log_det_Sigma_tilde = 0
  sum_Delta2              = 0
  
  vector_i = VECTORS[[i_node]]
  
  ind1 = which(vector_i==1) 
  ind0 = which(vector_i==0) 
  
  lambda_coup = lambda_coup_vec[i_node,1]
  lambda_snr  = lambda_snr_vec[i_node,1]
  
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
        mue_prior = matrix(0, length(parents)+1, 1) # pred x 1
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
      
      Sigma_tilde = diag(n_obs) + (t(X) %*% LAMBDA) %*% X  ### not sure if here %*% or * is correct.
      
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

  #########################
  log_prob_VECTOR = (n_parents + 1) * log(0.5)
  #########################
  
  log_score = log_prob_breaks + log_prob_graph + log_score_i + log_prob_VECTOR
  
  return (log_score)
} # COMPUTE_LOCAL_LOG_SCORE function ends
    




###########################################################################
## This is the same function as in INITIALIZE.R... right? why define twice here?
###########################################################################
COMPUTE_LOG_SCORE <- function(DATA, DAG, MATRIX, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS){
  
  #global Prior; Prior needs to be a global variable!
  
  log_prob_breaks <- 0
  
  n_nodes <- dim(MATRIX)[1]
  m <- dim(MATRIX)[2]
  
  for (i_node in 1:n_nodes){
    k <- length(DATA[[i_node]]) 
    log_prob_k <- log(dpois(k,1))
    k_cps <- k-1
    breakpoints <- which((MATRIX[i_node, 2:ncol(MATRIX)] - MATRIX[i_node, 1:(ncol(MATRIX)-1)])!=0) 
    
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
    k_i = length(DATA[[i_node]]) #k_i is the number of mixture components for the i-th node
    
    parents = which(DAG[ ,i_node] != 0)
    
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
        
        Sigma_tilde = diag(n_obs) + (t(X) %*% LAMBDA) %*% X  ### not sure if here %*% or * is correct.
        
        # pred * obs
        inv_Sigma_tilde = diag(n_obs) - t(X) %*% solve(solve(LAMBDA) + X %*% t(X)) %*% X
        
        # (1 x obs) * (obs x obs) * (obs x 1) 
        sum_Delta2 = sum_Delta2 + (t(y - m_tilde) %*% inv_Sigma_tilde %*% (y - m_tilde))
        
        sum_log_det_Sigma_tilde = sum_log_det_Sigma_tilde + log(det(Sigma_tilde))
        
        # pred * pred
        Sigma_inv = solve(LAMBDA) + X %*% t(X)  # pred * pred 
        
        # pred * 1
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
  
}# COMPUTE_LOG_SCORE function ends
