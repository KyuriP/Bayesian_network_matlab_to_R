
UPDATE <- function(DATA, DAG, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS){
  
  ################## MAKE THEM GLOBAL! 
  # global alpha_snr;
  # global beta_snr;
  # 
  # global alpha_coup;
  # global beta_coup;
  
  n_nodes = length(DATA)
  
  for (i_node in 1:n_nodes){
    
    vector_i = VECTORS[[i_node]]
    
    n_comps  = length(DATA[[i_node]])
    parents  = which(DAG[,i_node] != 0)
    
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
    
    alpha_sigma = nue_var/2
    beta_sigma  = nue_var/2
    
    for (component in 1:n_comps){
      data = DATA[[i_node]][component][[1]]  # to make data as a matrix [[1]]
      
      n_plus <- dim(data)[1]
      n_obs <- dim(data)[2]
      
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
      
      m_tilde  = t(X) %*% mue_prior # obs x 1
      
      # obs x obs       
      inv_Sigma_tilde = diag(n_obs) - t(X) %*% solve(solve(LAMBDA) + X %*% t(X)) %*% X
      
      Delta2       = t(y-m_tilde) %*% inv_Sigma_tilde %*% (y-m_tilde)
      
      alpha_sigma = alpha_sigma +  n_obs/2
      beta_sigma  = beta_sigma  + Delta2/2
      
      # pred x pred              
      Sigma_inv = solve(LAMBDA) + X %*% t(X)
      
      mue_apost  = solve(Sigma_inv) %*% (solve(LAMBDA) %*% mue_prior + X %*% y) 
    } # for loop component ends
    
    inv_var_all = rgamma(n = 1, shape = alpha_sigma, scale = (1/beta_sigma))
    var_all     = 1/inv_var_all
    
    ######################################################################
    ######################################################################  
    ######################################################################
    
    alpha_snr_i = alpha_snr
    beta_snr_i  = beta_snr
    
    alpha_coup_i = alpha_coup
    beta_coup_i  = beta_coup
    
    ######################################################################
    
    
    alpha_snr_i  = alpha_snr_i  + (length(ind1) * 1  +  length(ind0) * n_comps)/2
    alpha_coup_i = alpha_coup_i + (length(ind1) * (n_comps-1))/2
    
    ######################################################################
    
    for (component in 1:n_comps){
      
      data = DATA[[i_node]][component][[1]]  # to make data as a matrix [[1]]
      
      n_plus <- dim(data)[1]
      n_obs <- dim(data)[2]
      
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
      
      # pred x pred              
      Sigma_inv = solve(LAMBDA) + X %*% t(X)  # pred * pred 
      
      mue_apost  = solve(Sigma_inv) %*% (solve(LAMBDA) %*% mue_prior + X %*% y) 
      ########## no need to transpose, as it is already a column vector now (no need to specify "n" either, mvrnorm works similar way as mvnrnd in matlab and gives the same dimension as mue_apost)
      W_i = as.matrix(MASS::mvrnorm(mu = mue_apost, Sigma = var_all*solve(Sigma_inv))) 

      
      if (component==1){
        beta_snr_i  = beta_snr_i  + 0.5 * inv_var_all * t(W_i - mue_prior) %*% solve(diag(length(parents)+1)) %*% (W_i - mue_prior) 
      } else {
        
        indi = which(ifelse(vector_i %in% -1, 0, 1) !=0 )
        vec_i = vector_i[indi]
        
        ind0_new = which(vec_i==0)
        ind1_new = which(vec_i==1)
        
        W_i_snr  = W_i[ind0_new]
        W_i_coup = W_i[ind1_new]                   # HERE CHECK ALL IF THE INDEXING WORKED FINE!!!!! ALL VETORS if just [] works or do i need to specify [,] comma somwhere...
        
        mue_prior_snr  = mue_prior[ind0_new]
        mue_prior_coup = mue_prior[ind1_new]
        
        n0 = length(ind0)
        n1 = length(ind1)
        
        if(n0>0){
          beta_snr_i  = beta_snr_i  + 0.5 * inv_var_all * t(W_i_snr - mue_prior_snr) %*% solve(diag(n0)) %*% (W_i_snr - mue_prior_snr)
          
        }
        if(n1>0){
          beta_coup_i = beta_coup_i + 0.5 * inv_var_all * t(W_i_coup-mue_prior_coup) %*% solve(diag(n1)) %*% (W_i_coup - mue_prior_coup)
        }
        
      }
      
    } # component for-loop ends
    
    inv_lambda_coup = rgamma(n = 1, shape = alpha_coup_i, scale = (1/beta_coup_i))
    lambda_coup_vec[i_node, 1] = (1/inv_lambda_coup)
    
    inv_lambda_snr            = rgamma(n = 1, shape = alpha_snr_i, scale = (1/beta_snr_i))
    lambda_snr_vec[i_node, 1]  = (1/inv_lambda_snr)
    
    
  } # for loop i_node ends
  output <- list(lambda_snr_vec, lambda_coup_vec, VECTORS)
  return(output)  
}  # UPDATE function ends


