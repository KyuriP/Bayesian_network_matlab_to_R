
function [lambda_snr_vec, lambda_coup_vec, VECTORS] = UPDATE(DATA, DAG, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)

global alpha_snr;
global beta_snr;

global alpha_coup;
global beta_coup;

n_nodes = length(DATA);

for i_node=1:n_nodes
    
    vector_i = VECTORS{i_node};
    
    n_comps  = length(DATA{i_node}); 
    parents  = find(DAG(:,i_node));
    
    ind1 = find(vector_i==1);
    ind0 = find(vector_i==0);
    
    lambda_coup = lambda_coup_vec(i_node,1);
    lambda_snr  = lambda_snr_vec(i_node,1);

    LAMBDA_VEC       = vector_i;
    LAMBDA_VEC(ind0) = lambda_snr;
    LAMBDA_VEC(ind1) = lambda_coup;
    
    LAMBDA_MAT = diag(LAMBDA_VEC);
    LAMBDA_MAT = LAMBDA_MAT([1;parents+1],[1;parents+1]);
    
    %%% FOR THE FIRST SEGMENT:

    LAMBDA_VEC_first       = vector_i;
    LAMBDA_VEC_first(ind0) = lambda_snr; %%%
    LAMBDA_VEC_first(ind1) = lambda_snr; %%%

    LAMBDA_MAT_first = diag(LAMBDA_VEC_first);
    LAMBDA_MAT_first = LAMBDA_MAT_first([1;parents+1],[1;parents+1]);
    
    alpha_sigma = nue_var/2;
    beta_sigma  = nue_var/2;
       
     for component=1:n_comps
         
        data = DATA{i_node}{component};
    
        [n_plus, n_obs] =  size(data);
   
        X = [ones(1,n_obs);data(parents,:)]; % pred x obs 
        y = data(end,:)';                    % obs x 1
         
        if (component==1)
            mue_prior = zeros(length(parents)+1,1); % pred x 1  
            LAMBDA    = LAMBDA_MAT_first;
        else
              if (length(parents)>0)
                    mue_prior = vector_i([1;parents+1],1) .* mue_apost;
              else
                    mue_prior = vector_i(1,1)             .* mue_apost;
              end
              
              LAMBDA = LAMBDA_MAT;
        end
         
        m_tilde = X'*mue_prior; % obs x 1
                 
        % obs x obs       
        inv_Sigma_tilde = eye(n_obs) - X'*inv(inv(LAMBDA)+X*X')*X;
                
        Delta2       = (y-m_tilde)'*inv_Sigma_tilde*(y-m_tilde); 
                
        alpha_sigma = alpha_sigma +  n_obs/2;
        beta_sigma  = beta_sigma  + Delta2/2;
    
        % pred x pred              
        Sigma_inv = inv(LAMBDA) + X*X';                        
   
        mue_apost  = inv(Sigma_inv)*(inv(LAMBDA)*mue_prior+X*y);
          
     end
    
    inv_var_all = gamrnd(alpha_sigma,(1/beta_sigma));
    var_all     = 1/inv_var_all;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    alpha_snr_i = alpha_snr;
    beta_snr_i  = beta_snr;
    
    alpha_coup_i = alpha_coup;
    beta_coup_i  = beta_coup;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    alpha_snr_i  = alpha_snr_i  + (length(ind1) * 1  +  length(ind0) * n_comps)/2; 
    alpha_coup_i = alpha_coup_i + (length(ind1) * (n_comps-1))/2; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for component=1:n_comps
        
        data = DATA{i_node}{component};
        [n_plus, n_obs] =  size(data);
   
        X = [ones(1,n_obs);data(parents,:)]; % pred x obs 
        y = data(end,:)';                    % obs  x 1
      
         if (component==1)
            mue_prior = zeros(length(parents)+1,1); % pred x 1
            LAMBDA    = LAMBDA_MAT_first;
         else
             if (length(parents)>0)
                    mue_prior = vector_i([1;parents+1],1) .* mue_apost;
             else
                    mue_prior = vector_i(1,1) .* mue_apost;
             end
             
             LAMBDA = LAMBDA_MAT;
         end
        
            % pred x pred              
            Sigma_inv = inv(LAMBDA) + X*X';                        
   
            mue_apost  = inv(Sigma_inv)*(inv(LAMBDA)*mue_prior+X*y);
          
            W_i = mvnrnd(mue_apost,var_all*inv(Sigma_inv));
            W_i = W_i';
            
            if (component==1)
                beta_snr_i  = beta_snr_i  + 0.5 * inv_var_all*(W_i-mue_prior)'*inv(eye(length(parents)+1))*(W_i-mue_prior);
            else
                  
               indi = find(vector_i~=(-1));
               vec_i = vector_i(indi);
               
               ind0_new = find(vec_i==0);
               ind1_new = find(vec_i==1);

               W_i_snr  = W_i(ind0_new);
               W_i_coup = W_i(ind1_new);
                
               mue_prior_snr  = mue_prior(ind0_new);
               mue_prior_coup = mue_prior(ind1_new);
                
               n0 = length(ind0);
               n1 = length(ind1);
                
                if(n0>0)
                    beta_snr_i  = beta_snr_i  + 0.5 * inv_var_all*(W_i_snr -mue_prior_snr)' *inv(eye(n0)) *(W_i_snr-mue_prior_snr);
                end
                if(n1>0)
                    beta_coup_i = beta_coup_i + 0.5 * inv_var_all*(W_i_coup-mue_prior_coup)' *inv(eye(n1))*(W_i_coup-mue_prior_coup);
                end
            end
    end

    inv_lambda_coup           = gamrnd(alpha_coup_i,(1/beta_coup_i));
    lambda_coup_vec(i_node,1) = (1/inv_lambda_coup);
    
    inv_lambda_snr            = gamrnd(alpha_snr_i,(1/beta_snr_i));
    lambda_snr_vec(i_node,1)  = (1/inv_lambda_snr);
    
end


return
         
