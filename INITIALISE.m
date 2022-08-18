
function [Run] = INITIALISE(DATA_ALL, DAG, steps, nue_var, lambda_snr_vec, lambda_coup_vec, MATRIX, VECTORS)

[n_plus, m]= size(DATA_ALL{1});

n_nodes = length(DATA_ALL);

for node_i=1:n_nodes
    data = DATA_ALL{node_i}; 
    for component=1:max(MATRIX(node_i,:)) 
      DATA{node_i}{component} = data(:,find(MATRIX(node_i,:)==component)); 
    end  
end



[log_score] = COMPUTE_LOG_SCORE(DATA, DAG, MATRIX, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);

for i=1:(steps+1)
    Run.dag{i}             = 0;
    Run.matrix{i}          = 0;
    Run.Log_Scores(i)      = 0;
    Run.lambda_snr_vec{i}  = 0;
    Run.lambda_coup_vec{i} = 0;
    Run.VECTORS{i}         = 0;
end    
    
% Initialisation:
Run.dag{1}             = DAG;
Run.matrix{1}          = MATRIX;
Run.Log_Scores(1)      = log_score;
Run.steps(1)           = 1;
Run.lambda_snr_vec{1}  = lambda_snr_vec;
Run.lambda_coup_vec{1} = lambda_coup_vec;
Run.VECTORS{1}         = VECTORS;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [log_score] = COMPUTE_LOG_SCORE(DATA, DAG, MATRIX, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)

global Prior;

log_prob_breaks = 0;

[n_nodes, m]= size(MATRIX);

for i_node=1:n_nodes
 
	k = length(DATA{i_node});      % k is the number of mixture components 
	
	log_prob_k = log(poisspdf(k,1));

	k_cps = k-1;

	breakpoints = find(MATRIX(i_node,2:end)-MATRIX(i_node,1:end-1));
        
	if (length(breakpoints)==0)
   		log_prob_break = 0;
	else
            
    		breakpoints = [0,breakpoints,m];
            
    		log_prob_break = log(prod(1:(2*k_cps+1)))-log(prod(((m-1)-(2*k_cps+1)+1):(m-1)));
            
    		for i=2:length(breakpoints)
        		log_prob_break = log_prob_break + log(breakpoints(i)-breakpoints(i-1)-1);
    		end
	end

	log_prob_breaks = log_prob_breaks + log_prob_break + log_prob_k;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_prob_graph = 0;

for node=1:n_nodes
    log_prob_graph = log_prob_graph + Prior(length(find(DAG(:,node)))+1);  % this line ?
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_prob_data = 0;

for i_node=1:n_nodes
    
    k_i = length(DATA{i_node});
        
    parents = find(DAG(:,i_node));
    
    lambda_coup = lambda_coup_vec(i_node,1);
    lambda_snr  = lambda_snr_vec(i_node,1);
    
    sum_log_det_Sigma_tilde = 0;
    sum_Delta2              = 0;
    
    vector_i = VECTORS{i_node};
    
    ind1 = find(vector_i==1);
    ind0 = find(vector_i==0);
    
    LAMBDA_VEC       = vector_i;
    LAMBDA_VEC(ind0) = lambda_snr;
    LAMBDA_VEC(ind1) = lambda_coup;
    
    LAMBDA_MAT = diag(LAMBDA_VEC);
    LAMBDA_MAT = LAMBDA_MAT([1;parents+1],[1;parents+1]);    % how does this slicing work?
    
    
    %%% FOR THE FIRST SEGMENT:

    LAMBDA_VEC_first       = vector_i;
    LAMBDA_VEC_first(ind0) = lambda_snr; %%%
    LAMBDA_VEC_first(ind1) = lambda_snr; %%%

    LAMBDA_MAT_first = diag(LAMBDA_VEC_first);
    LAMBDA_MAT_first = LAMBDA_MAT_first([1;parents+1],[1;parents+1]);
    
    
    for component=1:k_i
        
        data = DATA{i_node}{component};
    
        [n_plus, n_obs] =  size(data);
            
        if(n_obs==0)
                % do nothing 
        else
             X = [ones(1,n_obs);data(parents,:)]; % pred x obs 
             y = data(end,:)';  % obs x 1

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
              
                m_tilde     = X'*mue_prior; % obs x 1
                    
                Sigma_tilde = eye(n_obs) + X'*LAMBDA*X;
                
                % pred x obs
                inv_Sigma_tilde = eye(n_obs) - X'*inv(inv(LAMBDA)+X*X')*X;
                
                                         %  (1 x obs) * (obs x obs)   * (obs x 1) 
                sum_Delta2 = sum_Delta2 + (y-m_tilde)'*inv_Sigma_tilde*(y-m_tilde); 
                
                sum_log_det_Sigma_tilde = sum_log_det_Sigma_tilde + log(det(Sigma_tilde));  
                                    
                Sigma_inv = inv(LAMBDA) + X*X';  % pred x pred 

                mue_apost  = inv(Sigma_inv)*(inv(LAMBDA)*mue_prior+X*y); % pred x 1
                           
        end
          
    end
    
    sum_1 = gammaln((m+nue_var)/2) - gammaln(nue_var/2);
                
    sum_2 = (nue_var/2)*log(nue_var) - (m/2)*log(pi) - 0.5 * sum_log_det_Sigma_tilde;
                
    sum_3 = -(m+nue_var)/2 * log(nue_var+sum_Delta2);
              
    log_score_i = sum_1 + sum_2 + sum_3;

    log_prob_data = log_prob_data + log_score_i;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global alpha_snr;
global beta_snr;
global alpha_coup;
global beta_coup;

log_prob_lambda = 0;

for i_node=1:n_nodes
    log_prob_lambda_snr_i  = log(gampdf(1/lambda_snr_vec(i_node,1), alpha_snr, (1/beta_snr)));
    log_prob_lambda_coup_i = log(gampdf(1/lambda_coup_vec(i_node,1),alpha_coup,(1/beta_coup)));
    log_prob_lambda        = log_prob_lambda + log_prob_lambda_snr_i + log_prob_lambda_coup_i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_prob_VECTOR = (sum(sum(DAG)) + n_nodes) * log(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_score = log_prob_breaks + log_prob_graph + log_prob_data + log_prob_lambda + log_prob_VECTOR;

return


