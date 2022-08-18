

function [Run] = START(DATA_ALL, steps, step_iterations, k_max, Run, k_transition, nue_var)

fan_in = 3;                           % maximal size of parent-sets                      

DAG             = Run.dag{1};         % the current DAG
MATRIX          = Run.matrix{1};      % the current allocation matrix
log_score       = Run.Log_Scores(1);  % the current log score
lambda_snr_vec  = Run.lambda_snr_vec{1};
lambda_coup_vec = Run.lambda_coup_vec{1};
VECTORS         = Run.VECTORS{1};

fprintf('\n###########################################################\n')
fprintf('An MCMC simulation for the EWC NH-DBN model has been started \n')
fprintf('###########################################################\n')
fprintf('Log-score of the initial graph: %1.5f \n',log_score)
fprintf('###########################################################\n\n')

%%% Start of the MCMC simulation

Counter = 2;

for i = 1:steps

[DAG, log_score, MATRIX, lambda_snr_vec, lambda_coup_vec, nue_var, VECTORS] = MCMC_INNER(step_iterations, DAG, fan_in, log_score, MATRIX, DATA_ALL, k_max, nue_var, lambda_snr_vec, lambda_coup_vec, k_transition, VECTORS);

    Run.dag{Counter}              = DAG;
    Run.matrix{Counter}           = MATRIX;
    Run.Log_Scores(Counter)       = log_score;
    Run.steps(1)                  = i+1;
    Run.lambda_snr_vec{Counter}   = lambda_snr_vec;
    Run.lambda_coup_vec{Counter}  = lambda_coup_vec;
    Run.VECTORS{Counter}          = VECTORS;
    Counter = Counter+1;

end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DAG, log_score, MATRIX, lambda_snr_vec, lambda_coup_vec, nue_var, VECTORS] = MCMC_INNER(step_iterations, DAG, fan_in, log_score, MATRIX, DATA_ALL, k_max, nue_var, lambda_snr_vec, lambda_coup_vec, k_transition, VECTORS)

[n_plus, m]  = size(DATA_ALL{1}); 
[n_parents_keep, n_nodes] = size(DAG);

for t=1:step_iterations
        
for i_node=1:n_nodes

    x = rand;

    DATA     = [];
    DATA_NEW = [];
    
    p_0 = 0.5; % structure MCMC
    p_1 = 0.7; % Birth Move
    p_2 = 0.9; % Death Move
    % else     % Reallocation Move
    
    if(x<=p_0) % PERFORM A STRUCTURE-MCMC MOVE
        
       y = rand(1);
      
       if (y<0.5) % Perform an addition/deletion move
        
            child_node    = i_node;
            
            old_parents = DAG(:,child_node);
            new_parents = old_parents;
            
            parents_card_old = sum(old_parents);
                
            parents_indicis_old = find(old_parents);
                          
                if (parents_card_old<fan_in)
            
                        neighbours_old = n_parents_keep-1;
                        
                        indicis = randperm(n_parents_keep);
                        x_ind = indicis(1);
                        
                        if (x_ind==child_node)
                            x_ind = indicis(2);
                        end

                        parent_index = x_ind; % delete or add this parent node
            
                        new_parents(parent_index,1) = 1 - new_parents(parent_index,1); 
                
                else % elseif (parent_card_old==fan_in)
                             
                    x_ind = randperm(fan_in);
                    x_ind = x_ind(1);
                
                    parent_index = parents_indicis_old(x_ind); % delete this parent node
                
                    new_parents(parent_index,1) = 0;
                
                    neighbours_old = parents_card_old; % = fan_in
                    end
                
                    parents_card_new    = sum(new_parents);
                
                    if (parents_card_new<fan_in)         
                        neighbours_new = n_parents_keep-1; 
                    else
                        neighbours_new = parents_card_new; % = fan_in
                    end
             
                    DAG_NEW = DAG;
                    DAG_NEW(:,child_node) = new_parents;
                
            
                        data = DATA_ALL{child_node};
                        k_i = max(MATRIX(child_node,:));
                        for component=1:k_i      
                            DATA{child_node}{component} = data(:,find(MATRIX(child_node,:)==component));
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        VECTORS_NEW  = VECTORS;   
                        vector_i_new = VECTORS_NEW{child_node};
                        
                       if (new_parents(parent_index,1)==1) % proposal is to add a parent
                                           
                                vector_i_new(parent_index+1,1) = (rand(1,1)<0.5);
                                log_HR_supplement = log(2);
                       else % proposal is to delete a parent
                                vector_i_new(parent_index+1,1) = -1;
                                log_HR_supplement = -log(2);
                       end
                       
                       VECTORS_NEW{child_node} = vector_i_new;
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
                local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG_NEW, MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW);
                local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG,     MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);
           
                A = exp(local_score_new - local_score_old + log(neighbours_old) - log(neighbours_new) + log_HR_supplement);

                u = rand(1); 
 
                if (u<A) % accept the move:
                    DAG       = DAG_NEW; 
                    VECTORS   = VECTORS_NEW;

                    log_score = log_score + local_score_new - local_score_old;
                end
         
                clear DATA;
                
       else % Perform an exchange move
              
                child_node    = i_node;

                old_parents = find(DAG(:,child_node));
                
                parents_card_old = length(old_parents);
                
                    if (parents_card_old==0)
            
                        % do nothing
                
                    else % perform an exchange move
                             
                        DAG_NEW = DAG;
                        
                        indicis        = randperm(parents_card_old);
                        index          = indicis(1);
                
                        parent_old_index = old_parents(index); % delete this parent node
                        
                        candidate_parents = find(DAG(:,child_node)==0);
                        
                        candidates_card = length(candidate_parents);
                        
                        indicis = randperm(candidates_card);
                        index   = indicis(1);
                        
                        parent_new_index = candidate_parents(index);
                        
                        
                        if (parent_new_index==child_node)
                            index = indicis(2);
                            parent_new_index = candidate_parents(index);
                        end
                        
                        DAG_NEW(parent_old_index,child_node) = 0;
                        DAG_NEW(parent_new_index,child_node) = 1;
                
                            data = DATA_ALL{child_node};
                            k_i = max(MATRIX(child_node,:));
                            
                            for component=1:k_i      
                                DATA{child_node}{component} = data(:,find(MATRIX(child_node,:)==component));
                            end
                      
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        VECTORS_NEW = VECTORS;
                        vector_i_new = VECTORS_NEW{child_node};
                           
                        vector_i_new(parent_old_index+1,1) = -1;
                        vector_i_new(parent_new_index+1,1) = (rand(1,1)<0.5);
                         
                        VECTORS_NEW{child_node} = vector_i_new;
                            
                        log_HR_supplement = 0;
                           
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG_NEW, MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW);
                        local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG,     MATRIX, child_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);
            
                        log_hastings = 0;
                        
                        A = exp(local_score_new - local_score_old + log_hastings + log_HR_supplement);
        
                        u = rand(1); 
 
                        if (u<A) % accept the move:
                            DAG       = DAG_NEW; 
                            VECTORS   = VECTORS_NEW;
                            log_score = log_score + local_score_new - local_score_old;  
                        end
         
                        clear DATA;
                    end
              
       end
                    
    elseif(x<=p_1) % PERFORM A BREAKPOINT BIRTH MOVE
       
                if(max(MATRIX(i_node,:))<k_max)
                
                NEW_CANDIDATES = ones(1,m-1);

                break_vec     = find(MATRIX(i_node,2:end)-MATRIX(i_node,1:(end-1)));
                break_vec_add = [0,break_vec,m];
                
                    for i=break_vec_add     %location for the new candidates
                        NEW_CANDIDATES(1,max([i-(k_transition-1),1]):min([i+(k_transition-1),m-1])) = 0;
                    end
            
                NEW_CANDIDATES = find(NEW_CANDIDATES==1);
          
                n_candidates = length(NEW_CANDIDATES);
               
                if (n_candidates>0)
            
                    indicis      = randperm(n_candidates);
                    index        = indicis(1);
                    
                    index        = NEW_CANDIDATES(index);
                
                    MATRIX_NEW = MATRIX;
                    MATRIX_NEW(i_node,(index+1):end) = MATRIX_NEW(i_node,(index+1):end)+1;
                
                    clear DATA;
                    clear DATA_NEW;

                    data = DATA_ALL{i_node};
                    
                    for component=1:max(MATRIX(i_node,:))   
                        DATA{i_node}{component} = data(:,find(MATRIX(i_node,:)==component));
                    end
                  
                    for component=1:max(MATRIX_NEW(i_node,:)) 
                        DATA_NEW{i_node}{component} = data(:,find(MATRIX_NEW(i_node,:)==component));
                    end
                   
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    vector_i_new                  = VECTORS{i_node};
                    
                    parents   = find(DAG(:,i_node));
                    n_parents = length(parents);
                    
                    vector_i_new([1;parents+1],1) = (rand(n_parents+1,1)<0.5);
                 
                    log_HR_supplement  = 0;
                      
                    VECTORS_NEW         = VECTORS;
                    VECTORS_NEW{i_node} = vector_i_new; 

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
                    local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_NEW, DAG, MATRIX_NEW, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW);
                    local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA    , DAG, MATRIX,     i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);
             
                    n_breakpoints_new = max(MATRIX_NEW(i_node,:));
                
                    log_hastings = log_HR_supplement + log(n_candidates)-log(n_breakpoints_new);
                
                    A = exp(local_score_new - local_score_old + log_hastings);
                
                    u = rand(1); 
              
                    if (u<A) % accept the move:
                        MATRIX    = MATRIX_NEW;
                        VECTORS   = VECTORS_NEW;
                        log_score = log_score + local_score_new - local_score_old;
                    end
                  
                end
                
            end
                
            clear DATA;
            clear DATA_NEW;
            

    elseif(x<=p_2)  % PERFORM A BREAKPOINT DEATH MOVE
      
          CANDIDATES = find(MATRIX(i_node,2:end)-MATRIX(i_node,1:end-1));
         
            if(length(CANDIDATES)==0) % then there is no breakpoint which can be removed
             
            else
                        
               n_candidates = length(CANDIDATES);
               indicis = randperm(n_candidates);
               index   = indicis(1);
               
               candidate = CANDIDATES(index);
               
               MATRIX_NEW                           = MATRIX;
               MATRIX_NEW(i_node,(candidate+1):end) = MATRIX_NEW(i_node,(candidate+1):end)-1;
               
               clear DATA;
               clear DATA_NEW;
               
               data = DATA_ALL{i_node};
               
               for component=1:max(MATRIX(i_node,:)) 
                    DATA{i_node}{component} = data(:,find(MATRIX(i_node,:)==component));
               end
               
               for component=1:max(MATRIX_NEW(i_node,:)) 
                    DATA_NEW{i_node}{component} = data(:,find(MATRIX_NEW(i_node,:)==component));
               end
             
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    vector_i_new                  = VECTORS{i_node};
                    
                    parents   = find(DAG(:,i_node));
                    n_parents = length(parents);
                    
                    vector_i_new([1;parents+1],1) = (rand(n_parents+1,1)<0.5);
                 
                    log_HR_supplement  = 0;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    VECTORS_NEW = VECTORS;
                    VECTORS_NEW{i_node} = vector_i_new; 
                  
             
               local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_NEW, DAG, MATRIX_NEW, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW);
               local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA    , DAG, MATRIX,     i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);
              
                
               BIRTH_CANDIDATES = ones(1,m-1);
               break_vec_new = find(MATRIX_NEW(i_node,2:end)-MATRIX_NEW(i_node,1:end-1));
            
               break_vec_new_add = [0,break_vec_new,m];
                
                for i=break_vec_new_add
                    BIRTH_CANDIDATES(1,max([i-(k_transition-1),1]):min([i+(k_transition-1),m-1]))=0;
                end
             
                BIRTH_CANDIDATES   = find(BIRTH_CANDIDATES==1);
                n_birth_candidates = length(BIRTH_CANDIDATES);
                
                log_hastings =  log_HR_supplement + log(n_candidates) - log(n_birth_candidates); 
                
                A = exp(local_score_new - local_score_old + log_hastings);

                u = rand(1); 
 
                if (u<A) % accept the move:
                    MATRIX      = MATRIX_NEW; 
                    VECTORS     = VECTORS_NEW;
                    log_score   = log_score + local_score_new - local_score_old;     
                end
   
            end
         
            clear DATA;
            clear DATA_NEW;
 
             
    else % PERFORM A BREAKPOINT REALLOCATION MOVE
         
          CANDIDATES   = find(MATRIX(i_node,2:end)-MATRIX(i_node,1:end-1)); 
          n_candidates = length(CANDIDATES);
        
          if(n_candidates>0) % then there are breakpoints which can be re-allocated
           
                indicis   = randperm(n_candidates);
                index_old = indicis(1);
                candidate = CANDIDATES(index_old);
                
                CANDIDATES_ADD = [0,CANDIDATES,m];
                
                index_candidate = find(CANDIDATES_ADD==candidate);
                
                k_minus = CANDIDATES_ADD(index_candidate-1);
                k_i     = CANDIDATES_ADD(index_candidate); 
                k_plus  = CANDIDATES_ADD(index_candidate+1);    
                    
                k_i_new_possible = (k_minus+k_transition):(k_plus-k_transition);
                n_candidates_new = length(k_i_new_possible);
                
                indicis_new = randperm(n_candidates_new);
                index_new   = indicis_new(1);
                
                candidate_new = k_i_new_possible(index_new);
                
                MATRIX_NEW                                = MATRIX;
                MATRIX_NEW(i_node, (candidate+1):end)     = MATRIX_NEW(i_node,(candidate+1):end) - 1;
                MATRIX_NEW(i_node, (candidate_new+1):end) = MATRIX_NEW(i_node,(candidate_new+1):end) + 1;
                
                clear DATA;
                clear DATA_NEW;
                
 
                data = DATA_ALL{i_node};
                
                for component=1:max(MATRIX(i_node,:))   
                    DATA{i_node}{component} = data(:,find(MATRIX(i_node,:)==component));
                end
          
                for component=1:max(MATRIX_NEW(i_node,:))  
                    DATA_NEW{i_node}{component} = data(:,find(MATRIX_NEW(i_node,:)==component));
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                vector_i_new                  = VECTORS{i_node};
                    
                parents   = find(DAG(:,i_node));
                n_parents = length(parents);
                    
                vector_i_new([1;parents+1],1) = (rand(n_parents+1,1)<0.5);
                 
                log_HR_supplement  = 0;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                VECTORS_NEW         = VECTORS;
                VECTORS_NEW{i_node} = vector_i_new;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_NEW, DAG, MATRIX_NEW, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW);
                local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA    , DAG, MATRIX,     i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);

                A = exp(local_score_new - local_score_old + log_HR_supplement);

                u = rand(1); 
 
                        if (u<A) % accept the move:
                        MATRIX      = MATRIX_NEW; 
                        VECTORS     = VECTORS_NEW;
                        log_score   = log_score + local_score_new - local_score_old; 
                        end
                    
                    end
                
            end
            
            clear DATA;
            clear DATA_NEW;
                
    end % MOVE-TYPE-LOOP
     
% Update hyperparameters:

clear DATA;
DATA = [];
       
for node_i=1:n_nodes
    data = DATA_ALL{node_i}; 
    for component=1:max(MATRIX(node_i,:)) 
        DATA{node_i}{component} = data(:,find(MATRIX(node_i,:)==component));
    end    
end
        
[lambda_snr_vec, lambda_coup_vec, VECTORS]  = UPDATE(DATA, DAG, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);

% ADDITIONAL VECTOR_MOVES

for i_node=1:n_nodes
   
    VECTORS_NEW = VECTORS;
    vector_i    = VECTORS{i_node};
    
    parents   = find(DAG(:,i_node));
    n_parents = length(parents); 
    
    indicis = randperm(n_parents+1);
    index   = indicis(1);
    
    coefficients = [1;parents+1];
    coefficient  = coefficients(index);
    
    vector_i_new = vector_i;
    vector_i_new(coefficient,1) = 1-vector_i_new(coefficient,1);
    
    VECTORS_NEW{i_node} = vector_i_new;
        
    local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS_NEW);
    local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);
              
    log_hastings = 0;
                
    A = exp(local_score_new - local_score_old + log_hastings);
    u = rand(1); 

                if (u<A) % accept the move:
                    VECTORS   = VECTORS_NEW;
                    log_score = log_score + local_score_new - local_score_old;     
                end     
end
           

[log_score] = COMPUTE_LOG_SCORE(DATA, DAG, MATRIX, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS);
     
  
end % ITERATION-STEPS-LOOP

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [log_score] = COMPUTE_LOCAL_LOG_SCORE(DATA, DAG, MATRIX, i_node, nue_var, lambda_snr_vec, lambda_coup_vec, VECTORS)

global Prior;

[n_nodes, m]= size(MATRIX);

k = length(DATA{i_node});      
		
log_prob_k = log(poisspdf(k,1));

k_cps = k-1;
	
breakpoints = find(MATRIX(i_node,2:end)-MATRIX(i_node,1:end-1));
        	
if (length(breakpoints)==0)
   	log_prob_break = 0;
	
else
                	
	breakpoints    = [0,breakpoints,m];
	log_prob_break = log(prod(1:(2*k_cps+1)))-log(prod(((m-1)-(2*k_cps+1)+1):(m-1)));
            
	for i=2:length(breakpoints)
		log_prob_break = log_prob_break + log(breakpoints(i)-breakpoints(i-1)-1);
	end

end

log_prob_breaks = log_prob_break + log_prob_k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_prob_graph =  Prior(length(find(DAG(:,i_node)))+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the local score for i_node:

k_i = length(DATA{i_node}); 
    
parents   = find(DAG(:,i_node));
n_parents = length(parents);

sum_log_det_Sigma_tilde = 0;
sum_Delta2              = 0;

vector_i = VECTORS{i_node};

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

for component=1:k_i
        
        data = DATA{i_node}{component};
    
        [n_plus, n_obs] =  size(data);
            
        if(n_obs==0)
                % do nothing 
        else   
            X = [ones(1,n_obs);data(parents,:)]; % pred x obs 
            y = data(end,:)'; % obs x 1
            
            if (component==1)
                mue_prior = zeros(length(parents)+1,1); % pred x 1  
                LAMBDA    = LAMBDA_MAT_first;
            else
                
                if (length(parents)>0)
                    mue_prior = vector_i([1;parents+1],1) .* mue_apost;
                else % only intercept
                    mue_prior = vector_i(1,1)             .* mue_apost;
                end
                
                LAMBDA = LAMBDA_MAT;
            end
            
             m_tilde     = X'*mue_prior; % obsx1
                    
             Sigma_tilde = eye(n_obs) + X'*LAMBDA*X;
                
             % pred x obs
             inv_Sigma_tilde = eye(n_obs) - X'*inv(inv(LAMBDA)+X*X')*X;
                
             % (1xobs) * (obsxobs) * (obsx1) 
             sum_Delta2 = sum_Delta2 + (y-m_tilde)'*inv_Sigma_tilde*(y-m_tilde); 
                
             sum_log_det_Sigma_tilde = sum_log_det_Sigma_tilde + log(det(Sigma_tilde));  
                                 
             % pred x pred              
             Sigma_inv = inv(LAMBDA) + X*X';
                
             % pred x 1 
             mue_apost  = inv(Sigma_inv)*(inv(LAMBDA)*mue_prior+X*y);
                         
        end
        
end
    
sum_1 = gammaln((m+nue_var)/2) - gammaln(nue_var/2);
                
sum_2 = (nue_var/2)*log(nue_var) - (m/2)*log(pi) - 0.5 * sum_log_det_Sigma_tilde;
                
sum_3 = -(m+nue_var)/2 * log(nue_var+sum_Delta2);
              
log_score_i = sum_1 + sum_2 + sum_3;


%%%%%%%%%%%%%%%%%%%%%%%%%

log_prob_VECTOR = (n_parents + 1) * log(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%

log_score = log_prob_breaks + log_prob_graph + log_score_i + log_prob_VECTOR;

return

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
    log_prob_graph = log_prob_graph + Prior(length(find(DAG(:,node)))+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_prob_data = 0;

for i_node=1:n_nodes
    
    k_i = length(DATA{i_node}); % k_i is the number of mixture components for the i-th node
        
    parents = find(DAG(:,i_node));
    
    sum_log_det_Sigma_tilde = 0;
    sum_Delta2              = 0;
     
    vector_i = VECTORS{i_node};

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

    for component=1:k_i
        
        data = DATA{i_node}{component};
    
        [n_plus, n_obs] =  size(data);
            
        if(n_obs==0)
                % do nothing 
        else
            
             X = [ones(1,n_obs);data(parents,:)]; % pred x obs 
             y = data(end,:)'; % obsx1
 
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
             
                m_tilde = X'*mue_prior; % obsx1
                    
                Sigma_tilde = eye(n_obs) + X'*LAMBDA*X;
                
                % pred x obs
                inv_Sigma_tilde = eye(n_obs) - X'*inv(inv(LAMBDA)+X*X')*X;
                
                % (1xobs) * (obsxobs) * (obsx1) 
                sum_Delta2 = sum_Delta2 + (y-m_tilde)'*inv_Sigma_tilde*(y-m_tilde); 
                
                sum_log_det_Sigma_tilde = sum_log_det_Sigma_tilde + log(det(Sigma_tilde));  
                      
                % pred x pred              
                Sigma_inv = inv(LAMBDA) + X*X';

                % pred x 1
                mue_apost  = inv(Sigma_inv)*(inv(LAMBDA)*mue_prior+X*y);
                
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
    log_prob_lambda        = log_prob_lambda+ log_prob_lambda_snr_i + log_prob_lambda_coup_i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_prob_VECTOR = (sum(sum(DAG)) + n_nodes) * log(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

log_score = log_prob_breaks + log_prob_graph + log_prob_data + log_prob_lambda + log_prob_VECTOR;

return

