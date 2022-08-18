

function [Run] = PROC_YEAST(data_1, data_2, steps, step_iterations, k_max, k_transition, lambda_coup, lambda_snr, nue_var, MATRIX)


global DATA_ALL;

[DATA_1] = SHIFT_DATA(data_1);
[DATA_2] = SHIFT_DATA(data_2);

for i=1:length(DATA_1)
   DATA_ALL{i} = [DATA_1{i},DATA_2{i}];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Prior;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Prior] = SET_HYPERPARAMETERS(DATA_ALL);

n_nodes = length(DATA_ALL);

DAG = zeros(n_nodes,n_nodes);

lambda_coup_vec = lambda_coup * ones(n_nodes,1);
lambda_snr_vec  = lambda_snr  * ones(n_nodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_node=1:n_nodes
    VECTORS{i_node} = [0;-1 * ones(n_nodes,1)]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Run] = INITIALISE(DATA_ALL, DAG, steps, nue_var, lambda_snr_vec, lambda_coup_vec, MATRIX, VECTORS);
[Run] = START(DATA_ALL, steps, step_iterations, k_max, Run, k_transition, nue_var);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DATA_SHIFTED] = SHIFT_DATA(data)

[n, t] = size(data);

n_plus = n+1; % Covariance matrices will be of size (n+1)x(n+1)

for node_i=1:n
            % For each variable X_i consider
            % data_1[X_1(t-1),...,X_(n)(t-1),...,X_1(t-slice),...,X_(n)(t-slice), X_(i)(t)]      
            obs                  = 1:(t-1);  
            data_new             = zeros(n_plus,t-1);    
            data_new(1:n,obs)    = data(1:n,1:(t-1));       
            data_new(n+1,obs)    = data(node_i,2:t);   
            DATA_SHIFTED{node_i} = data_new;      
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Prior] = SET_HYPERPARAMETERS(DATA_ALL)

[n_plus, n_obs] = size(DATA_ALL{1});

n_nodes = n_plus-1;

Prior = zeros(n_nodes+1,1); 
 
return;

