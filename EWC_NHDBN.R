library(R.matlab)

rm(list=ls())

# Load the true network ('TRUE.mat') from the current working directory:
true <- readMat('TRUE.mat')[[1]]   ## NOTE: in Matlab it is named "TRUE" // in R, TRUE is a logical value which is not recommended to assign as a variable name

# check the data
head(true); class(true)
# TRUE(i,j) = 1 means that there is an edge from node i to node j
# TRUE(i,j) = 0 means that there is no edge from node i to node j

#################################################
# PREPARATION OF THE MCMC-BASED MODEL INFERENCE #
#################################################
# SELECT THE NUMBER OF MCMC ITERATIONS
# NOTE:
#In total: 'steps*step_iterations' 
# MCMC iterations will be performed 
# where 'step_iterations' is the thin-out factor

# E.g. for performing 5000 MCMC iterations
# and keeping only every 5-th MCMC sample (thinning out by the factor 5),
# set:
steps = 1000   
step_iterations = 5 
# Total number of MCMC iterations: steps*step_iterations = 5000
# And every 5-th sample will be kept.

##################################################################
# SELECT THE HYPERPARAMETERS OF THE EWC NH-DBN MODEL (CF. PAPER) #
##################################################################

#WE HAVE DEFINED SOME HYPERPARAMETERS AS GLOBAL VARIABLES, SO THAT
#THEY ARE ACCESSIBLE FROM ALL FUNCTIONS AND DO NOT ALWAYS HAVE 
#TO BE PROVIDED AS INPUT ARGUMENTS FOR THEM.

#THE HYPERPARAMETERS FOR THE TWO LAMBDA PARAMETERS    
# (LAMBDA_U AND LAMBDA_C) ARE GLOBAL

# global alpha_snr;   ## --> They are meant to be global variables
# global beta_snr;
# global alpha_coup;
# global beta_coup;


# SET HYPERPARAMETERS OF LAMBDA_U 
alpha_snr = 2 
beta_snr = 0.2
# SET HYPERPARAMETERS OF LAMBDA_C
alpha_coup = 2
beta_coup = 0.2

# HYPERPARAMATER FOR THE NOISE VARIANCE PARAMETER (NOT GLOBAL)
nue_var = 0.01                ####### nue_var is not a global variable!!


##############################################################
# LOAD THE TWO YEAST DATA FROM THE CURRENT WORKING DIRECTORY #
##############################################################
data_on_original <- readMat('data_on_original.mat')[[1]]
data_off_original <- readMat('data_off_original.mat')[[1]]

# check the data
head(data_on_original);class(data_on_original); dim(data_on_original)
head(data_off_original); class(data_off_original); dim(data_off_original)

# In both data sets 
# - the 5 rows refer to the 5 genes, as explained above
# - the columns refer to the temporal observations

################################################################
# STANDARDIZE THE YEAST DATA (CF. SUPPLEMENTARY PAPER, PART C) #
################################################################

#REMOVE WASHING PERIOD POINTS (= REMOVE FIRST TIME POINT FROM BOTH)
data_on  <- data_on_original[,-1]
data_off <- data_off_original[,-1]

# DETERMINE NUMBER OF NODES 'n' 
# AND DETERMINE THE NUMBERS OF REMAINING TIME POINTS 'm_1' and 'm_2'
n <- dim(data_on)[1] # == dim(data_off)[1]
m_1 <- dim(data_on)[2] # n=5, m_1 = 16-1 = 15
m_2 <- dim(data_off)[2] # n=5, m_2 = 21-1 = 20

# APPLY ZSCORE STANDARDIZATION ON MERGED DATA
data <- t(apply(cbind(data_on, data_off), 1, scale))

# THEN SEPARATE THE DATA AGAIN
data_on <- data[,1:m_1]
data_off <- data[,(m_1+1):35]

# THIS IS REQUIRED, AS WE WANT TO SHIFT BOTH DATA SETS SEPARATELY.
# AFTER THE SHIFTS THE SHIFTED DATA WILL BE MERGED AGAIN
# THIS HAPPENS IN THE FUNCTION 'PROC_YEAST.m'
# FOR DETAILS PLEASE SEE SUPPLEMENTARY PAPER, PART C

###################################################################### 
# SPECIFY THE INITIAL DATA SEGMENTATION THROUGH AN ALLOCATION MATRIX #
######################################################################

# For all 'n' genes we initialize with H=1 (no changepoints).
MATRIX = matrix(1, nrow = n, ncol = m_1 + m_2 - 2)

# Note that we will lose two more observations. 
# For both data sets we lose one data point, 
# as we will shift both data sets separately.
# (Each edge interaction is subject to a time lag of 1.)

# HOW TO INTERPRET MATRIX?
# MATRIX(i,j) = k
# MEANS THAT THE 'j'-TH DATA POINT OF TARGET GENE 'i' BELONGS TO SEGMENT 'h'

# E.G. IF WE WANTED TO INITIALISE WITH THE TRUE DATA SEGMENTATION
# MATRIX = ones(n,m_1+m_2-2);
# MATRIX(:,m_1:end) = 2;
# Then for all genes (=all rows)
# - the first (m_1-1) data points belong to segment h=1
# - the last  (m_2-1) data points belong to segment h=2

###################################################################### 

# BEFORE RUNNING AN MCMC SIMULATION 
# SPECIFY SOME TUNING PARAMETERS
# AND INITIALIZE THE TWO LAMBDA PARAMETERS 
# Specify the maximal number of data segments per gene:
H_max = 10

# Because of the changepoint location prior, neighbouring changepoints
# must have at least the distance (\tau_{h+1}-tau_{h}= 2). 
# Give this piece of information to the MCMC algorithm: 
k_transition = 2

# Initial lambda parameters
lambda_snr  = 1 # corresponds to $\lambda_{u}$ in the paper
lambda_coup = 1 # corresponds to $\lambda_{c}$ in the paper

#################################
### START THE MCMC SIMULATION ###
#################################

# Depending on the selected number of MCMC iterations 
# (selected via 'steps' and 'step_iterations', as explained above)
# the MCMC simulation might take some time.

[Run] = PROC_YEAST(data_on, data_off, steps, step_iterations, H_max, k_transition, lambda_snr, lambda_coup, nue_var, MATRIX);

################################# 
### THE OUTPUT FILE 'Run.mat' ### 
################################# 

# Most important are the sampled network structures.
# From the sampled network structures, the marginal edge posterior 
# probabilities (=edge scores) can be computed.

# Number of sampled network structures:
n_samples = length(Run.dag) # n_samples = steps + 1;

# The MCMC sampled networks are
# Run.dag{2},...,Run.dag{n_samples}

# Note that Run.dag{1} only contains the initial network structure.
# Default initialisation: A network without any edges.

# Average the sampled network structures to get estimates of the 
# marginal edge posterior probabilities (= edge scores)

# Initialize:
NET = zeros(n,n)

# When ignoring the first 100 network structures 
# to take a burn-in phase into account: 

burn_in = 100 # 'burn_in' must be lower than 'steps'

# Add all networks matrices (sampled after burn-in) up:
for i_sample=(burn_in+1):n_samples
NET = NET + Run.dag{i_sample};
end

NET = NET/(n_samples - burn_in)

# Consider the matrix of edge scores:
NET

# DAG(i,j) is the score of the edge i->j

# E.g. the score of the edge 1->2 is
NET(1,2)
# In this application, it is the score for the edge CBF1 -> GAL4 


###################################################
### Network reconstruction accuracy (AUC value) ### 
###################################################
# Load the true network:
load('TRUE');

# And compute the AUC value:
[AUC] = COMPUTE_AUC(Run,TRUE,burn_in);

# The precision-recall AUC for this MCMC simulations is:
AUC


####################################################################### 
### Note that the 'Run' file also contains other interesting output ###
### Here we give a few example diagnostics                          ###
####################################################################### 

# [1] Make a plot of the logarithmic scores (log_likelihood + log_priors)
# along the MCMC iterations:

figure(1)
clf
plot(Run.Log_Scores)

###################################################################### 

# [2] The 10-th network structure in Run is

j = 10;

Run.dag{j}
# The corresponding data point segmentation (allocation) is:
Run.matrix{j}

###################################################################### 

# [3] Identifying coupled/uncoupled edges:
# To get insight, which edges of the j-th sampled network 
# were coupled/uncoupled, extract the information as follows:

# e.g. if j=10
j = 10;

OUT = [];

for i_node = 1:n
OUT = [OUT,Run.VECTORS{j}{i_node}];
end

# Remore the first row (as the first row corresponds to the intercept)
OUT_EDGES_ONLY = OUT(2:end,:);
OUT_EDGES_ONLY

# Compare with the j-th DAG: Run.dag{j}

Run.dag{j}


# OUT_EDGES(i,j) = -1 means that there is           no edge from i to j
# OUT_EDGES(i,j) =  0 means that there is an uncoupled edge from i to j
# OUT_EDGES(i,j) = +1 means that there is a coupled    edge from i to j

# For DAG(i,j) = 0 we have: OUT_EDGES(i,j) = -1 (no edge Z_i->Z_j)
# For DAG(i,j) = 1 we have: either OUT_EDGES(i,j) = 1 or OUT_EDGES(i,j) = 0
# There is an edge i->j and it is either coupled (1) or uncoupled (0).

###################################################################### 

# [4] To get the fractions of coupled/uncoupled edges, proceed as follows:

# The edge score for the edge 1 -> 2 is
parent_node = 1;
child_node  = 2;
edge_score = NET(parent_node,child_node);

edge_score

# Count the number of coupled and uncoupled edges

n_coupled   = 0;
n_uncoupled = 0;

for i_sample=(burn_in+1):n_samples

edge_status = Run.VECTORS{i_sample}{child_node}(1+parent_node); 
# add 1 to parent_node to jump over the first element for the intercept  

if(edge_status==0) % edge was uncoupled
n_uncoupled = n_uncoupled + 1; 
elseif(edge_status==1) % edge was coupled
n_coupled = n_coupled + 1;
else % edge was not even present
% so do nothing.
end
end

# Compute the fractions:
p_coupled   = n_coupled/(n_samples-burn_in);
p_uncoupled = n_uncoupled/(n_samples-burn_in);

% Note that
p_coupled + p_uncoupled
% is equal to
edge_score

% With which fractions was the edge 'parent_node -> child_node' 
% 'coupled' and 'uncoupled'? 
  
  % Simply compute the ratios:
  
  [p_coupled,p_uncoupled]/edge_score

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [5] Make trace plots of the sampled lambda parameters

% For example, let's make the trace plots for the third target node (=SWI5)

% Set:
node = 3;

% Collect the sampled parameters in:
lambda_u = []; % uncoupled
lambda_c = []; % coupled

% Let's here include the burn-in samples
for i_sample=1:n_samples % without burn_in:  for i_sample=(burn_in+1):n_samples
lambda_u = [lambda_u,Run.lambda_snr_vec{i_sample}(node)];
lambda_c = [lambda_c,Run.lambda_coup_vec{i_sample}(node)];
end

% Make trace plots
figure(2)
clf
subplot(2,1,1)
plot(lambda_u)
subplot(2,1,2)
plot(lambda_c)

% END OF EXAMPLE DIAGNOSTICS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW TO USE THE CODE WHEN THERE IS ONLY ONE SINGLE TIME SERIES?
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If there is only one single time series 'data' available, 
% then the function 'PROC_ONLY_ONE.m' 
% must be used instead of 'PROC_YEAST.m'

% For example:
  % Generate a random data set with 'n=5' nodes and 'm=21' time points:
  data = randn(5,21);

% For demonstrating purposes, let's introduce 
% a very strong uncoupled edge: '1->2'

data(2,2:11)  = (+2) * data(1,1:10)  + 0.01 * randn(1,10);
data(2,12:21) = (-2) * data(1,11:20) + 0.01 * randn(1,10);

% Alternatively, you can here use your own 'data' matrix instead

% E.g. via:
% data = ...
% load('data')

% From now on 'data' is a Matlab matrix.
% Every row corresponds to a network variable (e.g. a gene).
% Every column corresponds to a time point.


[n,m] = size(data);
% There are: 'n' network variables and 'm' time points.

% Initialize the data point segmentation matrix so that there are
% no changepoints (i.e. H=1 for all variables).
MATRIX = ones(n,m-1); 
% because of the time lag, we have to shift the data
% so that we lose one data point (m -> m-1)

% Start an MCMC simulation with 'steps*step_iterations' MCMC iterations
steps = 1000;
step_iterations = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And as before, we again set:
H_max = 10;
k_transition = 2;
lambda_snr = 1;
lambda_coup = 1;
nue_var = 0.01;
%global alpha_snr;
%global beta_snr;
%global alpha_coup;
%global beta_coup;
alpha_snr  = 2;
beta_snr   = 0.2;
alpha_coup = 2;
beta_coup  = 0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start the MCMC simulation"
% The 'data' will be shifted within the function 'PROC_ONLY_ONE.m':
[Run] = PROC_ONLY_ONE(data, steps, step_iterations, H_max, k_transition, lambda_snr, lambda_coup, nue_var, MATRIX);

% As before, compute the edge scores from the 'Run' file:

n_samples = length(Run.dag);

burn_in = 100;

NET = zeros(n,n);

for i_sample=(burn_in+1):n_samples
    NET = NET + Run.dag{i_sample};
end

% Compute the edge scores:
NET = NET/(n_samples-burn_in);

NET

% Has EWC NH-DBN  'found' the true edge (1->2)?
% The edge score is hopefully near 1:
NET(1,2)

% Also the other diagnostics [1-5] can be applied, as described above.
% For example:

% Make a trace plot of logarithmic scores:
figure(3)
clf
plot(Run.Log_Scores)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Has EWC NH-DBN inferred the changepoint?

% Let's have a look at the sampled allocation vectors for target node 2.

% Collect all allocation vectors for target node 2

node = 2;

ALLOCATION = [];

for i_sample=1:n_samples
ALLOCATION = [ALLOCATION;Run.matrix{i_sample}(node,:)];
end

% Every row of 'ALLOCATION' is a sampled allocation vector


% What are the probabilities of a changepoint a the different locations? 
  
  % Hopefully there is a high probability at location 10 
% Note that data point 10 refers to the original time point 11, 
% because of the time lag (data shift).

% Posterior probabilities per location:
  sum(ALLOCATION(:,2:end) - ALLOCATION(:,1:end-1))/n_samples

% Make a plot of those probabilities:
  figure(4)
clf
plot(sum(ALLOCATION(:,2:end) - ALLOCATION(:,1:end-1))/n_samples)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SCRIPT FILE 
% BY MARCO GRZEGORCZYK 
% BERNOULLI INSTITUTE, GRONINGEN UNIVERSITY, NL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

