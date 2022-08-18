

function [AUC] = COMPUTE_AUC(Run,TRUE,burn_in)

    % Compute edge scores:
    [SCORES] = DIR(Run, burn_in);
        
    % Compute AUC (are under the precision recall curve)
    [AUC] = AUPRC(TRUE, SCORES);
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SCORES] = DIR(Run, burn_in)

DAGS = Run.dag;

n_samples = length(DAGS);
n_nodes   = size(DAGS{1},1);

DIRECTED   = zeros(n_nodes,n_nodes);

for i = (burn_in+1):n_samples
    DIRECTED = DIRECTED +  DAGS{i};
end 

SCORES   = DIRECTED/(n_samples-burn_in);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [auprc_com] = AUPRC(True_Matrix, Aposteriori_Matrix)

n_nodes = length(True_Matrix);

for i=1:n_nodes
    True_Matrix(i,i)=-1;
end

n_edges = length(find(True_Matrix==1));

Aposterioris = [];

for i = 1:n_nodes
    for j = 1:n_nodes
        if(i~=j)
            Aposterioris = [Aposterioris, Aposteriori_Matrix(i,j)];
        end
    end
end

Aposterioris = sort(Aposterioris,'descend');

APO_values(1) = Aposterioris(1);

for i = 2:length(Aposterioris)
    if Aposterioris(i)==APO_values(end)
    % do nothing
    else
    APO_values(end+1) = Aposterioris(i);
    end
end


MATRIX = zeros(n_nodes,n_nodes);

TPS = [];
FPS = [];


for i = 1:length(APO_values)
    
    indicis = find(Aposteriori_Matrix>=APO_values(i));
    MATRIX(indicis) = 1;
    
    TP = length(find(MATRIX==1 & True_Matrix==1));
    FP = length(find(MATRIX==1 & True_Matrix==0));
    
    TPS = [TPS,TP];
    FPS = [FPS,FP];
end


for i=2:length(TPS)
    if (TPS(i)-TPS(i-1))>1
       
        NEW_TPS = [];
        NEW_FPS = [];
        
        for x = 1:(TPS(i)-TPS(i-1)-1)
            skew    = (FPS(i)-FPS(i-1))/(TPS(i)-TPS(i-1));
            NEW_TPS = [NEW_TPS,TPS(i-1)+x];
            NEW_FPS = [NEW_FPS,FPS(i-1)+ skew*x];
        end
        
        TPS = [TPS(1:i-1),NEW_TPS,TPS(i:end)];
        FPS = [FPS(1:i-1),NEW_FPS,FPS(i:end)];
        
    end
    
end
      
PRECISION = TPS(2:end)./(TPS(2:end)+FPS(2:end));
RECALL    = TPS(2:end)/n_edges;

PRECISION = [1,PRECISION];
RECALL    = [0,RECALL];


auprc_com = trapz(RECALL,PRECISION);

return;

