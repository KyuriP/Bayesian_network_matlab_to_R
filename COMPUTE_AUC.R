
COMPUTE_AUC <- function(Run,TRUE,burn_in){
  # Compute edge scores:
  SCORES <- DIR(Run, burn_in)
  
  # Compute AUC (are under the precision recall curve)
  AUC <- AUPRC(TRUE, SCORES)
  
  return(AUC)
} 

###############################################################################

DIR <- function(Run, burn_in){
  
  DAGS = Run$dag
  
  n_samples = length(DAGS)
  n_nodes = dim(DAGS[[1]])[1]
  
  DIRECTED   = matrix(0, n_nodes,n_nodes)
  for (i in (burn_in+1):n_samples){
    DIRECTED = DIRECTED +  DAGS[[i]]
  }
  
  SCORES   = DIRECTED/(n_samples - burn_in)
  
  return(SCORES)  
  
}
  

###############################################################################

AUPRC <- function(True_Matrix, Aposteriori_Matrix){
  
  n_nodes = dim(True_Matrix)[1]
  
  for (i in 1:n_nodes){
    True_Matrix[i,i] = -1
  }
  
  n_edges = length(which(True_Matrix==1))
  
  Aposterioris = matrix()
  for (i in 1:n_nodes){
    for (j in 1:n_nodes){
      if (ifelse(i %in% j, 0, 1)){
        Aposterioris = matrix(c(Aposterioris, Aposteriori_Matrix[i,j]), nrow=1)
      }
    }
  }

  Aposterioris = sort(Aposterioris, decreasing = TRUE)
  
  # define APO_values first  ######################## is it appropriate?
  APO_values <- numeric(0)
  APO_values[1] = Aposterioris[1]
  
  for (i in 2:length(Aposterioris)){
    if (Aposterioris[i] == APO_values[length(APO_values)]){
    # do nothing 
    }
    else{
      APO_values[length(APO_values) + 1] = Aposterioris[i]
    }

  }
  
  MATRIX =  matrix(0, n_nodes, n_nodes)
  
  TPS = matrix()
  FPS = matrix()

  for (i in 1:length(APO_values)){
    indicis = which(Aposteriori_Matrix >= APO_values[i])
    MATRIX[indicis] = 1
    
    TP = length(which(MATRIX==1 & True_Matrix==1))
    FP = length(which(MATRIX==1 & True_Matrix==0))
    
    TPS = matrix(c(TPS,TP), nrow = 1)
    FPS = matrix(c(FPS,FP), nrow = 1)
    ###### remove NA from TPS and FPS   ##########
    TPS = TPS[!is.na(TPS)]
    FPS = FPS[!is.na(FPS)]
  }
  
  
  for (i in 2:length(TPS)){
    if ((TPS[i] - TPS[i-1])>1){
      
      NEW_TPS = matrix()
      NEW_FPS = matrix()
      
      for (x in 1:(TPS[i] - TPS[i-1]-1)){
        skew    = (FPS[i] - FPS[i-1])/(TPS[i] - TPS[i-1])
        NEW_TPS = matrix(c(NEW_TPS, TPS[i-1]+x), nrow = 1)
        NEW_FPS = matrix(c(NEW_FPS, FPS[i-1] + skew*x), nrow = 1)
      }
      TPS = matrix(c(TPS[1:i-1], NEW_TPS, TPS[i:length(TPS)]), nrow = 1)
      FPS = matrix(c(FPS[1:i-1], NEW_FPS, FPS[i:length(FPS)]), nrow = 1)
    }

  }
  
  PRECISION = TPS[2:length(TPS)] / (TPS[2:length(TPS)] + FPS[2:length(FPS)])
  RECALL    = TPS[2:length(TPS)] / n_edges
  
  PRECISION = matrix(c(1,PRECISION), nrow = 1)
  RECALL    = matrix(c(0,RECALL), nrow=1)
  
  auprc_com = pracma::trapz(RECALL,PRECISION)
  
  return(auprc_com)

} # function ends


  

  
  