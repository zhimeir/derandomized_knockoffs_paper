vanilla_pfer_filter <- function(X,y, v0 = 1, knockoff_method = "gaussian",
                        knockoff_stat = stat.glmnet_coefdiff,seed = 24601,
                        mu = NULL,Sigma =NULL,#parameter for gaussian knockoff
                        pInit = NULL, Q = NULL,pEmit = NULL #parameter for hmm knockoff
){
  #check input
  #Initialization
  set.seed(seed)
  n = dim(X)[1]
  p = dim(X)[2]
  pi = rep(0,p)
  lambda = rep(0,p)
  if(knockoff_method == "gaussian"){
    diags = knockoff::create.solve_asdp(Sigma)
  }
  
  #Main part
  v = floor(v0)+rbinom(1,1,v0-floor(v0))
  if(knockoff_method == "gaussian"){
    Xk = create.gaussian(X,mu,Sigma,diag_s = diags)}
  if(knockoff_method == "hmm"){
    Xk = knockoffHMM(X, pInit, Q,pEmit,seed = seed)
  }
  W = knockoff_stat(X,Xk,y)
  order_w = order(abs(W),decreasing = TRUE)
  sorted_w = W[order_w]
  negid = which(sorted_w<0)
  S = c()
  if(v>0){ #if v=0 then the selection set is empty
    if(v> length(negid)){ #if the total number of negitives is less than v, select all positives
      S = which(W>0)
    }else{
      TT = negid[v]
      S = which(sorted_w[1:TT]>0)
      S = order_w[S]
    }
  }
  
return(list(S=S))
}


