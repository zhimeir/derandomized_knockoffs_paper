vanilla_fwer_filter <- function(X,y, k = 1, alpha = 0.1,
                                knockoff_method = "gaussian",
                                knockoff_stat = stat.glmnet_coefdiff,seed = 24601,
                                mu = NULL,Sigma =NULL,#parameter for gaussian knockoff
                                pInit = NULL, Q = NULL,pEmit = NULL #parameter for hmm knockoff
                                ){
  ## Initialization
  set.seed(seed)
  n <- dim(X)[1]
  p <- dim(X)[2]

  ## Determine v
  v_list <- seq(1,20,by=1)
  p_list <- pnbinom(k,v_list,.5,lower.tail = FALSE)+dnbinom(k,v_list,.5) 
  if(sum(p_list<=alpha)>0){
    v0 <- max(which((p_list<=alpha) == TRUE))
    prob <- (alpha-p_list[v0])/(p_list[v0+1] - p_list[v0])
  }else{
    v0 <- 0
    prob <- (alpha)/(p_list[1])
  }

  ## Janson and Su procedure
  v <- floor(v0)+rbinom(1,1,prob)
  if(knockoff_method == "gaussian"){
    diags <- knockoff::create.solve_asdp(Sigma)
    Xk <- create.gaussian(X,mu,Sigma,diag_s = diags)}
  if(knockoff_method == "hmm"){
    Xk <- knockoffHMM(X, pInit, Q,pEmit,seed = seed+2020)
  }
  W <- knockoff_stat(X,Xk,y)
  order_w <- order(abs(W),decreasing = TRUE)
  sorted_w <- W[order_w]
  negid <- which(sorted_w<0)
  S <- c()
  ## Only consider v>0 because if v=0 then the selection set is empty
  if(v>0){
    ## if the total number of negitives is less than v, select all positives
    if(v>length(negid)){
      S <- which(W>0)
    }else{
      TT <- negid[v]
      S <- which(sorted_w[1:TT]>0)
      S <- order_w[S]
    }
  }
  return(list(S=S))
}


