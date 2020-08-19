filter_ukb <- function(W0,k,alpha = 0.1,tau =0.5,seed = 24601){
#Initialization
  M = dim(W0)[1]
  p = dim(W0)[2]
  alpha = alpha  #account for multiple k
  
#Calculate v
    v0 = getV(k,alpha,tau,xi = 2*tau, nu = 1, h1 = k,h2 = 1/2*k)
    v0 = max(v0,tau)
  #v_candidate = k_candidate
  
  pi = rep(0,p)
  set.seed(seed)
  loc <- rep(0,M)
  loc[sample(1:M,ceiling(M*(v0-floor(v0))))] = 1
  
#k-FWER part
  for (m in 1:M){
      v = floor(v0)+loc[m]
      W = W0[m,]
      S = getRejectionSet(W,v)$S  
      pi[S] = pi[S]+1
  }
  pi = pi/M
  S = which(pi>=tau)
  return(list(S=S,pi=pi))
}


