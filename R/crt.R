crt <- function(X,Y,mu,Sigma,K,lambda,family = "gaussian"){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #run cross-validated lasso on the original dataset
  res <- glmnet(X,Y,family = family,lambda=lambda)
  beta <- res$beta
  crt_pval <- sapply(1:p,crt_marginal,X0=X,Y=Y,
                     mu = mu,Sigma=Sigma,K=K,
                     beta = beta,
                     lambda = lambda,
                     family = family) 
  return(crt_pval)
}

## Computing the CRT p-value for feature j
crt_marginal <- function(X0,Y,mu,Sigma,lambda,beta,j,K,family){
  n <- dim(X)[1]
  inv_sigma <- solve(Sigma[-j,-j])
  cond_mu <- Sigma[j,-j]%*%inv_sigma%*%t(X0[,-j] - t(replicate(n,mu[-j])))
  cond_sigma <- Sigma[j,j] - Sigma[j,-j]%*%inv_sigma%*%Sigma[-j,j]
  crt_stats <- sapply(X=1:K, FUN=crt_sampling,
                      X0=X0,Y=Y,j=j,
                      cond_mu = cond_mu,
                      cond_sigma = cond_sigma,
                      lambda = lambda,family = family)
  pval <- (1+sum(crt_stats>=abs(beta[j])))/(K+1)
  return(pval)
}

## Sample Xj conditional on other covariates and compute the test stat
crt_sampling <- function(X0,Y,j,cond_mu,cond_sigma,lambda = lambda,seed,family){
  set.seed(seed)
  new_xj <- sapply(cond_mu,rnorm,n=1,sd = cond_sigma)
  new_x <- X0
  new_x[,j] <- new_xj
  res <- glmnet(new_x,Y,lambda = lambda,family = family)
  beta <- res$beta
  return(abs(beta[j]))  
}
