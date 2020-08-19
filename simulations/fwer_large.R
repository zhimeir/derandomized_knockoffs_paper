#!/usr/bin/env Rscript

## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
amp <- as.integer(args[1])
ParamsRowIndex <- as.integer(args[2])
if(is.na(amp)==1){
  amp <- 4
}
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex <- 1
}
####################################
## Libraries and sources
####################################
suppressPackageStartupMessages(library("knockoff"))
suppressPackageStartupMessages(library("SNPknock"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("corpcor"))
suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("stabs"))

file_vec <- c("crt","fwer_filter","getV","vanilla_fwer_filter")
getfile <- sapply(paste0("../R/",file_vec,".R"),source)
settingName <- "fwer_large"

####################################
## Parameters
####################################
set.seed(24601)
n = 2000
p = 1000
k = 60
rho = 0.2
Sigma = toeplitz(rho^(0:(p-1)))
M=30 #Sample number
alpha=0.1
k_target= 3
nonzero = sample(1:p,k)
beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p,0,1))/sqrt(n)
y.sample = function(X) rbinom(1,1,exp(X%*%beta0)/(1+exp(X%*%beta0)))
seed = as.integer(amp*ParamsRowIndex)

####################################
## Generating data
####################################
set.seed(seed)
X = matrix(rnorm(n*p),n) %*% chol(Sigma)
y = apply(X,1,y.sample)

####################################
## derandomized knockoffs 
####################################
cat("Running derandomized knockoffs...")
tau = 0.5
res = fwer_filter(X,y,k=k_target, M=M,tau = tau,alpha = alpha,
                  knockoff_stat = stat.lasso_coefdiff_bin,
                  mu = rep(0,p),
                  Sigma = Sigma,seed = seed+24601)
rej = res$S
power = sum(beta0[rej]!=0)/k
V = sum(beta0[rej]==0)
pi = res$pi

save_res <- data.frame(power = power, V=V, method = "dkn")

savedir = paste0('../results/',settingName,'/pi_amp_',as.character(amp),"_run_",as.character(ParamsRowIndex),'.csv')
write.csv(pi,savedir)

savedir = paste0('../results/',settingName,'/mdl_amp_',as.character(amp),"_run_",as.character(ParamsRowIndex),'.csv')
write.csv(nonzero,savedir)

cat("done.\n")

####################################
## Vanilla knockoffs
####################################
cat("Running vanilla knockoffs...")
res = vanilla_fwer_filter(X,y,k=k_target,alpha=alpha,
                          knockoff_stat = stat.lasso_coefdiff_bin,
                          mu = rep(0,p),
                          Sigma = Sigma,seed = seed+24601)
rej = res$S
power = sum(beta0[rej]!=0)/k
V = sum(beta0[rej]==0)
save_res <- rbind(save_res,c(power,V,"vkn"))
cat("done.\n")
####################################
## Computing CRT p-values and apply Bonforroni
####################################
cat("Running Bonferroni...")
bin_pval = function(x,y){
  res = t.test(x[y==0],x[y==1],alternative = "two.sided")
  return(res$p.value)
}
p_val = apply(X=X,MARGIN = 2,FUN = bin_pval,y=y)
rej = which(p_val<=k_target*alpha/p)
power = sum(beta0[rej]!=0)/k
V = sum(beta0[rej]==0)
save_res <- rbind(save_res,c(power,V,"bonf"))
savedir = paste0('../results/',settingName,'/res_amp_',as.character(amp),"_run_",as.character(ParamsRowIndex),'.csv')
write.csv(save_res,savedir)
cat("Done.")
