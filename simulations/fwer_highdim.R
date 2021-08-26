#!/usr/bin/env Rscript

## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
amp <- as.integer(args[1])
ParamsRowIndex <- as.integer(args[2])
if(is.na(ParamsRowIndex)){
  ParamsRowIndex <- 1 
}
if(is.na(amp)){
  amp <- 35 
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
suppressPackageStartupMessages(library("R.matlab"))
suppressPackageStartupMessages(library("SNPknock"))
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("RCurl"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("httr"))
suppressPackageStartupMessages(library("stabs"))

source_gitfile <- function(filename){
  url = paste0("https://raw.githubusercontent.com/zhimeir/skn/master/R/",filename,".R")
  script = GET(url = url,authenticate("7cd7f35469c8b2759863e4be7cc95abf7fb98abf",""))
  script = content(script,"text")
  eval(parse(text = script),envir= .GlobalEnv)
}

file_vec <- c("crt","fwer_filter","getV","getOutput","vanilla_fwer_filter")
getfile <- sapply(file_vec,source_gitfile)

####################################
## Parameters
####################################
set.seed(24601)
n <- 200
p <- 2000
k <- 30
rho <- 0.2
Sigma <- toeplitz(rho^(0:(p-1)))
M <- 30 
alpha <- 0.1
k_target <- 2
nonzero <- sample(1:p,k)
beta0 <- amp * (1:p %in% nonzero)*sign(rnorm(p,0,1))/sqrt(n)
y.sample <- function(X) X%*%beta0 + rnorm(n,1)
#y.sample <- function(X) rbinom(1,1,exp(X%*%beta0)/(1+exp(X%*%beta0)))
seed <- as.integer(amp*ParamsRowIndex)

####################################
## Generating data
####################################
set.seed(seed)
X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
#y <- apply(X,1,y.sample)
y <- y.sample(X)

####################################
## stable knockoff 
####################################
tau <- .81
res <- fwer_filter(X, y, k = k_target, alpha = alpha,
                   M = M, tau = tau, v0 = 1,
                   knockoff_stat = stat.lasso_coefdiff,
                   mu = rep(0,p),Sigma = Sigma,seed = seed+24601)
rej <- res$S
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)
pi <- res$pi

savedir <- paste0('./results/',settingName,'/skn/res_',as.character(amp),"_",as.character(ParamsRowIndex),'.mat')
writeMat(savedir,power=power,V=V,pi = pi,nonzero = nonzero)


####################################
## Vanilla knockoff
####################################
res <- vanilla_fwer_filter(X, y, k = k_target, alpha = alpha,
                          knockoff_stat = stat.lasso_coefdiff,
                          mu = rep(0,p),
                          Sigma = Sigma,seed = seed+24601)
rej <- res$S
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)

savedir <- paste0('./results/',settingName,'/vkn/res_',as.character(amp),"_",as.character(ParamsRowIndex),'.mat')
writeMat(savedir,power=power,V=V,nonzero = nonzero)

####################################
## Stability Selection
####################################
res <- stabsel(X,y,fitfun = lars.lasso,
               cutoff = 0.75,PFER = alpha*k_target)

rej <- res$selected
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)

savedir = paste0('./results/',settingName,'/stabs/res_',as.character(amp),"_",as.character(ParamsRowIndex),'.mat')
writeMat(savedir,power=power,V=V,nonzero = nonzero)

####################################
## Computing CRT p-values and apply Bonforroni
####################################
## bin_pval = function(x,y){
##   res = t.test(x[y==0],x[y==1],alternative = "two.sided")
##   return(res$p.value)
## }
## p_val = apply(X=X,MARGIN = 2,FUN = bin_pval,y=y)
## rej = which(p_val<=k_target*alpha/p)
## power = sum(beta0[rej]!=0)/k
## V = sum(beta0[rej]==0)
## 
## savedir = paste0('./results/',settingName,'/crt/res_',as.character(amp),"_",as.character(ParamsRowIndex),'.mat')
## writeMat(savedir,power=power,V=V)
