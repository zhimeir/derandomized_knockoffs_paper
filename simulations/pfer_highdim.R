#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
amp <- as.integer(args[1])
ParamsRowIndex <- as.integer(args[2])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex <- 1 
}
if(is.na(amp)==1){
  amp <- 50
}
####################################
## Libraries and sources
####################################
library("knockoff")
library("SNPknock")
library("dplyr")
library("corpcor")
library("glmnet")
library("MASS")
library("R.matlab")
library("SNPknock")
library("devtools")
library("RCurl")
library("ggplot2")
library("httr")
library("stabs")

source_gitfile <- function(filename){
  url = paste0("https://raw.githubusercontent.com/zhimeir/skn/master/R/",filename,".R")
  script = GET(url = url,authenticate("513fae7de06c7682ecf82ac03513241987ef7926",""))
  script = content(script,"text")
  eval(parse(text = script),envir= .GlobalEnv)
}

file_vec <- c("crt","pfer_filter","getV","getOutput","vanilla_pfer_filter")
getfile <- sapply(file_vec,source_gitfile)

####################################
## Parameters
####################################
set.seed(24601)
n <- 200
p <- 2000
k <- 30
rho <- 0.5
Sigma <- toeplitz(rho^(0:(p-1)))
M <- 51
v0 <- 2
nonzero <- sample(1:p,k)
beta0 <- amp * (1:p %in% nonzero)*sign(rnorm(p,0,1))/sqrt(n)
y.sample <- function(X) X%*%beta0+ rnorm(n,1)
seed <- as.integer(amp*ParamsRowIndex)

####################################
## Generating data
####################################
set.seed(seed)
X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
y <- y.sample(X)

####################################
## stable knockoff 
####################################
tau <- 0.5
res <- pfer_filter(X,y,v0=v0, M=M,tau = tau, mu = rep(0,p),
                  Sigma = Sigma,seed = seed+24601)
rej <- res$S
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)
pi <- res$pi

savedir <- paste0('./results/',settingName,'/skn/res_',as.character(amp),"_",as.character(ParamsRowIndex),'.mat')
writeMat(savedir,power=power,V=V,pi = pi,nonzero = nonzero)

####################################
## Vanilla knockoff
####################################
res <- vanilla_pfer_filter(X,y,v0=v0,mu = rep(0,p),
                          Sigma = Sigma,seed = seed+24601)
rej <- res$S
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)

savedir <- paste0('./results/',settingName,'/vkn/res_',as.character(amp),"_",as.character(ParamsRowIndex),'.mat')
writeMat(savedir,power=power,V=V,nonzero = nonzero)


####################################
## Stability selection
####################################
res <- stabsel(X,y,fitfun = lars.lasso,cutoff = 0.75,PFER = v0)
rej <- res$selected
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)

savedir <- paste0('./results/',settingName,'/stabs/res_',as.character(amp),"_",as.character(ParamsRowIndex),'.mat')
writeMat(savedir,power=power,V=V)
