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
suppressPackageStartupMessages(library("here"))

## Get the current path
get_repo_path <- function(){
    repo <- here()
    if (grepl("derandomized_knockoffs_paper", repo)){
        root <- strsplit(repo, "derandomized_knockoffs_paper")[[1]][1]
        return(paste0(root, "derandomized_knockoffs_paper/"))
    }
    repo <- getwd()
    if (grepl("derandomized_knockoffs_paper", repo)){
        root <- strsplit(repo, "derandomized_knockoffs_paper")[[1]][1]
        return(paste0(root, "derandomized_knockoffs_paper/"))
    }
    repo <- Sys.getenv("derandomKnockpath")
    if (grepl("derandomized_knockoffs_paper", repo)){
      root <- strsplit(repo, "derandomized_knockoffs_paper")[[1]][1]
      return(paste0(root, "derandomized_knockoffs_paper/"))
    }
    stop("The derandomized_knockoffs_paper repo is not found!")
}
repo_path <- get_repo_path()
setwd(repo_path)

## Source utility functions
file_vec <- c("crt","pfer_filter","vanilla_pfer_filter")
getfile <- sapply(paste0("./R/",file_vec,".R"),source)
settingName <- "pfer_small"

## Read the pre-computed s for knockoffs
diags <- read.csv("./data/pfer_amp_gaussian.csv")

####################################
## Parameters
####################################
set.seed(24601)
n <- 200
p <- 100
k <- 30
rho <- 0.6
Sigma <- toeplitz(rho^(0:(p-1)))
M <- 31
n_crt <- 3*p
v0 <- 1
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
## derandomized knockoffs
####################################
cat("Running derandomized knockoffs...")
tau <- 0.5
res <- pfer_filter(X,y,v0=v0,M=M,tau=tau,mu = rep(0,p),
                  Sigma=Sigma,diags=diags[,1],
                  seed=seed+24601)
rej <- res$S
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)
pi <- res$pi

save_res <- data.frame(power = power, V=V, method = "dkn")

savedir <- paste0('./results/',settingName,'/pi_amp_',as.character(amp),"_run_",as.character(ParamsRowIndex),'.csv')
write.csv(pi,savedir)

savedir <- paste0('./results/',settingName,'/mdl_amp_',as.character(amp),"_run_",as.character(ParamsRowIndex),'.csv')
write.csv(nonzero,savedir)
cat("done.\n")

####################################
## Vanilla knockoffs
####################################
cat("Running vanilla knockoffs...")
res <- vanilla_pfer_filter(X,y,v0=v0,mu = rep(0,p),
                          Sigma = Sigma,diag=diags[,1],
                          seed = seed+24601)
rej <- res$S
power <- sum(beta0[rej]!=0)/k
V <- sum(beta0[rej]==0)
new_res <- data.frame(power = power, V = V, method = "vkn")
save_res <- rbind(save_res, new_res)
cat("done.\n")

####################################
## Computing CRT p-values and apply Bonforroni
####################################
cat("Running Bonferroni...")
nlambda <- 30
lambda_list <- rep(0,nlambda)
for(i in 1:nlambda){
  Xnew <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  Ynew <- y.sample(Xnew)
  res <- cv.glmnet(Xnew,Ynew,family = "gaussian")
  lambda_list[i] <- res$lambda.min
}
lambda <- median(lambda_list)
p_val = crt(X,y,mu = rep(0,p),Sigma=Sigma,K=n_crt,lambda=lambda,family="gaussian")
rej = which(p_val<=v0/p)
power = sum(beta0[rej]!=0)/k
V = sum(beta0[rej]==0)
new_res <- data.frame(power = power, V=V, method = "bonf")
save_res <- rbind(save_res, new_res)
cat("done.\n")

####################################
## Stability selection
####################################
cat("Running Stability Selection...")
res = stabsel(X,y,fitfun = lars.lasso,cutoff = 0.75,PFER = v0)
rej = res$selected
power = sum(beta0[rej]!=0)/k
V = sum(beta0[rej]==0)
new_res <- data.frame(power = power, V = V, method = "stabs")
save_res <- rbind(save_res, new_res)
cat("done.\n")

## Store the final results
savedir <- sprintf("./results/%s/res_amp_%d_run_%d.csv", settingName, amp, ParamsRowIndex)
cat(sprintf("Results stored to %s...", savedir))
write.csv(save_res,savedir)
cat("done.")
