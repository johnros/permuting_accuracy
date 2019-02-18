suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doMC))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(gap))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(Hotelling))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(RPushbullet))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plotly))
# suppressPackageStartupMessages(library(sparsediscrim))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(LiblineaR))
suppressPackageStartupMessages(library(data.table))
# source("https://bioconductor.org/biocLite.R")
# biocLite("globaltest")
suppressPackageStartupMessages(library('globaltest'))
suppressPackageStartupMessages(library(RobPer))
# suppressPackageStartupMessages(library(fungible))
suppressPackageStartupMessages(library(kernlab))
suppressPackageStartupMessages(library(energy))
suppressPackageStartupMessages(library(HDtest))
suppressPackageStartupMessages(library(dqrng))



balance.log <- 'balance_log.txt'
file.remove(balance.log)
file.create(balance.log)
the.message <- paste(file.name)





tr <- function(A) sum(diag(A))
Frob <- function(A) sum(A^2)

Euclid <- function(x) sqrt(sum(x^2))

gcd <- function(x,y) {
  r <- x%%y;
  return(ifelse(r, gcd(y, r), y))
}


matrix2list <- function(A){
  split(A, rep(1:ncol(A), each = nrow(A)))
}


ar1_cov <- function(n, rho, sigma=1){
  x <- diag(n) 
  x <- sigma * rho^abs(row(x)-col(x))
  x
  return(x)
}
## Testing
# lattice::levelplot(ar1_cov(10,0.8))

# The correlation implies by a brownian motion
browninan_cov <- function(n){
  Sigma <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      Sigma[i,j] <- min(i,j)
    }
  }
  D <- 1/sqrt(diag(Sigma))
  Corr <- diag(D) %*% Sigma %*% diag(D)
  return(Corr)
}
## Testing
# lattice::levelplot(browninan_cov(10,0.8))



browninan_cov2 <- function(n){
  Sigma <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      Sigma[i,j] <- min(i/n,j/n)
    }
  }
  return(Sigma)
}
## Testing
# lattice::levelplot(browninan_cov2(10))


seq_cov <- function(p){
  Sigma0 <- 1:p
  Sigma.0.trace <- sum(Sigma0)
  # Sigma <- diag(Sigma0/Sigma.0.trace)
  Sigma <- diag(Sigma0)
  return(Sigma)
}






balanced_folding <- function(labels, n.folds, balance){
  stopifnot(is.logical(balance))
  n <- length(labels)
  
  if(balance){
    folds.list <- createFolds(y=labels, n.folds)
    fold.inds <- rep(NA, n)
    for(v in seq_along(folds.list)) fold.inds[folds.list[[v]]] <- v 
  }
  else {
    fold.inds <- sample(rep(1:n.folds, length=n))
  }
  return(fold.inds)  
}


my.summary <- function(pvals){
  cat(
    sprintf('cost=%s, permutations=%d, replications=%d, n=%d, p=%d, train.prop=%s, effect=%s, n.folds=%d', 
            cost, n.permutations, n.replications, n, p, train.prop, effect, n.folds), 
    '\n',
    apply(pvals, 2, function(x) mean(x<0.05)))
}

my.ecdf <- function(x,t){
  # ecdf(-x)(-t) #+ 1/(length(x)+1)
  mean(x>=t)
}


randomizedTest <- function(alpha, pval, pval.strickt){
  stopifnot(pval>=pval.strickt)
  reject <- NULL
  if(pval.strickt> alpha) reject <- FALSE
  else if (pval<alpha) reject <- TRUE
  else {
    p.diff <- (alpha-pval.strickt)/(pval- pval.strickt)
    reject <- rbinom(1, 1, p.diff) %>% as.logical
  }
  return(reject)
}
## Testing:
# sum(replicate(1e3, randomizedTest(0.05, 0.06, 0.03)))


# statistic.levels <- c("Oracle", "Hotelling", "Hotelling.shrink", "Goeman", "sd", "MMD","dCOV", 
#                       "lda.CV.1", "lda.noCV.1", "svm.CV.1", "svm.CV.2", "svm.noCV.1", 
#                       "svm.noCV.2")


source('statistics.R')

makeSymmetric <- function(n){
  x <- matrix(rnorm(n*n), n) 
  ind <- lower.tri(x) 
  x[ind] <- t(x)[ind] 
  x
}
## Testing:
# makeSymmetric(3)


makePrecision <- function(X,p,SigmaRaw){
  # X <- noise.augment
  # SigmaRaw <- Sigma
  
  vars <- colnames(X)
  new.p <- p+1+choose(p,2)
  Sigma <- Matrix(0, nrow = new.p, ncol = new.p, sparse = TRUE)
  Sigma[1:p,1:p] <- SigmaRaw
  split ###
  
}