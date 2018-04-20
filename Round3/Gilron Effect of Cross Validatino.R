statistic.levels <- c("Oracle", "Hotelling", "Hotelling.shrink", "Goeman", "sd", "MMD","dCOV", 
                      "lda.CV.1", "lda.noCV.1", "svm.CV.1", "svm.CV.2", "svm.noCV.1", 
                      "svm.noCV.2")


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
suppressPackageStartupMessages(library(sparsediscrim))
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



balance.log <- 'balance_log.txt'
file.remove(balance.log)
file.create(balance.log)
the.message <- paste(file.name)




tr <- function(A) sum(diag(A))

Euclid <- function(x) sqrt(sum(x^2))

gcd <- function(x,y) {
  r <- x%%y;
  return(ifelse(r, gcd(y, r), y))
}


matrix2list <- function(A){
  split(A, rep(1:ncol(A), each = nrow(A)))
}


t_stat <- function(x,y){
  mean.x <- mean(x)
  mean.y <- mean(y)
  var.x <- var(x)
  var.y <- var(y)
  n.x <- length(x)
  n.x.1 <- n.x-1
  n.y <- length(y)
  n.y.1 <- n.y-1
  var.pooled <- (n.x.1 * var.x + n.y.1 * var.y)/(n.x.1+n.y.1)
  nom <- (mean.y-mean.x)^2
  denom <- var.pooled * (1/n.x + 1/n.y)
  t.stat <- nom / denom
  return(t.stat)
}
## Testing:
# t_stat(rnorm(10), rnorm(20,2))

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




# A cross validation wrapper.
# Constructed to allow validation within a permutation test (thus- "old.labels").
t_cv <- function(FUN, noise, new.labels, old.labels, fold.ids, ...){
  n.folds <- max(fold.ids)
  t.cv <- rep(NA, n.folds) # initialize output container
  
  for(v in 1:n.folds){
    # v <- 1
    test.ind <- fold.ids==v
    train.ind <- !test.ind
    
    train.noise <- noise[train.ind, ]
    train.labels <- new.labels[train.ind]
    test.noise <- noise[test.ind, ]
    test.labels <- old.labels[test.ind] 
    
    t.cv[v] <- FUN(train.noise, train.labels, test.noise, test.labels, ...)
  }
  result <- mean(t.cv)
  return(result) # average accuracy over folds
}







t_Oracle <- function(x,y, Sx, Sy){
  x.bar <- colMeans(x)
  y.bar <- colMeans(y)
  delta <- x.bar - y.bar
  
  nx <- nrow(x)
  ny <- nrow(y)
  n <- nx + ny
  nx1 <- nx
  ny1 <- ny
  n1 <- n
  
  S <- (nx1*Sx + ny1*Sy)/n1
  S.inv <- solve(S)
  
  T2 <- delta %*% S.inv %*% delta
  return(T2)
}
### Testing:
# t_Oracle(x = rmvnorm(1e2, rep(0,1e1)),
#          y = rmvnorm(1e2, rep(0,1e1)),
#          Sx=diag(1e1), Sy=diag(1e1))




t_Hotelling <- function(x,y, shrinkage){
  T2 <- hotelling.stat(x,y, shrinkage=shrinkage)$statistic
  return(T2)
}
## Testing:
# t_Hotelling(x = rmvnorm(1e2, rep(0,1e1)), 
#             y = rmvnorm(1e2, rep(0,1e1)), 
#             shrinkage = TRUE)

# Goeman's high-dim test
t_goeman <- function(x,y){
  X <- rbind(x,y)
  dimnames(X) <- list(NULL,LETTERS[1:ncol(X)])
  y <- as.matrix(c(rep(FALSE, nrow(x)),rep(TRUE,nrow(y))))
  result <- globaltest::gt(y,X)
  result@result[,'Statistic'] %>% unname()
}
## Testing
# t_goeman(x1,x2)




# Linear SVM. Arguments self explanatory.
t_svm <- function(train.noise, train.labels, test.noise, test.labels, cost, type){
  svm.1 <- svm(x=train.noise, y=train.labels, type='C-classification', kernel='linear', cost=cost)
  accuracy <- mean(predict(svm.1, newdata=test.noise)==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}

# Compute cross validated test statistics
# noise: the predictors 
# labels: 
# fold.ids: the asignemt of observations to folds
# cost: the svm cost parameter. 
t_svm_cv <- function(noise, new.labels, old.labels, fold.ids, cost, type){
  t_cv(FUN = t_svm, noise, new.labels, old.labels, fold.ids, cost, type)
}




# Linear SVM. Arguments self explanatory.
t_svml2 <- function(train.noise, train.labels, test.noise, test.labels, cost, type){
  # svm.1 <- glmnet(x=train.noise, y=train.labels, family = 'binomial', alpha = 0)
  svm.1 <- LiblineaR(data=train.noise, target =train.labels, type=1, cost=cost)
  predict.labels <- predict(svm.1, newx=test.noise, type='class')$predictions
  accuracy <- mean(predict.labels==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}

# Compute cross validated test statistics
# noise: the predictors 
# labels: 
# fold.ids: the asignemt of observations to folds
# cost: the svm cost parameter. 
t_svml2_cv <- function(noise, new.labels, old.labels, fold.ids, cost, type){
  t_cv(FUN = t_svml2, noise, new.labels, old.labels, fold.ids, cost, type)
}








# Linear Discriminant analysis
t_lda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- lda(x=train.noise, grouping =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}

t_lda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_lda, noise, new.labels, old.labels, fold.ids, type)
}


t_dlda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- dlda(x=train.noise, y =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}

t_dlda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_dlda, noise, new.labels, old.labels, fold.ids, type)
}





t_hdrda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- hdrda(x=train.noise, y =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}



t_hdrda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_hdrda, noise, new.labels, old.labels, fold.ids, type)
}





t_sdlda <- function(train.noise, train.labels, test.noise, test.labels, type){
  lda.1 <- sdlda(x=train.noise, y =train.labels)
  predictions <- predict(lda.1, newdata=test.noise)$class
  accuracy <- mean(predictions==test.labels)
  
  if(type==1){
    statistic <- accuracy
  }
  else if(type==2){
    .p <- mean(test.labels)
    p <- max(.p,1-.p)
    statistic <- abs(accuracy-p)/sqrt(p*(1-p))
  }
  return(statistic)
}



t_sdlda_cv <- function(noise, new.labels, old.labels, fold.ids, type){
  t_cv(FUN = t_sdlda, noise, new.labels, old.labels, fold.ids, type)
}





t_SD <- function(x,y){
  x.bar <- colMeans(x)
  y.bar <- colMeans(y)
  delta <- x.bar - y.bar
  
  Sx <- cov(x)
  Sy <- cov(y)
  
  nx <- nrow(x)
  ny <- nrow(y)
  n <- nx + ny
  nx1 <- nx-1
  ny1 <- ny-1
  n1 <- n-2
  
  S <- (nx1*Sx + ny1*Sy)/n1
  
  D <- diag(S)
  D.inv <- if(length(D)>1) diag(1/D) else 1/D
  
  R <- sqrt(D.inv) %*% S %*% sqrt(D.inv)
  # R <- cor(rbind(x,y))
  p <- ncol(x)
  tr.R2 <- tr(R %*% R)
  d <- 1 + tr.R2/p^(3/2)
  
  nominator <- 1/(1/nx+1/ny) * delta %*% D.inv %*% delta - p
  denominator <- sqrt( 2 * d *(tr.R2 - p^2/n1))
  c(nominator/denominator)
  # list(nom=nominator, den=denominator, t=nominator/denominator)
}
## Testing:
# .p <- 1e0
# .nx <- 1e1
# .ny <- 1e1
# .x <- rmvnorm(.nx, rep(0,.p))
# .y <- rmvnorm(.ny, rep(0,.p))
# t2013(.x,.y)






  





# A Bootstrapping wrapper.
# Sketch:
## For B bootstrap samples:
## Compute statistic
## Average over samples
t_boot <- function(FUN, noise, labels, B, type2, ...){
  
  t.boot <- rep(NA, B) # initialize output container
  n.samples <- nrow(noise)
  
  for(b in 1:B){
    train.ind <- sample(1:n.samples, replace = TRUE)
    
    train.noise <- noise[train.ind, ]
    train.labels <- labels[train.ind]
    test.noise <- noise[-train.ind,]
    test.labels <- labels[-train.ind] 
    
    t.boot[b] <- FUN(train.noise, train.labels, test.noise, test.labels, ...)
  }
  error.boot <- mean(t.boot)
  
  if(type2==2) {
    result <- error.boot
  }
  if(type2==1) {
    error.resubstitute <- FUN(noise, labels, noise, labels, ...)
    result <- 0.368 * error.resubstitute + 0.632 * error.boot
  }
  
  return(result) 
}



t_svm_boot <- function(noise, labels, B, cost, type2, type){
  t_boot(FUN = t_svm, noise = noise, labels = labels, B=B, type2 = type2, cost, type)
}

t_lda_boot <- function(noise, labels, B, type2, type){
  t_boot(FUN = t_lda, noise = noise, labels = labels, B = B, type2, type)
}


t_sdlda_boot <- function(noise, labels, B, type2, type){
  t_boot(FUN = t_sdlda, noise = noise, labels = labels, B = B, type2=type2, type=type)
}




t_svm_highdim_boot <- function(noise, labels, B, cost, type2, type){
  t_boot(FUN = t_svml2, noise = noise, labels = labels, B=B, type2 = type2, cost=cost, type=type)
}


t_kmmd <- function(x,y,...){
  kmmd(x,y,...)@mmdstats[[1]]
}

t_dcov <- function(x,y){
  unname(eqdist.e(rbind(x,y), sizes = c(nrow(x), nrow(y))))
}
## Testing:
# t_dcov(rmvnorm(n = 20,rep(10,10)),rmvnorm(n = 20,rep(0,10)))



t_simes <- function(x,y){
  stopifnot(ncol(y)==ncol(x))
  p <- ncol(x)
  group <- as.factor(c(rep('x', nrow(x)), rep('y',nrow(y))))
  my.t <- function(z) t.test(z~group)$p.value
  p.vals <- apply(rbind(x,y), 2, my.t) # Compute variable-wise pvalues
  p.Simes <- p * min(sort(p.vals)/seq_along(p.vals)) # Compute the Simes statistic
  return(c(pvalue=p.Simes))
}
## Testing
# library(mvtnorm)
# t_simes(rmvnorm(n = 20,rep(0,10)),rmvnorm(n = 30,rep(0,10)))



statistics <- function(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2){
  result <- c(
    Oracle=t_Oracle(x1, x2, Sigma, Sigma),
    Hotelling=t_Hotelling(x1, x2, FALSE),
    Schafer=t_Hotelling(x1, x2, TRUE),
    Goeman=t_goeman(x1, x2),
    Srivastava=t_SD(x1, x2),
    Gretton=t_kmmd(x1, x2),
    dCOV=t_dcov(x1,x2),
    Simes=t_Simes(x1,x2),
    lda.CV.1=t_lda_cv(noise, labels, labels, fold.ids, type=1),
    # lda.CV.2=t_lda_cv(noise, labels, labels, fold.ids, type=2),
    lda.noCV.1=t_lda(noise, labels, noise, labels, type=1),
    # lda.noCV.2=t_lda(noise, labels, noise, labels, type=2),
    svm.CV.1=t_svm_cv(noise, labels, labels, fold.ids, cost=cost.1, type=1),
    svm.CV.2=t_svm_cv(noise, labels, labels, fold.ids, cost=cost.2, type=1),
    # svm.CV.3=t_svm_cv(noise, labels, labels, fold.ids, cost=cost.1, type=2),
    # svm.CV.4=t_svm_cv(noise, labels, labels, fold.ids, cost=cost.2, type=2),
    svm.noCV.1=t_svm(noise, labels, noise, labels, cost=cost.1, type=1),
    svm.noCV.2=t_svm(noise, labels, noise, labels, cost=cost.2, type=1)
    # svm.noCV.3=t_svm(noise, labels, noise, labels, cost=cost.1, type=2),
    # svm.noCV.4=t_svm(noise, labels, noise, labels, cost=cost.2, type=2)
  )
  return(result)
}
## Testing:


statistics.2 <- function(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2){
  c(
    statistics(x1,x2,Sigma,noise,labels,fold.ids,cost.1,cost.2),
    svm.Boot.1=t_svm_boot(noise, labels, B=10, type2=2, cost=cost.1, type=1),
    svm.Boot.2=t_svm_boot(noise, labels, B=10, type2=2, cost=cost.2, type=1),
    svm.Boot.3=t_svm_boot(noise, labels, B=50, type2=2, cost=cost.1, type=1),
    svm.Boot.4=t_svm_boot(noise, labels, B=50, type2=2, cost=cost.2,  type=1),
    lda.Boot.1=t_lda_boot(noise, labels, B=10, type2=2, type=1)
  )
}
