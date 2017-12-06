library(magrittr)
library(mvtnorm)
library(foreach)
library(doMC)
library(lattice)
library(gap)
library(e1071)
library(Hotelling)
library(MASS)
library(ggplot2)
# library(caret)
library(RPushbullet)
# library(doRNG)
library(reshape2)
library(dplyr)
library(plotly)
library(sparsediscrim)
library(glmnet)
library(LiblineaR)

statistic.filter1 <- c("svm.CV.3","svm.CV.4","svm.noCV.3","svm.noCV.4","lda.noCV.2","lda.CV.2")

tr <- function(A) sum(diag(A))

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

ar1_cov <- function(n, rho){
  result <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      result[i,j] <- round(rho^(abs(i-j)) ,2)
    }
  }
  return(result)
}
## Testing
# ar1_cov(10,0.8)


antiAr1_cov <- function(n, rho){
  result <- matrix(NA, n, n)
  for(i in 1:n){
    for(j in 1:n){
      result[i,j] <- round(rho^(n-abs(i-j)) ,2)
    }
  }
  return(result)
}
## Testing
# lattice::levelplot(antiAr1_cov(10,0.8))




comSym_cov <- function(n, rho){
  result <- diag(n)
  del <- n/2
  for(i in 1:(n/2)) {
    result[i,i+del] <- rho
    result[i+del,i] <- rho
  }
  return(result)
}
## Testing:
# image(comSym_cov(1e1,0.9))


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
sum(replicate(1e3, randomizedTest(0.05, 0.06, 0.03)))



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









t_Hotelling <- function(x,y, shrinkage){
  
  T2 <- hotelling.stat(x,y, shrinkage=shrinkage)$statistic
  return(T2)
  
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
  S.inv <- solve(S)
  
  T2 <- delta %*% S.inv %*% delta
  return(T2)
}
## Testing:
# .p <- 1e0
# .nx <- 1e1
# .ny <- 1e1
# .x <- rmvnorm(.nx, rep(0,.p))
# .y <- rmvnorm(.ny, rep(0,.p))
# t2013(.x,.y)



# Linear SVM. 
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

# Compute cross validated accuracy with l2 regulirized SVM.
t_svm_cv <- function(noise, new.labels, old.labels, fold.ids, cost, type){
  t_cv(FUN = t_svm, noise, new.labels, old.labels, fold.ids, cost, type)
}




# Linear SVM. Arguments self explanatory.
t_svml2 <- function(train.noise, train.labels, test.noise, test.labels, cost, type){
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

# l2 regulirized svm using Liblinear, and K-fold accuracy estimate
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
  lda.1 <- sdlda(x=train.noise, y=train.labels)
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






  





# A Bootstrapped accuracy estimator.
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
## Testing:



t_svm_boot <- function(noise, labels, B, cost, type2, type){
  t_boot(FUN = t_svm, noise = noise, labels = labels, B=B, type2 = type2, cost, type)
}

t_lda_boot <- function(noise, labels, B, type2, type){
  t_boot(FUN = t_lda, noise = noise, labels = labels, B = B, type2, type)
}


t_sdlda_boot <- function(noise, labels, B, type2, type){
  t_boot(FUN = t_sdlda, noise = noise, labels = labels, B = B, type2=type2, type=type)
}
## Testing:
# t_sdlda_boot(noise, labels, B=50, type2=1, type=1)



# Bootstrapped accuracy estimate with l2 regulirized SVM
t_svm_highdim_boot <- function(noise, labels, B, cost, type2, type){
  t_boot(FUN = t_svml2, noise = noise, labels = labels, B=B, type2 = type2, cost=cost, type=type)
}


