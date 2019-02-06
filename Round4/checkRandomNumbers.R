RhpcBLASctl::blas_set_num_threads(1)
registerDoMC(cores=100)
library(doRNG)
registerDoRNG()

data <- foreach(j=1:1e3, .combine=rbind) %dopar%{
  runif(1e0)
}

dim(data)
plot(data, type='l', cex=0.1)
matplot(data, type='l')
    