n <- 1e2
p <- 2
noise <- matrix(rnorm(n*p), ncol=p, nrow=n)
# noise.augment <- augmentDesign(noise)
b0 <- -2
beta0 <- rep(0,p)
beta <- beta0 
B0 <- matrix(1, ncol = p, nrow = p)/sqrt(p)
B <- B0*2
# plogis(colSums(c(0,0) * (B %*% c(0,0))))
link <- b0 + noise%*%beta + colSums(t(noise) * (B %*% t(noise)))
probs <- plogis(link)
labels <- c(TRUE,FALSE)[rbinom(n,1, probs)+1]
plot(probs,type='h', col=labels+1)
plot(noise, col=labels+1, pch=ifelse(labels,15,9), xlab=expression(x[1]), ylab=expression(x[2]))

# pdf(file = '~/workspace/permuting_accuracy/Round4/Output/quadratic-form.pdf')
# plot(noise, col=labels+1, pch=ifelse(labels,15,9), xlab=expression(x[1]), ylab=expression(x[2]))
# dev.off()

