n <- 1e2
p <- 2
noise <- matrix(rnorm(n*p, sd=1), ncol=p, nrow=n)
noise.augment <- augmentDesign(noise)
B <- matrix(1, ncol = p, nrow = p)*2
# plogis(colSums(c(0,0) * (B %*% c(0,0))))
link <- -2+colSums(t(noise) * (B %*% t(noise)))
probs <- plogis(link)
plot(probs,type='h')
labels <- rbinom(n, 1, probs)
plot(noise, col=labels+1, pch=ifelse(labels,15,9), xlab=expression(x[1]), ylab=expression(x[2]))

# pdf(file = '~/workspace/permuting_accuracy/Round4/Output/quadratic-form.pdf')
# plot(noise, col=labels+1, pch=ifelse(labels,15,9), xlab=expression(x[1]), ylab=expression(x[2]))
# dev.off()

