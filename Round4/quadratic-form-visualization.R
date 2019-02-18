n <- 100
p <- 2
noise <- matrix(rnorm(n*p, sd=10), ncol=p, nrow=n)
link <- apply(noise, 1, crossprod)
plot(link, type='h')
link2 <- link-median(link)
plot(link2,type='h')
probs <- plogis(link2)
plot(probs,type='h')
labels <- rbinom(n, 1, probs)

plot(noise, col=labels+1, pch=ifelse(labels,15,9), xlab=expression(x[1]), ylab=expression(x[2]))

pdf(file = '~/workspace/permuting_accuracy/Round4/Output/quadratic-form.pdf')
plot(noise, col=labels+1, pch=ifelse(labels,15,9), xlab=expression(x[1]), ylab=expression(x[2]))
dev.off()

