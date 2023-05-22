source('ConstrLasso.R')
n <- 50
p <- 100
s <- 6
sig <- 0.5

x <- matrix(rnorm(n*p), n, p)
x <- scale(x,center=rexp(p)/10, scale=c(1,5,0.1,10,2,6,runif(p-6)))
bet0 <- c(seq(from=-2,to=2,length.out=s),rep(0,p-s))
y <- x%*%bet0 + rnorm(n)*sig
C <- matrix(1,p,1)
C <- cbind(C, c(-1,1,1,-1, rep(0,p-4)))
t(C)%*%bet0

res1 <- ConstrLasso(y,x,C)
t(C)%*%res1$bet

res2 <- ConstrLassoCrossVal(y=y,x=x,C=C)
t(C)%*%res1$bet
plot(bet0,res2$bet[,which.min(res2$cvm)])
plot(res2$cvm)
