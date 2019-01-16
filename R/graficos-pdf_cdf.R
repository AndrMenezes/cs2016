require(pscl)

y     <- seq(0.01, 5, by=0.01)
alpha <- c(3, 4, 5, 6)
beta  <- c(2, 3, 4, 5)

densigamma
pigamma
qigamma
rigamma

par(mfrow=c(2,2), mar = c(3, 3, 1, 1))
# Alpha variando ----------------------------------------------------------
# PDF ---------------------------------------------------------------------
for (j in 1:length(alpha)) {
  fy  <- densigamma(x=y, alpha=alpha[j], beta=3)
  if(j==1){
    plot(y, fy, type="l", col=j, lwd=2,xlab="",ylab='',ylim=c(0,2.2))
    mtext("y", side = 1, line = 2)
    mtext("f(y)", side = 2, line = 1.7)
  }else{lines(y, fy, col=j, lwd=2)}
}
legend("topright", paste("GI(", alpha, ";3)",sep=""),lwd=2,inset=0.04,
       col=1:4)
# CDF ---------------------------------------------------------------------
for (j in 1:length(alpha)) {
  fy  <- pigamma(q=y, alpha=alpha[j], beta=3)
  if(j==1){
    plot(y, fy, type="l", col=j, lwd=2,xlab="",ylab='')
    mtext("y", side = 1, line = 2)
    mtext("F(y)", side = 2, line = 1.7)
  }else{lines(y, fy, col=j, lwd=2)}
}
legend("bottomright", paste("GI(", alpha, ";3)",sep=""),lwd=2,inset=0.04,
       col=1:4)


# Beta variando -----------------------------------------------------------
# PDF ---------------------------------------------------------------------
for (j in 1:length(alpha)) {
  fy  <- densigamma(x=y, alpha=3, beta=beta[j])
  if(j==1){
    plot(y,fy,type="l",col=j,lwd=2,xlab="",ylab="",ylim=c(0,1.2))
    mtext("y", side = 1, line = 2)
    mtext("f(y)", side = 2, line = 1.7)
  }else{lines(y, fy, col=j, lwd=2)}
}
legend("topright", paste("GI(3;",beta,")",sep=""),lwd=2,inset=0.04,
       col=1:4)

# CDF ---------------------------------------------------------------------
for (j in 1:length(alpha)) {
  fy  <- pigamma(q=y, alpha=3, beta=beta[j])
  if(j==1){
    plot(y,fy,type="l",col=j,lwd=2,xlab="",ylab="")
    mtext("y", side = 1, line = 2)
    mtext("F(y)", side = 2, line = 1.7)
  }else{lines(y, fy, col=j, lwd=2)}
}
legend("bottomright", paste("GI(3;",beta,")",sep=""),lwd=2,inset=0.04,
       col=1:4)



# MLE ---------------------------------------------------------------------
x = rgamma(100, 2, 3)
y = 1/x
lh <- function(par, x){
  alpha <- par[1]
  beta <- par[2]
  # sum(log(densigamma(x,alpha,beta)))
  n <- length(x)
  n*alpha*log(beta) - n*log(gamma(alpha)) + beta*sum(1/x) - (alpha+1)*sum(log(x));
}
optim(par=c(10,3), lh, x=y, control = list(fnscale=-1))



x = rigamma(100, 2, 3)
U <- function(par, x){
  alpha <- par[1]
  beta <- par[2]
  n <- length(x)
  I <- matrix(nrow=2)  
  I[1,] <- n*log(beta) - n*digamma(alpha) - sum(log(x))
  I[2,] <- n*alpha/(beta) + sum(1/x)
  return(I)  
}

library(rootSolve)
multiroot(f = U, start = c(2,3), x=x)







