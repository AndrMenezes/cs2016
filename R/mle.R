library(fitdistrplus)

logLH <- function(par, x){
  alpha <- par[1]
  beta <- par[2]
  n <- length(x)
  return(n*(alpha*log(beta) - log(gamma(alpha))) - beta*sum(1/x) - (alpha + 1)*sum(log(x)))
}

set.seed(1502)
y <- rinvgama(100, 5, 3)

optim(par = c(5, 3), fn = logLH, x=y, control = list(fnscale=-1))
fitdist(data = y, distr = "invgama", start = c(5, 3), method = "mle")