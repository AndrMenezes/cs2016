MPS <- function(par, x, p){
  xin <- sort(x)
  n <- length(xin)
  Di <- diff(c(0, p(xin, par[1], par[2]),1))
  H <- (1/(n+1)) * sum(log(Di))
  return(H)
}

set.seed(1502)
y <- rinvgama(100, 5, 3)

optim(par = c(5, 3), fn = MPS, p=pinvgama, x=y, control = list(fnscale=-1))