PCE <- function(par, x, q){
  xin <- sort(x)
  n <- length(xin)
  i <- 1:n
  pi <- i / (n+1)  
  P <- sum((xin - q(pi, par[1], par[2]))^2)
  return(P)
}

set.seed(1502)
y <- rinvgama(100, 5, 3)

optim(par = c(5, 3), fn = PCE, x=y, q=qinvgama)