WLS <- function(par, x, p){
  xin <- sort(x)
  n <- length(xin)
  i <- 1:n
  W <- sum(i*(n-i+1)/((n+1)^2*(n+2))*((p(xin,par[1], par[2]))-i/(n+1))^2)
  return(W)
}

set.seed(1502)
y <- rinvgama(50, 5, 3)

optim(par = c(5, 3), fn = WLS, x=y, p=pinvgama)