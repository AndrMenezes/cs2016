OLS <- function(par, x, p){
  xin <- sort(x)
  n <- length(xin)
  i <- 1:n
  S <- sum(((p(xin, par[1], par[2])) - i/(n+1))^2)
  return(S)
}

set.seed(1502)
y <- rinvgama(50, 5, 3)

optim(par = c(5, 3), fn = OLS, x=y, p=pinvgama)