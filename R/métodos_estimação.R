# Funcao densidade de probabilidade 
dinvgama <- function(x, alpha, beta){
  (beta^alpha)/gamma(alpha)*x^(-alpha - 1)*exp(-beta/x)
}

# Funcao distribuicao acumulada
pinvgama <- function(q, alpha, beta){
  1 - pgamma(1/q, alpha, beta)
}

# Funcao quantil
qinvgama <- function (p, alpha, beta) {
  1/qgamma(1 - p, alpha, beta)
}

# Funcao variate
rinvgama <- function(n, alpha, beta){
  1/rgamma(n, alpha, beta)
}

###################### Métodos de Estimação #############################

# Method of Maximum Likelihood --------------------------------------------
logLH <- function(par, x){
  alpha <- par[1]
  beta <- par[2]
  n <- length(x)
  return(n*(alpha*log(beta) - log(gamma(alpha))) - beta*sum(1/x) - (alpha + 1)*sum(log(x)))
}

# Method of Maximum Product of Spacings -----------------------------------
MPS <- function(par, x, p){
  xin <- sort(x)
  n <- length(xin)
  Di <- diff(c(0, p(xin, par[1], par[2]),1))
  G <- (1/(n+1)) * sum(log(Di))
  return(G)
}

# Method of Percentiles ---------------------------------------------------
PCE <- function(par, x, q){
  xin <- sort(x)
  n <- length(xin)
  i <- 1:n
  pi <- i / (n+1)  
  P <- sum((xin - q(pi, par[1], par[2]))^2)
  return(P)
}

# Ordinary Least-Squares --------------------------------------------------
OLS <- function(par, x, p){
  xin <- sort(x)
  n <- length(xin)
  i <- 1:n
  S <- sum(((p(xin, par[1], par[2])) - i/(n+1))^2)
  return(S)
}

# Weighted Least-Squares --------------------------------------------------
WLS <- function(par, x, p){
  xin <- sort(x)
  n <- length(xin)
  i <- 1:n
  W <- sum(i*(n-i+1)/((n+1)^2*(n+2))*((p(xin,par[1], par[2]))-i/(n+1))^2)
  return(W)
}


# Exemplos ----------------------------------------------------------------
set.seed(1502)
y <- rinvgama(50, 5, 3)
optim(par = c(5, 3), fn = logLH, x=y, control = list(fnscale=-1))
optim(par = c(5, 3), fn = MPS, x=y, p=pinvgama, control = list(fnscale=-1))
optim(par = c(5, 3), fn = PCE, x=y, q=qinvgama)
optim(par = c(5, 3), fn = OLS, x=y, p=pinvgama)
optim(par = c(5, 3), fn = WLS, x=y, p=pinvgama)

rm(y)
