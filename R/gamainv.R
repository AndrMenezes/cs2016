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