source("gamainv.R")
source("métodos_estimação.R")

B = 10000
alpha = c(1, 2, 3, 4)
beta = 1
ni = c(10, 20, 40, 60, 80, 100)

est_alpha <- matrix(nrow=5, ncol=B)
est_beta <- matrix(nrow=5, ncol=B)

mat <- matrix(nrow=5, ncol=length(ni), 
               dimnames = list(c("MLE","MPS","PCE","OLS","WLS"), ni))
vies_alpha <- list(mat, mat, mat, mat)
names(vies_alpha) <- paste(alpha, "e" ,beta)
eqm_beta <- eqm_alpha <- vies_beta <- vies_alpha


inicio <- proc.time()
for(l in 1:length(alpha)) {
  k <- 1
  ini <- c(alpha[l], beta)
  set.seed(1502)
  dados <- rinvgama(B*max(ni), alpha[l], beta)
  X  <- matrix(dados, nrow=max(ni), ncol=B)
  for(i in ni){
    for(j in 1:B){
      xx  <- X[(1:i), j]
      aux <- matrix(c(optim(par=ini,fn=logLH,x=xx,control=list(fnscale = -1))$par,
                      optim(par=ini,fn=MPS,x=xx,p=pinvgama,control=list(fnscale = -1))$par,
                      optim(par=ini,fn=PCE,x=xx,q=qinvgama)$par,
                      optim(par=ini,fn=OLS,x=xx,p=pinvgama)$par,
                      optim(par=ini,fn=WLS,x=xx,p=pinvgama)$par),nrow=5, ncol=2,byrow=T)
      est_alpha[,j] <- aux[,1]
      est_beta[,j] <- aux[,2]
    }
    vies_alpha[[l]][,k] <- rowMeans(est_alpha) - alpha[l]
    vies_beta[[l]][,k] <- rowMeans(est_beta) - beta
    eqm_alpha[[l]][,k] <- rowMeans((est_alpha - alpha[l])^2)
    eqm_beta[[l]][,k] <- rowMeans((est_beta - beta)^2)
    k <- k+1;  
  }
}
(fim <- proc.time() - inicio)

saveRDS(vies_alpha, "vies_alpha.rds")
saveRDS(vies_beta, "vies_beta.rds")
saveRDS(eqm_alpha, "eqm_alpha.rds")
saveRDS(eqm_beta, "eqm_beta.rds")
