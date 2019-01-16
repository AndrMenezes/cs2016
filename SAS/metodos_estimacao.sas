proc iml;
/*Função Densidade de Probabilidade*/
start dinvgama(x, alpha, beta);
	pdf = (beta##alpha/gamma(alpha))#x##(-alpha-1)#exp(-beta/x);
	return(pdf);
finish;
/*Função Distribuição Acumulada*/
start pinvgama(q,alpha,beta);
	cdf = 1 - cdf("GAMMA",1/q,alpha,1/beta);
	return(cdf);
finish;
/*Função Quantil*/
start qinvgama(p,alpha,beta);
	qf = 1/quantile("GAMMA",1-p,alpha,1/beta);
	return(qf);
finish;
/*Função Variate*/
start rinvgama(n,alpha,beta);
	aux	= j(n,1);                
	call randgen(aux, "Gamma", alpha, 1/beta);
	rg	= 1/aux; 		
	return(rg);
finish;

************ Métodos de Estimação *********************;

/*Method of Maximum Likelihood*/
start MLE(par) global(x);
   	alpha	= 	par[1];
   	beta 	= 	par[2];
   	n		=	nrow(x);
	l		=	n*(alpha*log(beta) - log(gamma(alpha))) - beta*sum(1/x) - (alpha + 1)*sum(log(x));
   	return (l);
finish;
/*Method of Maximum Product of Spacings*/
start MPS(par) global(x);
	alpha 	= 	par[1];
	beta	= 	par[2];
	xin 	= 	{0} // x // {1};
	n 		= 	nrow(xin);
	Di 		= 	abs(dif(pinvgama(xin, alpha, beta))[3:n]);
	G 		= 	(1/(n+1)) * sum(log(Di));
	return(G);
finish;
/*Method of Percentiles*/
start PCE(par) global(x);
	alpha 	= 	par[1];
	beta	= 	par[2]; 
	n 		= 	nrow(x);
	i 		= 	1:n;
	pi		= 	i/(n+1);
	P		=	sum((t(x) - qinvgama(pi,alpha,beta))##2);	
	return(P);
finish;
/*Ordinary Least-Squares*/
start OLS(par) global (x);
	alpha 	= 	par[1];
	beta	= 	par[2]; 
	n 		= 	nrow(x);
	i 		= 	1:n;
	Q 		=	sum(((pinvgama(x, alpha, beta)) - t(i/(n+1)) )##2);
	return(Q);
finish;
/*Weighted Least-Squares*/
start WLS(par) global (x);
	alpha 	= 	par[1];
	beta	= 	par[2]; 
	n 		= 	nrow(x);
	i 		= 	1:n;
	W 		=	sum(i#(n-i+1)/((n+1)##2#(n+2))#(( t(pinvgama(x, alpha, beta))) - i/(n+1) )##2);
	return(W);
finish;
/*Gravando os módulos*/
store _all_ module = _all_;
quit;

*************** Exemplos ***********************;
%let alpha = 5;
%let beta = 2;
proc iml;
load module = _all_;

/*Simulando dados*/
call randseed(1515);
x = rinvgama(100,&alpha,&beta); 
call sort(x);

*Passo1: Especificar o suporte dos parâmetros;
sup = { 1e-10   1e-10,  /* limite inferior: 0 < alpha; 0 < beta */
        .   		.}; /* limite superior:  alpha < infty; beta < infty */

*Passo 2: Chamar a rotina para nlpnra;
ini = {&alpha &beta};	/* chutes iniciais */
max = {1};	    /* encontra o máximo da função objetiva*/
min = {0};		/* encontra o máximo da função objetiva*/
call nlpnra(it, est_mle, "MLE", ini, max, sup);
call nlpnra(it, est_mps, "MPS", ini, max, sup);
call nlpnra(it, est_pce, "PCE", ini, min, sup);
call nlpnra(it, est_ols, "OLS", ini, min, sup);
call nlpnra(it, est_wls, "WLS", ini, min, sup);
print est_mle, est_mps, est_pce, est_ols, est_wls;  /*Resultados*/
quit;
