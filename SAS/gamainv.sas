proc iml;
/*Funcao Densidade de Probabilidade*/
start dinvgama(x, alpha, beta);
	pdf = ((beta##alpha)/gamma(alpha)) # x##(-alpha - 1) # exp(-beta/x);
	return(pdf);
finish;
/*Funcao Distribuicao Acumulada*/
start pinvgama(q,alpha,beta);
	cdf = 1 - cdf("GAMMA",1/q,alpha,1/beta);
	return(cdf);
finish;
/*Funcao Quantil*/
start qinvgama(p,alpha,beta);
	qf 	= 1/quantile("GAMMA",1-p,alpha,1/beta);
	return(qf);
finish;
/*Funcao Variate*/
start rinvgama(n,alpha,beta);
	aux	= j(n,1);                
	call randgen(aux, "Gamma", alpha, 1/beta);
	rg	= 1/aux; 		
	return(rg);
finish;
quit;
