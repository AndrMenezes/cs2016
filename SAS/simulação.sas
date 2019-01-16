%include "C:\<path>\functions.sas";
%include "C:\Users\André Felipe\Dropbox\UEM\3° Série\Estatística Computacional II\Distribuição Gama Inversa\Scripts\metodos_estimacao.sas";

%let B = 100;
%let nmax = 100;
%let alpha = 5;
%let beta  = 2;


proc iml;
load module = _all_;
sup = {1e-10   1e-10,  
        .   	.};
ini = {&alpha &beta};
max = {1, 0};
min = {0, 0};
call randseed(1212);
dados 	= rinvgama(100#&B, &alpha, &beta); 
ysim	= t(shape(dados, &B, &nmax));
est_alpha 	= j(4, &B, 0);
est_beta 	= j(4, &B, 0);
vies_alpha 	= j(4, 10, 0);
vies_beta 	= j(4, 10, 0);	
	do n = 10 to &nmax by 10;
		do j=1 to &B;
			x = ysim[(1:n),j];
			call nlpnra(it, est_mle, "MLE", ini, max, sup);
			call nlpnra(it, est_mps, "MPS", ini, max, sup);
			call nlpnra(it, est_ols, "OLS", ini, min, sup);
			call nlpnra(it, est_wls, "WLS", ini, min, sup);
			mat = est_mle // est_mps // est_ols // est_wls;
			est_alpha[,j] = mat[,1];
			est_beta[,j]  = mat[,2];
		end;
			vies_alpha[,k] 	= mean(est_alpha[(1:4),]) - &alpha;
			vies_beta[,k]	= mean(est_beta[(1:4),]) - &beta;	
	end;
quit;


%let B = 100;
%let nmax = 100;
%let alpha = 5;
%let beta  = 2;
proc iml;
load module = _all_;
sup = {1e-10   1e-10,  
        .   	.};
ini = {&alpha &beta};
max = {1, 0};
min = {0, 0};
est_alpha 	= j(5, &B, 0);
est_beta 	= j(5, &B, 0);
call randseed(1212);
do j=1 to &B;
		x = rinvgama(100, &alpha, &beta); 
		call nlpnra(it, est_mle, "MLE", ini, max, sup);
		call nlpnra(it, est_mps, "MPS", ini, max, sup);
		call nlpnra(it, est_pce, "PCE", ini, min, sup);
		call nlpnra(it, est_ols, "OLS", ini, min, sup);
		call nlpnra(it, est_wls, "WLS", ini, min, sup);
		mat = est_mle // est_mps // est_pce // est_ols // est_wls;
		est_alpha[,j] = mat[,1];
		est_beta[,j]  = mat[,2];
end;
print est_alpha, est_beta;
quit;

%let B = 100;
%let nmax = 100;
%let alpha = 5;
%let beta  = 2;
proc iml;
load module = _all_;
sup = {0   0,  
       .   .};
ini = {&alpha &beta};
max = {1, 0};
min = {0, 0};
est_alpha 	= j(1, &B, 0);
est_beta 	= j(1, &B, 0);
call randseed(1212);
do j=1 to &B;
		x = rinvgama(100, &alpha, &beta); 
		call nlpnra(it, est_mps, "MPS", ini, max, sup);
		est_alpha[,j] = est_mps[,1];
		est_beta[,j]  = est_mps[,2];
end;
print est_alpha, est_beta;
quit;




















data simulados(keep=y);
	call streaminit(222);
	do i=1 to 100;
		x = rand("Gamma",2,3);
		y = 1/x;
		output;
	end;
run;

proc nlp data=simulados tech=nrr;
	max ll;
	bounds alpha > 1e-10, beta > 1e-10;
	parms alpha=2, beta=3;
	ll = alpha*log(beta) + log(gamma(alpha)) - beta*(1/y) - (alpha+1)*log(y);
run;
