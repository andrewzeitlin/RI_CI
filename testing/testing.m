%{
	
2020.07.30
AFZ

Testing behavior of CIs for boundary cases.

%}


%--------------------------------------------------------------------------%
%%  Preliminaries 
%--------------------------------------------------------------------------%

rng('default')



%--------------------------------------------------------------------------%
%%  Construct data  
%--------------------------------------------------------------------------%

%  Parameters
N = 1000;
sigma2 = 1 ; % variance of the error term
R = 1000; % number of randomizations/permutations to consider 
ptreat = 0.5 % fraction treated 

se_analytic = sqrt(sigma2 / (N*ptreat*(1-ptreat))) ; 
tau = 1.96*se_analytic ; 

%  Potential randomizations 
T0 = double(tiedrank(rand(N,R))/N > ptreat ) ; % tiedrank operates within columns as required

%  DGP for actual sample 
e = randn(N,1) * sqrt(sigma2) ; 
t =  double(tiedrank(rand(N,1)) / N > ptreat ); % assign binary treatment 
y = tau*t + e ; % outcomes determined

data = array2table([y,t],'VariableNames',{'y','t'}); 

%  Run regression and check analytical standard errors for closeness to p-value 
mdl = fitlm(data,'y ~ t')


