%{
	
2020.07.30
AFZ

Testing behavior of CIs for boundary cases.

%}


%--------------------------------------------------------------------------%
%%  Preliminaries 
%--------------------------------------------------------------------------%

%  Clear graphics 

%  Set seed for replicability 
rng(12345)

%  File paths 
clear ri_ci 
addpath('../m/'); 

%  Parallel processing 
R = 100 ; % number of randomizations/permutations to consider 
RunParallel = true ; 
Noisily 	= true ; 

if RunParallel 
	pool = gcp('nocreate') ; 
	if isempty(pool) 
		mycluster = parcluster('local') ; 
		if isunix 
            mycluster.NumWorkers = str2double(getenv('NSLOTS'));
		else 
            mycluster.NumWorkers = 16;
		end
		parpool(mycluter, mycluster.NumWorkers) ; 
	end
end

%--------------------------------------------------------------------------%
%%  Construct data  
%--------------------------------------------------------------------------%

%  Parameters
N = 2000;
sigma2 = 1 ; % variance of the error term
ptreat = 0.5 % fraction treated 

%  Set treatment effect to deliver analytic p-value of 0.05.
se_analytic = sqrt(sigma2 / (N*ptreat*(1-ptreat))) ; 
tau = 1.96*se_analytic ; 
TAU = transpose([-.25:.05:.25]*se_analytic + tau) ; % set of possible values to consider 

RESULTS = array2table([TAU,NaN(length(TAU),4)],'VariableNames',{'tau','tau_hat','p_analytic', 'p_ri', 'ci_lower' }); % storing analytic results

%  Potential randomizations 
T0 = double(tiedrank(rand(N,R))/N > ptreat ) ; % tiedrank operates within columns as required

%  DGP for actual sample 
e = randn(N,1) * sqrt(sigma2) ; 
t =  double(tiedrank(rand(N,1)) / N > ptreat ); % assign binary treatment 
y = tau * t + e ; % outcomes determined

%  Data as table 
data = array2table([y,t,e],'VariableNames',{'y','t','e'}); 
data.constant = ones(N,1); 

%  Run regression and check analytical standard errors for closeness to p-value 
mdl = fitlm(data,'y ~ t')

%  First pass with ri_ci() 
% [b,N,pval] = ri_ci( ...
% 	data, ...
% 	'y', ...
% 	't', ...
% 	'T0', T0, ... 
% 	'P', R, ...
% 	... 'Controls', {'constant'}, ... 
% 	'TestZero',true, ...
% 	'FindCI', false ...
% 	) 

%  Simulate across alternative DGPs with varying treatment effect sizes
for tt = 1 :length(TAU) 
	y = TAU(tt)*t + e ;  
	data.y = y ; 

	%  Regression estimates and analytic standard errors 
	mdl = fitlm(data,'y ~ t') ; 
	RESULTS(tt,'tau_hat') = mdl.Coefficients('t','Estimate') ;
	RESULTS(tt,'p_analytic') = mdl.Coefficients('t','pValue') ; 

	%  RI p-values 
	[b,~,pval, ci] = ri_ci( ...
		data ...
		, 'y' ...
		, 't' ...
		, 'T0', T0 ...
		, 'P', R ...
		, 'TestZero', true ...
		, 'FindCI', true ...
		, 'RunParallel', RunParallel ...
		, 'Noisily', Noisily ... 
		)
	RESULTS(tt,'p_ri') = array2table(pval) ; 
	RESULTS(tt,'ci_lower') = array2table(ci(1)); 
end 


figure(1)
clf 
hold on 
plot(RESULTS.tau,RESULTS.p_analytic, 'DisplayName', 'Analytic')
plot(RESULTS.tau,RESULTS.p_ri, 'DisplayName', 'RI') 
yline(0.05) 
xlabel('Treatment effect (DGP)') 
ylabel('p-value') 
legend('Analytic', 'RI')  
hold off 
