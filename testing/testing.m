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
rng('default')

%  File paths 
clear ri_ci 
addpath('../m/'); 

%  Parallel processing 
R = 1000 ; % number of randomizations/permutations to consider 
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
		parpool(mycluster, mycluster.NumWorkers) ; 
	end
end

%--------------------------------------------------------------------------%
%%  Construct data  
%--------------------------------------------------------------------------%

%  Parameters
N = 2000;
ptreat = 0.5 % fraction treated 
t =  double(tiedrank(rand(N,1)) / N > ptreat ); % assign binary treatment 
SST_t = sum((t - mean(t)).^2); 
sigma2 = SST_t / (1.96^2) ; % Set error variance so that coefficient of 1 has expected analytic p-value of 0.05 

%  DGP for actual sample 
e = randn(N,1) * sqrt(sigma2) ; 
b_e = [ones(size(t)),t]\e ; % regress error terms on treatment 
ehat = [ones(size(t)), t] * b_e ;
e = e - ehat ; %  Residualize so that error terms no longer correlated with treatment, and have mean zero.

%  Base case -- confirming p-values are in the right ballpark
tau = 1 ; % 1.96*se_analytic ; 
y = tau * t + e ; % outcomes determined

%  Data as table 
data = array2table([y,t,e],'VariableNames',{'y','t','e'}); 
data.constant = ones(N,1); 
mdl = fitlm(data,'y ~ t')

%  Potential randomizations 
T0 = double(tiedrank(rand(N,R))/N > ptreat ) ; % tiedrank operates within columns as required

%  Comparison across alternative treatment magnitudes 
TAU = transpose([1.02:0.005:1.05]) ; % transpose([-0.5:.1:0.5]*se_analytic + tau) ; % set of possible values to consider 

RESULTS = array2table([TAU,NaN(length(TAU),5)],'VariableNames',{'tau','tau_hat','p_analytic', 'p_ri', 'ci_lower', 'ci_upper' }); % storing analytic results

%  Demo:  Run regression and check analytical standard errors for closeness to p-value 

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
 
for tt = 1:length(TAU) 
	
	tic
	
	y = TAU(tt)*t + e ;  
	data.y = y ; 

	%  Regression estimates and analytic standard errors 
	mdl = fitlm(data,'y ~ t') 
	RESULTS(tt,'tau_hat') = mdl.Coefficients('t','Estimate') ;
	RESULTS(tt,'p_analytic') = mdl.Coefficients('t','pValue') ; 

	%  RI p-values 
	[b,~,pval ] = ri_ci( ...  %,pval, ci 
		data ...
		, 'y' ...
		, 't' ...
		, 'T0', T0 ...
		, 'P', R ...
		, 'TestZero', true ...
		, 'FindCI', false ...
		, 'RunParallel', RunParallel ...
		, 'Noisily', Noisily ... 
		, 'MaxQueries', 100 ... 
		)
	RESULTS(tt,'p_ri') = array2table(pval) ; 

	%  RI confidence intervals: initial guess, low number of repetition
	[~,~,~,ci_init] = ri_ci( ...
		data ...
		, 'y' ...
		, 't' ...
		, 'T0', T0 ...
		, 'P', 100 ...
		, 'TestZero', false ...
		, 'FindCI', true ...
		, 'RunParallel', RunParallel ...
		, 'MaxQueries', 25 ... 
		, 'MinStepSize', 0.01 ...
		)

	% Formulate guess based on small test run 
	CIguess = [ ...
		ci_init(1) + 0.5*(ci_init(1)-b) ... 
		ci_init(1) - 0.25*(ci_init(1)-b) ... 
		ci_init(2) - 0.25*(ci_init(2)-b) ... 
		ci_init(2) + 0.5*(ci_init(2)-b) ... 
		] 

	%  Proper CI run over smaller search space 
	[~,~,~,ci] = ri_ci( ...
		data ...
		, 'y' ...
		, 't' ...
		, 'T0', T0 ...
		, 'P', R ...
		, 'TestZero', false ...
		, 'FindCI', true ...
		, 'CIguess', CIguess ...
		, 'RunParallel', RunParallel ...
		, 'MaxQueries', 100 ... 
		, 'MinStepSize', 0.01 ...
		)
	
	RESULTS(tt,'ci_lower') = array2table(ci(1)); 
	RESULTS(tt,'ci_upper') = array2table(ci(2)); 
	
	toc

end 
 

%  Visualize correspondence between estimated treatment effects and p-values
figure(1)
clf 
hold on 
plot(RESULTS.tau,RESULTS.p_analytic, 'DisplayName', 'Analytic')
plot(RESULTS.tau,RESULTS.p_ri, 'DisplayName', 'RI') 
yline(0.05) ;
xlabel('Treatment effect (DGP)') 
ylabel('p-value') 
legend('Analytic', 'RI')  
hold off 


%  Visualize correspondence between estimated p-values and CI lower bounds (note tau > 0 so ignoring upper bound for now). 
figure(2) 
clf 
plot(RESULTS.p_ri, RESULTS.tau_hat)
hold on 
plot(RESULTS.p_ri, RESULTS.ci_lower)
plot(RESULTS.p_ri, RESULTS.ci_upper)
scatter(RESULTS.p_ri, RESULTS.ci_lower)  
scatter(RESULTS.p_ri, RESULTS.ci_upper)
patch([RESULTS.p_ri', fliplr(RESULTS.p_ri')], [RESULTS.ci_lower', fliplr(RESULTS.ci_upper')], 'b','FaceAlpha',0.2)   
yline(0) ;
xline(0.05) ; 
xlabel('RI p-value') 
ylabel('Estimates') 
hold off 
