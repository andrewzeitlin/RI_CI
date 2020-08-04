


%--------------------------------------------------------------------------%
%%  Preliminaries 
%--------------------------------------------------------------------------%

%  Clear graphics 

%  Set seed for replicability 
rng('default')

%  File paths 
clear ri_ci 
clear ri_estimates 
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
N = 30000;
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

%%  Demo:  Run regression and check analytical standard errors for closeness to p-value 
%  Analytic benchmark

%  First pass with ri_ci() 
tic 
[b,N,pval] = ri_ci( ...
	data ...
	, 'y' ...
	, 't' ...
	, 'T0', T0 ... 
	, 'P', 1000 ...
	, 'TestZero',true ...
	, 'FindCI', false ...
	, 'TestSide', 'twosided' ... 
	, 'tau0', [0] ... 
	, 'Noisily', true ...
	) 

toc 