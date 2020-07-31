function [beta, N, pvalue, CI, varargout] = ri_ci(DATA, outcome, txvars, varargin) % model, stat, varargin

	%  Function to conduct RI, including (optionally) confidence intervals and test of the no-effect null.

	%  TODO (0): add point estimates to set of outputs, as default.
	%  TODO (1): add functionality to allow additional estimation commands:
	% 	- rereg()
	%   - fitlme()
	%	- clusterreg() 
	%
	%  TODO (2):  add plug-in principle functionality to KS test --  DGP

	%  Container for outputs:
	varargout = cell(1, nargout - 4);
	
	%  Parse inputs
	params = inputParser ; 
	addOptional(params,'T0',{}); %  Array of potential treatments
	addOptional(params,'P',{}); %  Number of permutations of potential treatments to consider.
	addOptional(params,'Controls',{}); % observational RHS variables
	addOptional(params,'Model','lm');  			% model type
	addOptional(params,'TestType','tStat');		%  name of test statistic
	addOptional(params,'TestSide','twosided'); 	% test type for the primary test
	addOptional(params,'TestZero',true); 			% to add a test of the zero null. Defaults on.
	addOptional(params,'FindCI',false);  			%  Switch: locate confidence interval
	addOptional(params,'CIguess',[0 0 0 0 ]);		% initial guess for confidence interval. Specified as a range around the lower bound (1st and second arguments) and a range around the upper bound (2nd and third arguments)
	addOptional(params,'CiSearchSize', 20) ; 		% default SDs of search width on EITHER side of the parameter estimate to search for CI (if CIguess not specified). 
	addOptional(params,'MaxQueries',100); 			% maximum number of trials to find each end of the confidence interval
	addOptional(params,'MinStepSize',0); 			% minimum step size for stopping rule.
	addOptional(params,'SignificanceLevel',0.05); 	% alpha for CI
	addOptional(params,'ShowMainEstimates',false); 	% to replay results of primary estimate.
	addOptional(params,'CheckBoundaries',true); 	% check boundaries for CI esitmation.
	addOptional(params,'GroupVar',{}); 			% group variable for random effects or clustered estimates
	addOptional(params,'WeightVar',{}); 			% group variable for random effects or clustered estimates
	addOptional(params,'TheTx',{}); 				% for vector-valued treatments, this cell array contains the name of the treatment of interest.
	addOptional(params,'tau0',[] ) ; 				%  null to test if reporting a p-value. 
	addOptional(params,'PlugIn',true); 		%  use plug-in estimate of nuisance variables for test of zero null and construction of CIs. Overrides specified tau0 for nuisance treatments
	addOptional(params,'RunParallel',false) ;  % RI with parfor loops in ri_estimates.
	addOptional(params,'Support',[-inf, inf]); % boundary values for potential outcomes.
	addOptional(params,'Noisily',false); 
	
	parse(params,varargin{:}); 

	T0 = params.Results.T0; 
	P = params.Results.P; 

	model = {params.Results.Model}; 
	TestType = {params.Results.TestType}; 
	TestSide = {params.Results.TestSide};
	Controls = params.Results.Controls;
	Support = sort(params.Results.Support); 
	TheTx = params.Results.TheTx ; 

	TestZero = params.Results.TestZero; 
	FindCI = params.Results.FindCI ;
	CIguess = sort(params.Results.CIguess); % forcing ascending order.
	PlugIn = params.Results.PlugIn; 

	MaxQueries = params.Results.MaxQueries; 
	MinStepSize = params.Results.MinStepSize; 
	SignificanceLevel = params.Results.SignificanceLevel;  
	GroupVar = params.Results.GroupVar; 
	WeightVar = params.Results.WeightVar; 
	tau0 = params.Results.tau0 ; 
	RunParallel = params.Results.RunParallel; 
	Noisily = params.Results.Noisily; 
	CiSearchSize = params.Results.CiSearchSize ; 

	%  If txvars or TheTx is specified as a character vector, convert to cell array.
	if ischar(txvars), txvars = {txvars}; end 
	if ischar(TheTx), TheTx = {TheTx}; end 

	%  If treatment assignment is vector-valued, need to decide which variable will be basis for the test for CI purposes.
	if length(TheTx) == 0 
		TheTx = txvars(1) ; % if no treatment variable specified, assume this is the first entryy in txvars
	end
	%  Require GroupVar to be specified for re model
	if (strcmp(model,'re') || strcmp(model,'lme')) && isempty(GroupVar)
		error('Random effects model requires option GroupVar to be specified.')
	end
	%  Populate a weight variable (as constants) if none exists.
	if strcmp(model,'lm') && isempty(WeightVar)
		DATA.uniform_weight_vector = ones(size(DATA,1),1);
		WeightVar = 'uniform_weight_vector' ; 
	end

	%  Require parameter TestSide set to right for KS test
	if strcmp(model,'ks') && ~strcmp(TestSide,'right')
		error('For the KS test, must specify a right-tailed p-value')
	end
	%  Default null hypothesis tau0 = 0
	if length(tau0) == 0 
		tau0 = zeros(length(txvars),1);
	end
	%  List of nuisance treatment variables
	if length(txvars) > 1 
		nuisanceTx = txvars(~strcmp(txvars,TheTx)) ; 
	else
		nuisanceTx = {}; 
	end

	%  CHECK FOR MISSINGS IN ANALYTIC VARIABLES. 
	% IF ANY, REMOVE FROM DATA AND POTENTIAL ASSIGNMENTS (as appropriate) 
    if sum(max(ismissing(DATA(:,[outcome txvars Controls GroupVar WeightVar])),[],2)) > 0 
        tokeep = min( ~ismissing(DATA(:,[outcome txvars Controls GroupVar WeightVar])), [], 2);
        DATA = DATA(tokeep,:) ; 
        if TestZero || FindCI % either of these requires the set of potential assignments
        	T0 = T0(tokeep,:,:) ;
        end
    end

    %  Extract number of observations.
    N = size(DATA,1); 
	
	%  If control set not empty, check for categorical variables, and expand to a set of dummies. 
	%  Then replace variable names in Controls
	if length(Controls) > 0 
		Controls_Temp = Controls ; % making changes in a temporary list, so that substitutions don't change looping over control variiables.
		for k = 1:length(Controls)
			xvar = Controls(k);
			if iscategorical(DATA{:,xvar})
				Dx = array2table(dummyvar(DATA{:,xvar})); % dummy variables 
				for ll = 1:size(Dx,2)
					Dx.Properties.VariableNames(ll) = strcat(xvar, '_', num2str(ll));
				end
				Dx = Dx(:,2:end);  %  Leave one category out.
				DATA = [DATA,Dx];  %  Append these indicators to the table.

				%  Now remove xvar from Controls and add VariableNames from table Dx instead.
				Controls_Temp = Controls_Temp(~strcmp(Controls_Temp,xvar));
				Controls_Temp = [Controls_Temp, Dx.Properties.VariableNames];
			end
		end
		Controls = Controls_Temp; 
	end

	%  Estimate the model using the actual assignment.  Save point estimate as beta; obtain test statistic.  
	if strcmp(model,'lm')
		lm = fitlm(DATA(:,[txvars Controls outcome]),'Weights',table2array(DATA(:,WeightVar))) ; 
		% TEST1 = table2array(lm.Coefficients([TheTx],[TestType]))'; 
		if Noisily  % under Noisily mode, replay the model 
			lm
		end 
		estimates = lm.Coefficients; 
		%  Check LM results for indicator variables where _1 suffix has been added; strip this.
		indicators = find(contains(estimates.Properties.RowNames,'_1'))' ; 
		for ii = indicators 
			%  Check if abbreviated form exists in list of covariates in model. If so, strip suffix.
			if sum(strcmp(strrep(estimates.Properties.RowNames(ii),'_1',''), txvars )) > 0 
				estimates.Properties.RowNames(ii) = strrep(estimates.Properties.RowNames(ii),'_1','') ; 
			end
		end

		beta = table2array(estimates(TheTx,'Estimate')); 
		%  Inputs into starting values for search for 95% CI
		if FindCI 
			se   = table2array(estimates(TheTx,'SE')); 
		end
		%  Coefficients on nuisance treatments for default behavior of p-values and CIs
		if length(txvars) > 1 
			b_nuisance = table2array(estimates(nuisanceTx,'Estimate')); 
		end

	elseif strcmp(model,'lme')
		lme = fitlmematrix( ...
			[ones(size(DATA,1),1), table2array(DATA(:, [txvars, Controls]))] ...  % FE design matrix
			, table2array(DATA(:,outcome)) ...  % outcome
			, ones(size(DATA,1),1) ... % RE design matrix
			, table2array(DATA(:, GroupVar)) ... % grouping variable(s) 
			, 'FixedEffectPredictors', [ 'Intercept' txvars Controls ] ...
			, 'RandomEffectPredictors', GroupVar ...
			) ; 

		%  Extract coefficients 
		%  Extracting parameter estimates and t statistics 
		[~,~,estimates] = fixedEffects(lme); 
		estimates = dataset2table(estimates);
		%  Remove the '_1' that fitlme() adds to indicator variable names, so that they match the supplied parameters exactly, if those parameters match the data.
		indicators = find(contains(estimates.Name,'_1'))' ; 
		for ii = indicators 
			estimates.Name(ii) = strrep(estimates.Name(ii),'_1','');
		end
		estimates.Properties.RowNames = estimates.Name; % add row names to make lookups faster & more transparent

		beta = table2array(estimates(TheTx,'Estimate')); 

		%  Inputs into starting values for search for 95% CI
		if FindCI 
			se   = table2array(estimates(TheTx,'SE')); 
		end

		%  Coefficients on nuisance treatments for default behavior of p-values and CIs
		if length(txvars) > 1 
			for kk = 1:length(nuisanceTx)
				b_nuisance(kk) = table2array(estimates(nuisanceTx{kk},'Estimate')); 
			end
		end
		if Noisily  % under Noisily mode, replay the model 
			lme   
		end 

	elseif strcmp(model,'re') 
		re_model = rereg(DATA,outcome,[txvars Controls] ,GroupVar ) ; 
		% TEST1 = table2array(re_model([TheTx],[TestType]));
		beta = table2array(re_model(TheTx,{'Estimate'}));
		if FindCI 
			se = table2array(re_model(TheTx,{'SE'}));
		end
		if length(txvars) > 1 
			b_nuisance = table2array(re_model(nuisanceTx,'Estimate')) ; 
		end
		if Noisily 
			re_model 
		end 
	elseif strcmp(model,'ks')
		if length(Controls)> 0 
			y = table2array(DATA(:,outcome)) ; 
			rhs = ones(size(y)) ; % initializing the covariate set 
			%  Residualize controls based on treatments and then add them to the covariate set.
			itx = [ones(size(y)),DATA(:,txvars)] ; % treatments and a constant. 
			for k = 1:length(Controls) 
				x = DATA(:,Controls(k)); 
				bx = itx \ x;  % residualizing control x
				xdd = x - itx*bx ; 
				rhs = [rhs, xdd] ; 
			end
			%  Now, residualize outcome wrt controls 
			bb = rhs\y ;        % regression coefficient
			ydd = y - rhs*bb ;  % residuals 
		else 
			ydd = table2array(DATA(:,outcome));
			ydd = ydd - mean(ydd); 
		end

		%  Replace outcome variable in DATA with residualized version
		%  and get rid of Controls 
		DATA(:,outcome) = array2table(ydd);
		Controls = {}; %  These have already been used to residualize.

		%  use simple linear regression to initialize search region, return point estimate as <beta>
		lm = fitlm(DATA(:,[txvars outcome])); 
		estimates = lm.Coefficients; 
		%  Check LM results for indicator variables where _1 suffix has been added; strip this.
		indicators = find(contains(estimates.Properties.RowNames,'_1'))' ; 
		for ii = indicators 
			%  Check if abbreviated form exists in list of covariates in model. If so, strip suffix.
			if sum(strcmp(strrep(estimates.Properties.RowNames(ii),'_1',''), txvars )) > 0 
				estimates.Properties.RowNames(ii) = strrep(estimates.Properties.RowNames(ii),'_1','') ; 
			end
		end

		%  Extract model results of interest
		beta = table2array(estimates(TheTx,'Estimate')); 
		%  Inputs into starting values for search for 95% CI
		if FindCI 
			se   = table2array(estimates(TheTx,'SE')); 
		end
		%  Coefficients on nuisance treatments for default behavior of p-values and CIs
		if length(txvars) > 1 
			b_nuisance = table2array(estimates(nuisanceTx,'Estimate')); 
		end
	end

	%----------------------------------------------------------------------%
	%--  P-VALUES  --------------------------------------------------------%
	%----------------------------------------------------------------------%

	%  Conduct RI on the model, using the point estimmate as a point of comparison. Two-sided test, null hypothesis as specified by user (parameter tau0). 
	if TestZero 
		if PlugIn && length(txvars) > 1
			tau0(~strcmp(txvars,TheTx)) = b_nuisance ; 
		end
		[ pvalue TEST0 TEST1 ] = ri_estimates( ...
			DATA ...
			, outcome ...
			, txvars ...
			, tau0 ...
			, model ...
			, T0 ...
			, P ...
			, 'TheTx', TheTx ...
			, 'Controls', Controls ... 
			, 'TestType', TestType ...
			, 'TestSide',TestSide ...
			, 'GroupVar', GroupVar ...
			, 'WeightVar', WeightVar ... 
			, 'RunParallel', RunParallel ...
			, 'ShowWaitBar', Noisily ... 
			, 'WaitMessage', ['Outcome ' outcome ': Now conducting test of zero null'] ... 
			) ; 
	else 
		pvalue = {}; 
	end

	%----------------------------------------------------------------------%
	%--  CONFIDENCE INTERVALS  --------------------------------------------%
	%----------------------------------------------------------------------%
	if FindCI 

		%  Confirm that p-value at boundaries of search region are below significance threshold
		if params.Results.CheckBoundaries 
			if mean(CIguess==0) == 1  % case where no search region has been specified
				lb = beta - CiSearchSize*1.96*se ; % 
				ub = beta + CiSearchSize*1.96*se ; 
			else 
				lb = CIguess(1);
				ub = CIguess(4); 
			end

			tau_prime = NaN(length(txvars),1); 
			if length(txvars) > 1
				if PlugIn 
					tau_prime(find(~strcmp(txvars,TheTx))) = b_nuisance ; 
				else 
					tau_prime = tau0; % allow manually specifying values for nuisance parameters
				end
			end
			tau_prime_ub = tau_prime ;
			tau_prime_ub(find(strcmp(txvars,TheTx))) = ub; 
			tau_prime_lb = tau_prime ; 
			tau_prime_lb(find(strcmp(txvars,TheTx))) = lb; 

			if strcmp(model,'ks')
				p_ub = ri_estimates( ...
					DATA,outcome,txvars, tau_prime_ub ...  % TODO: Update for vector-valued treatments with KS
					, model ...
					, T0, P ...
					, 'Controls', Controls ... 
					, 'TestSide','right' ...
					, 'RunParallel', RunParallel ...
					) ; 
				p_lb = ri_estimates( ...
					DATA,outcome,txvars,tau_prime_lb ...  %  TODO: Update for vector-valued treatments with KS
					, model ...
					, T0, P ...
					, 'Controls', Controls ... 
					, 'TestSide','right' ...
					, 'RunParallel', RunParallel ...
					) ; 
			else
				p_ub = ri_estimates( ...
					DATA, outcome,txvars,tau_prime_ub ...
					, model ...
					, T0, P  ...
					, 'TheTx', TheTx ... 
					, 'Controls', Controls ... 
					, 'TestType', TestType ...
					, 'TestSide' , 'twosided' ... 
					, 'RunParallel', RunParallel ...
					, 'GroupVar',GroupVar ...
					, 'WeightVar', WeightVar ...
					); 
				p_lb = ri_estimates( ...
					DATA, outcome,txvars,tau_prime_lb ...
					, model ...
					, T0, P ...
					, 'TheTx', TheTx ... 
					, 'Controls', Controls ... 
					, 'TestType', TestType ... 
					, 'TestSide', 'twosided' ... % change back to 'left'
					, 'GroupVar',GroupVar ...
					, 'WeightVar', WeightVar ...
					, 'RunParallel', RunParallel ...
					);
			end
			if p_ub > SignificanceLevel/2
				%sprintf('You had a test value of %0.2f', TEST1)
				%ksdensity(TEST0)
				sprintf('Initial attempt at upper bound %0.2f yielded p-value of %0.12f',ub,p_ub)
				error('Initial value for upper bound of CI not high enough.')
			end
			if p_lb > SignificanceLevel/2 
				sprintf('Initial attempt at lower bound %0.2f yielded p-value of %0.12f',lb,p_lb)
				error('Initial value for lower bound of CI not low enough.')
			end 				
		end

		%  Set up diagnostic figure.  
		if Noisily 
			figure; clf ;
			figno = get(gcf,'Number') ; % store figure numer 
			axis([lb ub 0 0.25]);
			hold on 
			figwidth = ub - lb ;
			rectangle('Position',[lb, 0, figwidth, 0.05],'EdgeColor','none','FaceColor',[0.75 0.75 0.75]);
			line([beta beta] , [0 1],'LineStyle','--','Color','b') % vertical line at point estimate 
			xlabel('Coefficient estimate')
			ylabel('Two-sided p-value') 
		end 

		%  CI.1.  Find upper boundary of CI. -------------------------------------%
		if mean(CIguess==0) == 1  % case where no search region has been specified
			lb = beta ;
			ub = beta + CiSearchSize*1*se ; 
		else 
			lb = CIguess(3);
			ub = CIguess(4); 
		end
		middle = (lb + ub) / 2 ; 

		tau_prime = NaN(length(txvars),1);  %  initializing. This will always hold the current trial value.

		QUERIES_UB = NaN(MaxQueries,2); 
		q = 1 ; 
		stepsize = MinStepSize + 1 ;
		while q <= MaxQueries & stepsize >= MinStepSize % TODO: add step-size constraint.
			% fprintf('This is counter number %i \n', q)
			QUERIES_UB(q,1) = middle;
			if strcmp(model,'ks')
				p = ri_estimates( ...
					DATA,outcome,txvars,middle ...
					, model ...
					, T0, P ...
					, 'Controls', Controls ... 
					, 'TestSide','right' ...
					, 'RunParallel', RunParallel ...
					) ; 
			else 
				if length(txvars) > 1
					if PlugIn 
						tau_prime(find(~strcmp(txvars,TheTx))) = b_nuisance ; 
					else 
						tau_prime = tau0;
					end
				end
				tau_prime(find(strcmp(txvars,TheTx))) = middle; 

				p = ri_estimates( ...
					DATA ...
					, outcome ...
					, txvars ...
					, tau_prime ...
					, model ...
					, T0 ...
					, P ... 
					, 'TheTx', TheTx ... 
					, 'Controls', Controls ... 
					, 'TestType', TestType ... 
					, 'TestSide', 'twosided' ... % 
					, 'GroupVar',GroupVar ...
					, 'WeightVar', WeightVar ...
					, 'RunParallel', RunParallel ...
					); 
			end
			QUERIES_UB(q,2) = p;

			if Noisily
				% figure(figno)
				h1=scatter(QUERIES_UB(:,1),QUERIES_UB(:,2),'MarkerEdgeColor','b');
			end 

			%  If reject, move left.  Otherwise, move right.
			if q > 1, m0 = middle ; end % store previous tested value if appropriate
			if p <= SignificanceLevel  % /2  Two-sided test doubles the p-value in the relevant tail
				ub = middle;
				middle = (ub + lb ) / 2;  
			else
				lb = middle;
				middle = (ub + lb ) / 2;  
			end
			if q > 1 
				stepsize = abs(m0 - middle) ;
			else 
				stepsize = MinStepSize + 1 ; % don't let this constraint bind on first try
			end 

			%  Update counter
			q = q+1;
		end
		%  Remove unused entries from queries table
		QUERIES_UB = QUERIES_UB(~isnan(QUERIES_UB(:,1)),:);
		CI_UB = max(QUERIES_UB(QUERIES_UB(:,2) >= SignificanceLevel, 1)); %SignificanceLevel/2, 1)) ;   <- Two-sided test doubles the p-value in the relevant tail

		%  Find lower boundary of CI. -------------------------------------%
		if mean(CIguess==0) == 1  % case where no search region has been specified
			ub = beta ;
			lb = beta - CiSearchSize*1*se ; 
		else 
			ub = CIguess(2);
			lb = CIguess(1);
		end
		middle = (lb + ub) / 2 ; 

		QUERIES_LB = NaN(MaxQueries,2); 
		q = 1 ; 
		stepsize = MinStepSize + 1 ;
		while q <= MaxQueries & stepsize >= MinStepSize 
			%  fprintf('This is counter number %i \n', q)
			QUERIES_LB(q,1) = middle;
			if strcmp(model,'ks')
				p = ri_estimates( ...
					DATA,outcome,txvars,middle ...
					, model ...
					, T0 ...
					, P ...
					, 'Controls', Controls ... 
					, 'TestSide','right' ...
					, 'RunParallel', RunParallel ...
					) ; 
			else
				if length(txvars) > 1
					if PlugIn 
						tau_prime(find(~strcmp(txvars,TheTx))) = b_nuisance ; 
					else 
						tau_prime = tau0;
					end
				end
				tau_prime(find(strcmp(txvars,TheTx))) = middle; 

				p = ri_estimates( ...
					DATA ...
					, outcome ...
					, txvars ...
					, tau_prime ...
					, model ...
					, T0 ...
					, P ...
					, 'TheTx', TheTx ... 
					, 'Controls', Controls ... 
					, 'TestType',TestType ...
					, 'TestSide', 'twosided' ... % 
					, 'GroupVar',GroupVar ...
					, 'WeightVar', WeightVar ...
					, 'RunParallel', RunParallel ...
					); 
			end 
			QUERIES_LB(q,2) = p;
			if Noisily
				% figure(figno) 
				h2=scatter(QUERIES_LB(:,1),QUERIES_LB(:,2),'MarkerEdgeColor', 'b'); 
			end 


			%  If reject, move left.  Otherwise, move right.
			if q > 1, m0 = middle ; end % store previous tested value if appropriate
			if p <= SignificanceLevel  %/2    <- Two-sided test doubles the p-value in the relevant tail
				lb = middle;
				middle = (ub + lb ) / 2;  
			else
				ub = middle;
				middle = (ub + lb ) / 2;  
			end

			%  Update counter(s)
			if q > 1 
				stepsize = abs(m0 - middle) ;
			else 
				stepsize = MinStepSize + 1 ; % don't let this constraint bind on first query
			end 
			q = q+1;
		end
		QUERIES_LB = QUERIES_LB(~isnan(QUERIES_LB(:,1)),:);
		CI_LB = min(QUERIES_LB( QUERIES_LB(:,2) >= SignificanceLevel, 1 )); % SignificanceLevel/2, 1));   <- Two-sided test doubles the p-value in the relevant tail

		%  Release the figure 
		if Noisily
			% figure(figno)
			hold off ; 
		end 

		%  Return CI as a vector
		CI = [CI_LB , CI_UB];
	else 
		CI = {} ; % to avoid needing to change output arguments when CI is not requested
	end

	%--------------------------------------------------------------------------%
	%%  WRITE FUNCTION OUTPUTS 
	%--------------------------------------------------------------------------%
	% Write function outputs
	% beta 		<- already defined above.
	% N 		<- already defined above.
	% pvalue 	<- defined above
	% CI 		<- defined above 
	if length(varargout) >= 1 , varargout{1} = TEST1 ; end 
	if length(varargout) >= 2 , varargout{2} = TEST0 ; end 
	if length(varargout) >= 3 , varargout{3} = QUERIES_UB ; end
	if length(varargout) >= 4 , varargout{4} = QUERIES_LB ; end 
end
