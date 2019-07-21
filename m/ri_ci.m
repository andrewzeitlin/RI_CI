function varargout = ri_ci(DATA, outcome, txvars, T0, P, varargin) % model, stat, varargin

	%  Function to conduct RI, including (optionally) confidence intervals and test of the no-effect null.

	%  TODO (0): add point estimates to set of outputs, as default.
	%  TODO (1): add functionality to allow additional estimation commands:
	% 	- rereg()
	%   - fitlme()
	%	- clusterreg() 
	%
	%  TODO (2):  add plug-in principle functionality to KS test --  DGP

	%  Container for outputs:
	varargout = cell(1, nargout);
	
	%  Parse inputs
	params = inputParser ; 
	addOptional(params,'Controls',{}); % observational RHS variables
	addOptional(params,'Model','lm');  			% model type
	addOptional(params,'TestType','tStat');		%  name of test statistic
	addOptional(params,'TestSide','twosided'); 	% test type for the primary test
	addOptional(params,'TestZero',true); 			% to add a test of the zero null. Defaults on.
	addOptional(params,'FindCI',false);  			%  Switch: locate confidence interval
	addOptional(params,'CIguess',[0 0 0 0 ]);		% initial guess for confidence interval. Specified as a range around the lower bound (1st and second arguments) and a range around the upper bound (2nd and third arguments)
	addOptional(params,'CiSearchSize', 20) ; 		% default SDs of search width on EITHER side of the parameter estimate to search for CI (if CIguess not specified). 
	addOptional(params,'MaxQueries',10); 			% maximum number of trials to find each end of the confidence interval
	addOptional(params,'MinStepSize',0); 			% minimum step size for stopping rule.
	addOptional(params,'SignificanceLevel',0.05); 	% alpha for CI
	addOptional(params,'ShowMainEstimates',false); 	% to replay results of primary estimate.
	addOptional(params,'CheckBoundaries',true); 	% check boundaries for CI esitmation.
	addOptional(params,'GroupVar',{''}); 			% group variable for random effects or clustered estimates
	addOptional(params,'TheTx',{}); 				% for vector-valued treatments, this cell array contains the name of the treatment of interest.
	addOptional(params,'tau0',[] ) ; 				%  null to test if reporting a p-value. 
	addOptional(params,'PlugIn',true); 		%  use plug-in estimate of nuisance variables for test of zero null and construction of CIs. Overrides specified tau0 for nuisance treatments
	addOptional(params,'RunParallel',false) ;  % RI with parfor loops in ri_estimates.
	addOptional(params,'Support',[-inf, inf]); % boundary values for potential outcomes.
	addOptional(params,'Noisily',false); 
	parse(params,varargin{:}); 

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
	groupvar = params.Results.GroupVar; 
	tau0 = params.Results.tau0 ; 
	RunParallel = params.Results.RunParallel; 
	Noisily = params.Results.Noisily; 
	CiSearchSize = params.Results.CiSearchSize ; 

	%  If treatment assignment is vector-valued, need to decide which variable will be basis for the test for CI purposes.
	if length(TheTx) == 0 
		TheTx = txvars(1) ; % if no treatment variable specified, assume this is the first entryy in txvars
	end
	tx = table2array(DATA(:,TheTx)) ; 
	if size(tau0,1) > size(tau0,2) 
		error('null treatment vector tau0 should be (1 x K) not (K x 1)');
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
	
	%  If control set not empty, check for categorical variables, and expand to a set of dummies. 
	%  Then replace variable names in Controls
	if length(Controls) > 0 
		sprintf('Checking controls for categorical variables')
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

	%  Containers for results 
	TEST0 = NaN(P,length(txvars)); % to hold null distribution for test statistic.

	%  Estimate the model using the actual assignment. Obtain test statistic.  
	if strcmp(model,'lm')
		lm = fitlm(DATA(:,[txvars Controls outcome])) 
		TEST1 = table2array(lm.Coefficients([TheTx],[TestType]))'; 

		%  Inputs into starting values for search for 95% CI
		if FindCI 
			beta = table2array(lm.Coefficients(TheTx,'Estimate')); 
			se   = table2array(lm.Coefficients(TheTx,'SE')); 
		end

		%  Coefficients on nuisance treatments for default behavior of p-values and CIs
		if length(txvars) > 1 
			b_nuisance = table2array(lm.Coefficients(nuisanceTx,'Estimate')); 
		end

	elseif strcmp(model,'re') 
		sprintf('Now estimating RE model')
		result = rereg(DATA,outcome,[txvars Controls] ,groupvar )
		TEST1 = table2array(result([txvars],[TestType]));
		if FindCI 
			beta = table2array(result(TheTx,{'beta'}));
			se = table2array(result(TheTx,{'SE'}));
		end
		if length(txvars) > 1 
			b_nuisance = table2array(lm.Coefficients(nuisanceTx,'Estimate'))
		end
	elseif strcmp(model,'ks')
		% TODO:  Allow residualization of outcome at this stage, then used residualized outcome for all subsequent analysis.
		if length(Controls)> 0 | length(txvars) > 1
			y = table2array(DATA(:,outcome)) ; 
			rhs = [ones(size(y)),table2array(DATA(:,[~strcmp(txvars,tx),Controls]))];
			bb = rhs\y ;        % regression coefficient
			ydd = y - rhs*bb ;  % residuals 
		else 
			ydd = table2array(DATA(:,outcome));
			ydd = ydd - mean(ydd); 
		end
		[~,~,TEST1] = kstest2(ydd(tx==0),ydd(tx==1)); % two-sample KS stat

		%  For the KS test with any subsequent testing
		%  replace outcome variable in DATA with residualized version
		%  and get rid of Controls 
		if TestZero || FindCI 
			DATA(:,outcome) = array2table(ydd);
			Controls = {}; %  These have already been used to residualize.
		end

		if FindCI 
			%  use simple linear regression to initialize search region.
			lm = fitlm(tx,ydd);
			beta = table2array(lm.Coefficients(TheTx,'Estimate')); 
			se   = table2array(lm.Coefficients(TheTx,'SE')); 
		end
	end

	%  Conduct RI on the model, using the point estimmate as a point of comparison. Two-sided test, null hypothesis as specified by user (parameter tau0). 
	if TestZero 
		if PlugIn && length(txvars) > 1
			tau0(~strcmp(txvars,TheTx)) = b_nuisance ; 
		end
		[ pvalue TEST0 ] = ri_estimates(DATA,outcome,txvars,tau0 ...
			, model,T0,P ...
			, 'TheTx', TheTx ...
			, 'Controls', Controls ... 
			, 'TestType', TestType ...
			, 'TestSide',TestSide ...
			, 'GroupVar',groupvar ...
			, 'RunParallel', RunParallel ...
			, 'ShowWaitBar', Noisily ... 
			, 'WaitMessage', ['Outcome ' outcome ': Now conducting test of zero null'] ... 
			) ; 
	end

	%----------------------------------------------------------------------%
	%  If requested, find confidence interval
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
					, 'ks' ...
					, T0, P ...
					, 'Controls', Controls ... 
					, 'TestSide','right' ...
					, 'RunParallel', RunParallel ...
					) ; 
				p_lb = ri_estimates( ...
					DATA,outcome,txvars,tau_prime_lb ...  %  TODO: Update for vector-valued treatments with KS
					, 'ks' ...
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
					, 'GroupVar',groupvar ...
					); 
				p_lb = ri_estimates( ...
					DATA, outcome,txvars,tau_prime_lb ...
					, model ...
					, T0, P ...
					, 'TheTx', TheTx ... 
					, 'Controls', Controls ... 
					, 'TestType', TestType ... 
					, 'TestSide', 'twosided' ... % change back to 'left'
					, 'GroupVar',groupvar ...
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
			rectangle('Position',[lb, 0, figwidth, 0.025],'EdgeColor','none','FaceColor',[0.75 0.75 0.75]);
			line([beta beta] , [0 1],'LineStyle','--','Color','b') % vertical line at point estimate 
			xlabel('Coefficient estimate')
			ylabel('p-value') 
		end 

		%  CI.1.  Find upper boundary of CI. -------------------------------------%
		if mean(CIguess==0) == 1  % case where no search region has been specified
			lb = beta ;
			ub = beta + CiSearchSize*1.96*se ; 
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
					, 'ks' ...
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
					DATA, outcome,txvars,tau_prime, model ...
					,T0,P ... 
					, 'TheTx', TheTx ... 
					, 'Controls', Controls ... 
					, 'TestType', TestType ... 
					, 'TestSide', 'twosided' ... % 'right' ...
					, 'GroupVar',groupvar ...
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
			if p <= SignificanceLevel/2
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
		CI_UB = max(QUERIES_UB(QUERIES_UB(:,2) > SignificanceLevel/2, 1)) ; 

		%  Find lower boundary of CI. -------------------------------------%
		if mean(CIguess==0) == 1  % case where no search region has been specified
			ub = beta ;
			lb = beta - CiSearchSize*1.96*se ; 
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
					DATA,outcome,txvars,middle,Controls ...
					, 'ks' ...
					, T0, P ...
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

				p = ri_estimates(DATA, outcome,txvars,tau_prime, model,T0,P ...
					, 'TheTx', TheTx ... 
					, 'Controls', Controls ... 
					, 'TestType',TestType ...
					, 'TestSide','twosided' ...
					, 'GroupVar',groupvar ...
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
			if p <= SignificanceLevel/2
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
		CI_LB = min(QUERIES_LB( QUERIES_LB(:,2) > SignificanceLevel/2, 1)); 

		%  Release the figure 
		if Noisily
			% figure(figno)
			hold off ; 
		end 

		%  Return CI as a vector
		CI = [CI_LB , CI_UB];
	end


	% Write function outputs
	if TestZero , varargout{1} = pvalue ; end
	if FindCI , varargout{2} = CI ; end
	if nargout >= 3 , varargout{3} = TEST1 ; end 
	if nargout >= 4 , varargout{4} = TEST0 ; end 
	if nargout >= 5 , varargout{5} = QUERIES_UB ; end
	if nargout >= 6 , varargout{6} = QUERIES_LB ; end 
end
