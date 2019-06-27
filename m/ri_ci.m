function varargout = ri_ci(DATA, outcome, txvars, T0, P, varargin) % model, stat, varargin

	%  Function to conduct RI, including (optionally) confidence intervaals and test of the no-effect null.

	%  TODO (1): add functionality to allow additional estimation commands:
	%	- clusterreg() 
	%   - fitlme()
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
	addOptional(params,'MaxQueries',10); 			% maximum number of trials to find each end of the confidence interval
	addOptional(params,'MinStepSize',0); 			% minimum step size for stopping rule.
	addOptional(params,'SignificanceLevel',0.05); 	% alpha for CI
	addOptional(params,'ShowMainEstimates',false); 	% to replay results of primary estimate.
	addOptional(params,'CheckBoundaries',true); 	% check boundaries for CI esitmation.
	addOptional(params,'GroupVar',{''}); 			% group variable for random effects or clustered estimates
	addOptional(params,'TheTx',{}); 				% for vector-valued treatments, this cell array contains the name of the treatment of interest.
	addOptional(params,'tau0',[] ) ; 				%  null to test if reporting a p-value. 
	addOptional(params,'PlugIn',true); 		%  use plug-in estimate of nuisance variables for test of zero null and construction of CIs. Overrides specified tau0 for nuisance treatments
	addOptional(params,'Noisily',false); 
	parse(params,varargin{:}); 

	model = {params.Results.Model}; 
	TestType = {params.Results.TestType}; 
	TestSide = {params.Results.TestSide};
	xvars = params.Results.Controls;
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
	Noisily = params.Results.Noisily; 

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
	

	%  Containers for results 
	TEST0 = NaN(P,length(txvars)); % to hold null distribution for test statistic.

	%  Estimate the model using the actual assignment. Obtain test statistic.  
	if strcmp(model,'lm')
		lm = fitlm(DATA(:,[txvars xvars outcome])) 
		TEST1 = table2array(lm.Coefficients([TheTx],[TestType]))'

		%  Starting values for search for 95% CI
		if FindCI 
			beta = table2array(lm.Coefficients(TheTx,'Estimate')); 
			se   = table2array(lm.Coefficients(TheTx,'SE')); 
		end

		%  Coefficients on nuisance treatments for default behavior of p-values and CIs
		if length(txvars) > 1 
			b_nuisance = table2array(lm.Coefficients(nuisanceTx,'Estimate'))
		end

	elseif strcmp(model,'re') 
		sprintf('Now estimating RE model')
		result = rereg(DATA,outcome,[txvars xvars] ,groupvar )
		TEST1 = table2array(result([txvars],[TestType]));
		if FindCI 
			beta = table2array(result(txvars,{'beta'}));
			se = table2array(result(txvars,{'SE'}));
		end
		if length(txvars) > 1 
			b_nuisance = table2array(lm.Coefficients(nuisanceTx,'Estimate'))
		end
	elseif strcmp(model,'ks')
		% TODO:  Allow residualization of outcome at this stage, then used residualized outcome for all subsequent analysis.
		if length(xvars)> 0 | length(txvars) > 1
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
		%  and get rid of xvars 
		if TestZero || FindCI 
			DATA(:,outcome) = array2table(ydd);
			xvars = {}; %  These have already been used to residualize.
		end

		if FindCI 
			%  use simple linear regression to initialize search region.
			lm = fitlm(tx,ydd);
			beta = table2array(lm.Coefficients(2,'Estimate')); 
			se   = table2array(lm.Coefficients(2,'SE')); 
		end
	end

	%  Conduct RI on the model, using the point estimmate as a point of comparison. Two-sided test, null hypothesis as specified by user (parameter tau0). 
	if TestZero 
		if PlugIn && length(txvars) > 1
			tau0(~strcmp(txvars,TheTx)) = b_nuisance ; 
		end
		[ pvalue TEST0 ] = ri_estimates(DATA,outcome,txvars,tau0 ...
			,xvars, model,T0,P ...
			, 'TheTx', TheTx ...
			, 'TestType', TestType ...
			, 'TestSide',TestSide ...
			, 'TestValue', TEST1 ...
			, 'GroupVar',groupvar ...
			) ; 
	end



	if params.Results.ShowMainEstimates 
		sprintf('__RESULTS OF ANALYTICAL MODEL:__')
		lm 
		sprintf('The p-value from Randomization Inference for the hypothesis that tau = %02.2f is %0.2f',tau0, pvalue)
	end

	%----------------------------------------------------------------------%
	%  If requested, find confidence interval
	if FindCI 

		%  Find upper boundary of CI. -------------------------------------%
		if mean(CIguess==0) == 1  % case where no search region has been specified
			lb = beta ;
			ub = beta + 10*1.96*se ; 
		else 
			lb = CIguess(3);
			ub = CIguess(4); 
		end
		middle = (lb + ub) / 2 ; 

		tau_prime = NaN(length(txvars),1);  %  initializing. This will always hold the current trial value.

		%  Confirm that p-value at upper bound of search region is below significance threshold
		if params.Results.CheckBoundaries 
			if length(txvars) > 1
				if PlugIn 
					tau_prime(find(~strcmp(txvars,TheTx))) = b_nuisance ; 
				else 
					tau_prime = tau0;
				end
			end
			tau_prime(find(strcmp(txvars,TheTx))) = ub; 

			if strcmp(model,'ks')
				p = ri_estimates( ...
					DATA,outcome,txvars,tau_prime,xvars ...  % TODO: Confirm vector-valued treatment extension hasn't broken simple case
					, model ... % 'ks' ...
					, T0, P ...
					, 'TestSide','right' ...
					) ; 
			else
				[p] = ri_estimates( ...
					DATA, outcome,txvars,tau_prime,xvars ...
					, model ...
					, T0, P  ...
					, 'TheTx', TheTx ... 
					, 'TestType', TestType ...
					, 'TestSide' , 'twosided' ... 
					, 'TestValue', TEST1 ...
					, 'GroupVar',groupvar ...
					); 
			end
			if p > SignificanceLevel/2
				sprintf('You had a test value of %0.2f', TEST1)
				ksdensity(TEST0)
				sprintf('Initial attempt at upper bound %0.2f yielded p-value of %0.12f',ub,p)
				error('Initial value for upper bound of CI not high enough.')
			end
		end

		QUERIES_UB = NaN(MaxQueries,2); 
		q = 1 ; 
		stepsize = MinStepSize + 1 ;
		while q <= MaxQueries & stepsize >= MinStepSize % TODO: add step-size constraint.
			% fprintf('This is counter number %i \n', q)
			QUERIES_UB(q,1) = middle;
			if strcmp(model,'ks')
				p = ri_estimates( ...
					DATA,outcome,txvars,middle,xvars ...
					, 'ks' ...
					, T0, P ...
					, 'TestSide','right' ...
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

				[p , ~ , ~ ] = ri_estimates( ...
					DATA, outcome,txvars,tau_prime,xvars, model ...
					,T0,P ... 
					, 'TheTx', TheTx ... 
					, 'TestType', TestType ... 
					,'TestSide', 'twosided' ... % 'right' ...
					, 'TestValue',TEST1 ...
					,'GroupVar',groupvar ...
					); 
			end
			QUERIES_UB(q,2) = p;

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
		CI_UB = max(QUERIES_UB(QUERIES_UB(:,2) > SignificanceLevel/2, 1)) ; 

		%  Find lower boundary of CI. -------------------------------------%
		if mean(CIguess==0) == 1  % case where no search region has been specified
			ub = beta ;
			lb = beta - 10*1.96*se ; 
		else 
			ub = CIguess(2);
			lb = CIguess(1);
		end
		middle = (lb + ub) / 2 ; 

		%  Confirm that p-value at lower bound of search region is below significance threshold
		if params.Results.CheckBoundaries 
			if strcmp(model,'ks')
				p = ri_estimates( ...
					DATA,outcome,txvars,lb,xvars ...
					, 'ks' ...
					, T0, P ...
					, 'TestSide','right' ...
					) ; 
			else
				if length(txvars) > 1
					if PlugIn 
						tau_prime(find(~strcmp(txvars,TheTx))) = b_nuisance ; 
					else 
						tau_prime = tau0;
					end
				end
				tau_prime(find(strcmp(txvars,TheTx))) = lb; 
				[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,tau_prime, xvars, model,T0,P ...
					, 'TheTx', TheTx ... 
					, 'TestType', TestType ... 
					, 'TestSide', 'twosided' ... % change back to 'left'
					, 'TestValue',TEST1 ...
					, 'GroupVar',groupvar);
			end 
			if p > SignificanceLevel/2
				error('Initial value for lower bound of CI not low enough.')
			end
		end

		QUERIES_LB = NaN(MaxQueries,2); 
		q = 1 ; 
		stepsize = MinStepSize + 1 ;
		while q <= MaxQueries & stepsize >= MinStepSize 
			%  fprintf('This is counter number %i \n', q)
			QUERIES_LB(q,1) = middle;
			if strcmp(model,'ks')
				p = ri_estimates( ...
					DATA,outcome,txvars,middle,xvars ...
					, 'ks' ...
					, T0, P ...
					, 'TestSide','right' ...
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

				[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,tau_prime, xvars, model,T0,P ...
					, 'TheTx', TheTx ... 
					, 'TestType',TestType ...
					, 'TestSide','twosided' ...
					, 'TestValue',TEST1 ...
					, 'GroupVar',groupvar ...
					); 
			end 
			QUERIES_LB(q,2) = p;

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
		CI_LB = min(QUERIES_LB( QUERIES_LB(:,2) > SignificanceLevel/2, 1)); 

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
