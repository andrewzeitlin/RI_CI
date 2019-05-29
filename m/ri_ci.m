

function [pvalue CI QUERIES_UB QUERIES_LB TEST1 TEST0 y0 ] = ri_ci(DATA, outcome, txvars, tau0 , T0, P, varargin) % model, stat, varargin

	%  Function to conduct RI, including (optionally) confidence intervaals and test of the no-effect null.

	if size(tau0,1) > size(tau0,2) 
		error('null treatment vector tau0 should be (1 x K) not (K x 1)');
	end

	%  Parse inputs
	params = inputParser ; 
	addOptional(params,'Controls',{'unspecified'});  % observational RHS variables
	addOptional(params,'Clusters',{'unspecified'});  % for clustered standard errors in linear model
	addOptional(params,'Model',{'lm'});  % model type
	addOptional(params,'TestType',{'tStat'});
	addOptional(params,'TestSide',{'twosided'}); % test type for the primary test
	addOptional(params,'FindCI',false);  %  Switch: locate confidence interval
	addOptional(params,'CIguess',[0 0]); % initial guess for confidence interval. 
	addOptional(params,'MaxQueries',10); % maximum number of trials to find each end of the confidence interval
	addOptional(params,'MinStepSize',0); % minimum step size for stopping rule.
	addOptional(params,'SignificanceLevel',0.05); % alpha for CI
	parse(params,varargin{:}); 

	g = params.Results.Clusters;
	model = params.Results.Model; 
	TestType = params.Results.TestType; 
	TestSide = params.Results.TestSide;
	xvars = params.Results.Controls;
	if strcmp(xvars{1}, 'unspecified'), xvars = {}; end
	FindCI = params.Results.FindCI ;
	MaxQueries = params.Results.MaxQueries; 
	SignificanceLevel = params.Results.SignificanceLevel;  

	if FindCI & length(txvars)> 1
		error('Confidence interval search currently supports only one-dimensional treatment.')
	end


	%  Containers for results 
	TEST0 = NaN(P,length(txvars)); % to hold null distribution for test statistic.

	%  Estimate the model using the actual assignment
	sprintf('__RESULTS OF ANALYTICAL MODEL:__');
	if strcmp(model,'lm')
		lm = fitlm(DATA(:,[txvars xvars outcome]))
		TEST1 = table2array(lm.Coefficients([txvars],[TestType]))';

		%  Starting values for search for 95% CI
		if FindCI 
			beta = table2array(lm.Coefficients(2,'Estimate')); 
			se   = table2array(lm.Coefficients(2,'SE')); 
		end
	end

	[ pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0,xvars, model,T0,P,TestType,TestSide,TEST1) ; 

	sprintf('The p-value from Randomization Inference for the hypothesis that tau = %02.2f is %0.2f',tau0, pvalue)

	%  If requested, find confidence interval
	if FindCI 

		%  Find upper boundary of CI. -------------------------------------%
		lb = beta ;
		ub = beta + 10*1.96*se ; 
		middle = (lb + ub) / 2 ; 

		%  Confirm that p-value at upper bound of search region is below significance threshold
		[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,ub,xvars, model,T0,P,TestType,"righttail",TEST1); 
		if p > SignificanceLevel/2
			error('Initial value for upper bound of CI not big enough.')
		end

		QUERIES_UB = NaN(MaxQueries,2); 
		q = 1 ; 
		while q <= MaxQueries % TODO: add step-size constraint.
			fprintf('This is counter number %i \n', q)
			QUERIES_UB(q,1) = middle;
			[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,middle,xvars, model,T0,P,TestType,"righttail",TEST1); 
			QUERIES_UB(q,2) = p;

			%  If reject, move left.  Otherwise, move right.
			if p <= SignificanceLevel/2
				ub = middle;
				middle = (ub + lb ) / 2;  
			else
				lb = middle;
				middle = (ub + lb ) / 2;  
			end

			%  Update counter
			q = q+1;
		end
		CI_UB = min(QUERIES_UB(QUERIES_UB(:,2) > SignificanceLevel/2),1)

		%  Find lower boundary of CI. -------------------------------------%
		ub = beta ;
		lb = beta - 10*1.96*se ; 
		middle = (lb + ub) / 2 ; 

		%  Confirm that p-value at lower bound of search region is below significance threshold
		[p , ~,  ~ ] = ri_estimates(DATA, outcome,txvars,lb,xvars, model,T0,P,TestType,"lefttail",TEST1); 
		if p > SignificanceLevel/2
			error('Initial value for lower bound of CI not low enough.')
		end

		QUERIES_LB = NaN(MaxQueries,2); 
		q = 1 ; 
		while q <= MaxQueries % TODO: add step-size constraint.
			fprintf('This is counter number %i \n', q)
			QUERIES_LB(q,1) = middle;
			[p , ~ , ~ ] = ri_estimates(DATA, outcome,txvars,middle,xvars, model,T0,P,TestType,"lefttail",TEST1); 
			QUERIES_UB(q,2) = p;

			%  If reject, move left.  Otherwise, move right.
			if p <= SignificanceLevel/2
				lb = middle;
				middle = (ub + lb ) / 2;  
			else
				ub = middle;
				middle = (ub + lb ) / 2;  
			end

			%  Update counter
			q = q+1;
		end
		CI_LB = max(QUERIES_UB(QUERIES_UB(:,2) > SignificanceLevel/2),1); 

		%  Return CI as a vector
		CI = [CI_LB , CI_UB];
	end

end


%  Subfunction to conduct RI for a particular value
function [pvalue TEST0 y0 ] = ri_estimates(DATA,outcome,txvars,tau0,xvars, model, T0, P,TestType,TestSide,TEST1)
	%  p-value for hypothesized sharp null.
	%  Run DGP in reverse to get y0
	y0 = table2array(DATA(:,outcome)) - table2array(DATA(:,txvars)) * tau0' ;

	%  Now, loop over feasible randomizations, impose treatment effect, re-estimate, and extract test statistic
	x = table2array(DATA(:,xvars)); % for speed.
	for pp = 1 : P

		%  Impose hypothesized DGP
		ystar = y0 + T0(:,pp,:) * tau0' ; % accmmodates possibliity of multiple treatment variables

		%  Estimate model and collect test statistic
		if strcmp(model,'lm')
			lm = fitlm([permute(T0(:,pp,:),[1 3 2]) , x], ystar);
			testStat = table2array( ...
				lm.Coefficients(2:1+length(txvars) ... % leaving room for constant term
				, [TestType]) ...
				)';
			TEST0(pp,:) = testStat;
		end
	end

	%  get p-value
	if ~strcmp(TestSide, "righttail") , p_left = mean(TEST0 < repmat(TEST1,P,1)) ; end
	if ~strcmp(TestSide, "lefttail") , p_right = mean(TEST0 > repmat(TEST1,P,1)) ; end
	
	if strcmp(TestSide, "righttail"), pvalue = p_right;
	elseif strcmp(TestSide,"lefttail"), pvalue = p_left; 
	else pvalue = min(2*min(p_left,p_right),1) ;
	end

end

